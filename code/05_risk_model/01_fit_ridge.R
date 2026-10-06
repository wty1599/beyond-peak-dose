source("code/00_setup/paths.R")
options(stringsAsFactors = FALSE, warn = 1)

project_library <- NULL
project_root <- normalizePath(".",winslash="/")
input_file <- release_path("mimic_model_input")
output_directory <- private_output("risk_model")
helper_file <- normalizePath("code/05_risk_model/model_helpers.R",winslash="/")
knot_file <- file.path(output_directory, "sap3_rcs_knots_v1_3.csv")
fold_file <- file.path(output_directory, "sap3_main_patient_fold_map_v1_3.csv")
snapshot_file <- file.path(output_directory, "sap3_input_snapshot_v1_3.csv")
expected_sha256 <- digest::digest(file=input_file,algo="sha256")

suppressPackageStartupMessages({
  library(data.table, lib.loc = project_library)
  library(digest, lib.loc = project_library)
  library(Hmisc, lib.loc = project_library)
  library(mice, lib.loc = project_library)
  library(glmnet, lib.loc = project_library)
})

source(helper_file, local = FALSE)

required_inputs <- c(input_file, helper_file, knot_file, fold_file, snapshot_file)
if (!all(file.exists(required_inputs))) {
  stop("Missing required Stage 3 input(s)")
}
if (!dir.exists(output_directory)) {
  stop("Preflight output directory does not exist")
}

stage3_outputs <- c(
  "sap3_main_model_artifact_v1_3.rds",
  "sap3_main_mids_v1_3.rds",
  "sap3_main_coefficients_v1_3.csv",
  "sap3_main_lambdas_v1_3.csv",
  "sap3_main_cv_deviance_v1_3.csv",
  "sap3_main_cv_fold_log_v1_3.csv",
  "sap3_main_predictions_v1_3.csv",
  "sap3_main_apparent_performance_v1_3.csv",
  "sap3_main_mice_logged_events_v1_3.csv",
  "sap3_main_mice_convergence_v1_3.csv",
  "sap3_main_qc_v1_3.csv"
)
existing_outputs <- file.exists(file.path(output_directory, stage3_outputs))
if (any(existing_outputs)) {
  stop(
    "Refusing to overwrite Stage 3 outputs: ",
    paste(stage3_outputs[existing_outputs], collapse = ", ")
  )
}

actual_sha256 <- digest(file = input_file, algo = "sha256")
if (!identical(actual_sha256, expected_sha256)) {
  stop("Input SHA-256 mismatch")
}

data <- fread(input_file)
sap_assert_main_data(data)
knots <- sap_load_knots(knot_file)
fold_map <- fread(fold_file)
row_fold <- sap_fold_id_for_rows(data, fold_map)
fold_membership_check <- data.table(
  subject_id = data$subject_id,
  fold = row_fold
)[, uniqueN(fold), by = subject_id]
stopifnot(
  length(row_fold) == SAP_EXPECTED_ROWS,
  fold_membership_check[, max(V1)] == 1L
)

cat("Stage 3: nested MICE/CV lambda selection started\n")
cv_result <- sap_select_lambdas_nested_mi(
  dt = data,
  knots = knots,
  fold_map = fold_map,
  seed_base = 202607310L,
  m = SAP_M
)

cat("Stage 3: full-data MICE started\n")
main_mice <- sap_run_mice(
  dt = data,
  seed = 202607320L,
  ignore = rep(FALSE, nrow(data)),
  neutral_outcome_rows = integer(0),
  m = SAP_M,
  maxit = SAP_MAXIT
)

y <- as.numeric(data[[SAP_OUTCOME]])
expected_main_row_key <- sprintf("sap_row_%08d", seq_len(nrow(data)))
for (j in seq_len(SAP_M)) {
  stopifnot(
    nrow(main_mice$completed[[j]]) == SAP_EXPECTED_ROWS,
    identical(
      attr(main_mice$completed[[j]], "sap_row_key"),
      expected_main_row_key
    ),
    sum(y) == SAP_EXPECTED_EVENTS,
    !anyNA(main_mice$completed[[j]][SAP_PREDICTORS])
  )
  design_check <- sap_build_design(main_mice$completed[[j]], knots)
  stopifnot(
    nrow(design_check) == SAP_EXPECTED_ROWS,
    ncol(design_check) == 27L
  )
}

cat("Stage 3: fitting 20 independently tuned ridge models\n")
main_fit <- sap_fit_completed_models(
  completed = main_mice$completed,
  y = y,
  knots = knots,
  lambdas = cv_result$lambda
)

pooled_probability <- main_fit$pooled_probability
if (length(pooled_probability) != SAP_EXPECTED_ROWS ||
    anyNA(pooled_probability) ||
    any(!is.finite(pooled_probability)) ||
    any(pooled_probability <= 0 | pooled_probability >= 1)) {
  stop("Invalid pooled main-model probabilities")
}

apparent_performance <- sap_performance(y, pooled_probability)
apparent_performance[
  ,
  `:=`(
    analysis = "FULL_COHORT_APPARENT",
    landmarks = SAP_EXPECTED_ROWS,
    patients = SAP_EXPECTED_PATIENTS,
    events = SAP_EXPECTED_EVENTS
  )
]
data.table::setcolorder(
  apparent_performance,
  c(
    "analysis", "landmarks", "patients", "events",
    "c_statistic", "calibration_intercept", "calibration_slope"
  )
)

coefficient_table <- sap_pool_coefficient_table(main_fit$coefficients)
lambda_table <- data.table(
  imputation_stream = seq_len(SAP_M),
  selected_lambda = cv_result$lambda,
  selected_lambda_index = cv_result$selected_index,
  minimum_mean_binomial_deviance = vapply(
    seq_len(SAP_M),
    function(j) cv_result$mean_deviance[j, cv_result$selected_index[j]],
    numeric(1)
  )
)
cv_deviance <- data.table(
  imputation_stream = rep(
    seq_len(SAP_M),
    times = length(SAP_LAMBDA_GRID)
  ),
  lambda_index = rep(
    seq_along(SAP_LAMBDA_GRID),
    each = SAP_M
  ),
  lambda = rep(
    SAP_LAMBDA_GRID,
    each = SAP_M
  ),
  mean_binomial_deviance = as.vector(cv_result$mean_deviance)
)

prediction_table <- data.table(
  model_row_id = data$model_row_id,
  subject_id = data$subject_id,
  outcome_restart_72h = y,
  pooled_probability = pooled_probability
)
stopifnot(
  identical(prediction_table$model_row_id, data$model_row_id),
  nrow(prediction_table) == SAP_EXPECTED_ROWS,
  sum(prediction_table$outcome_restart_72h) == SAP_EXPECTED_EVENTS
)

mice_logged <- sap_mice_log_table(main_mice, "full_data_main_mice")
mice_convergence <- tryCatch(
  data.table::as.data.table(mice::convergence(main_mice$mids)),
  error = function(e) data.table(
    convergence_error = conditionMessage(e)
  )
)

qc_table <- data.table(
  qc_number = seq_len(12L),
  check_item = c(
    "Input SHA-256",
    "Main rows preserved",
    "Unique model rows preserved",
    "Patients preserved",
    "Primary events preserved",
    "Imputation streams",
    "MICE iterations",
    "Design columns nonintercept",
    "Lambda selected independently",
    "No lambda at grid boundary",
    "Pooled predictions complete",
    "Lactate absent from main model"
  ),
  observed = c(
    actual_sha256,
    nrow(data),
    uniqueN(data$model_row_id),
    uniqueN(data$subject_id),
    sum(y),
    length(main_mice$completed),
    SAP_MAXIT,
    length(sap_design_column_names),
    length(cv_result$lambda),
    sum(cv_result$selected_index %in% c(1L, length(SAP_LAMBDA_GRID))),
    sum(is.finite(pooled_probability)),
    as.integer(!("predictor_peak_lactate_12h" %in% SAP_PREDICTORS))
  ),
  expected = c(
    expected_sha256,
    SAP_EXPECTED_ROWS,
    SAP_EXPECTED_ROWS,
    SAP_EXPECTED_PATIENTS,
    SAP_EXPECTED_EVENTS,
    SAP_M,
    SAP_MAXIT,
    27,
    SAP_M,
    0,
    SAP_EXPECTED_ROWS,
    1
  )
)
qc_table[, pass := observed == expected]
if (!all(qc_table$pass)) {
  stop("Stage 3 QC failed")
}

model_artifact <- list(
  protocol_version = "original_method_minimal_correction_20260920",
  input_file = input_file,
  input_sha256 = actual_sha256,
  data_ids = data[, .(
    model_row_id,
    subject_id,
    outcome_restart_72h
  )],
  knots = knots,
  design_columns = sap_design_column_names,
  fold_map = fold_map,
  lambda_grid = SAP_LAMBDA_GRID,
  selected_lambda = cv_result$lambda,
  coefficient_matrix = main_fit$coefficients,
  pooled_coefficients = rowMeans(main_fit$coefficients),
  completed_data = main_mice$completed,
  prediction_matrix = main_fit$predictions,
  pooled_probability = pooled_probability,
  apparent_performance = apparent_performance,
  cv_fold_log = cv_result$fold_log,
  main_mice_logged_events = main_mice$logged_events,
  random_seeds = list(
    main_fold = 2026072803L,
    nested_cv_mice_base = 202607310L,
    full_data_mice = 202607320L
  ),
  package_versions = c(
    R = paste(R.version$major, R.version$minor, sep = "."),
    mice = as.character(packageVersion("mice")),
    glmnet = as.character(packageVersion("glmnet")),
    Hmisc = as.character(packageVersion("Hmisc")),
    data.table = as.character(packageVersion("data.table"))
  )
)

saveRDS(
  model_artifact,
  file.path(output_directory, "sap3_main_model_artifact_v1_3.rds"),
  compress = "xz"
)
saveRDS(
  main_mice$mids,
  file.path(output_directory, "sap3_main_mids_v1_3.rds"),
  compress = "xz"
)
fwrite(
  coefficient_table,
  file.path(output_directory, "sap3_main_coefficients_v1_3.csv")
)
fwrite(
  lambda_table,
  file.path(output_directory, "sap3_main_lambdas_v1_3.csv")
)
fwrite(
  cv_deviance,
  file.path(output_directory, "sap3_main_cv_deviance_v1_3.csv")
)
fwrite(
  cv_result$fold_log,
  file.path(output_directory, "sap3_main_cv_fold_log_v1_3.csv")
)
fwrite(
  prediction_table,
  file.path(output_directory, "sap3_main_predictions_v1_3.csv")
)
fwrite(
  apparent_performance,
  file.path(output_directory, "sap3_main_apparent_performance_v1_3.csv")
)
fwrite(
  mice_logged,
  file.path(output_directory, "sap3_main_mice_logged_events_v1_3.csv")
)
fwrite(
  mice_convergence,
  file.path(output_directory, "sap3_main_mice_convergence_v1_3.csv")
)
fwrite(
  qc_table,
  file.path(output_directory, "sap3_main_qc_v1_3.csv")
)

cat("Stage 3 complete; model results written but not interpreted\n")
cat("Rows/patients/events:", nrow(data), uniqueN(data$subject_id), sum(y), "\n")
cat("Independent lambdas:", length(cv_result$lambda), "\n")
cat("All QC passed:", all(qc_table$pass), "\n")
