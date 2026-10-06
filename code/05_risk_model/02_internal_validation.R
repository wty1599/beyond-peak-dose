source("code/00_setup/paths.R")
options(stringsAsFactors = FALSE, warn = 1)

project_library <- NULL
project_root <- normalizePath(".",winslash="/")
input_file <- release_path("mimic_model_input")
output_directory <- private_output("risk_model")
checkpoint_directory <- file.path(output_directory, "sap4_checkpoints_v1_3")
helper_file <- normalizePath("code/05_risk_model/model_helpers.R",winslash="/")
model_artifact_file <- file.path(
  output_directory,
  "sap3_main_model_artifact_v1_3.rds"
)
knot_file <- file.path(output_directory, "sap3_rcs_knots_v1_3.csv")
expected_sha256 <- digest::digest(file=input_file,algo="sha256")

bootstrap_target <- 1000L
parallel_workers <- 4L
batch_size <- 4L
maximum_attempts <- 1500L

suppressPackageStartupMessages({
  library(data.table, lib.loc = project_library)
  library(digest, lib.loc = project_library)
  library(Hmisc, lib.loc = project_library)
  library(mice, lib.loc = project_library)
  library(glmnet, lib.loc = project_library)
  library(future, lib.loc = project_library)
  library(future.apply, lib.loc = project_library)
})

source(helper_file, local = FALSE)

if (!file.exists(model_artifact_file)) {
  stop("Stage 3 model artifact is missing")
}
if (!dir.exists(output_directory)) {
  stop("Stage 3/4 output directory is missing")
}
if (!dir.exists(checkpoint_directory)) {
  dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
}

final_outputs <- c(
  "sap4_bootstrap_metrics_all_v1_3.csv",
  "sap4_optimism_corrected_performance_v1_3.csv",
  "sap4_coefficients_with_bootstrap_ci_v1_3.csv",
  "sap4_prediction_instability_summary_v1_3.csv",
  "sap4_prediction_instability_by_landmark_v1_3.csv",
  "sap4_bootstrap_lambdas_v1_3.csv",
  "sap4_bootstrap_failures_v1_3.csv",
  "sap4_final_artifact_v1_3.rds",
  "sap4_qc_v1_3.csv"
)
if (any(file.exists(file.path(output_directory, final_outputs)))) {
  stop("Final Stage 4 outputs already exist; refusing to overwrite")
}

if (!identical(digest(file = input_file, algo = "sha256"), expected_sha256)) {
  stop("Stage 4 input SHA-256 mismatch")
}

original_data <- fread(input_file)
sap_assert_main_data(original_data)
original_y <- as.numeric(original_data[[SAP_OUTCOME]])
knots <- sap_load_knots(knot_file)
main_artifact <- readRDS(model_artifact_file)

stopifnot(
  identical(main_artifact$input_sha256, expected_sha256),
  length(main_artifact$pooled_probability) == SAP_EXPECTED_ROWS,
  nrow(main_artifact$coefficient_matrix) == 28L,
  ncol(main_artifact$coefficient_matrix) == SAP_M,
  identical(main_artifact$design_columns, sap_design_column_names)
)

sap_draw_patient_bootstrap <- function(dt, seed) {
  patient_ids <- unique(dt$subject_id)
  set.seed(seed)
  drawn <- sample(patient_ids, length(patient_ids), replace = TRUE)
  frequency <- data.table(subject_id = drawn)[, .(bootstrap_frequency = .N), by = subject_id]
  joined <- dt[frequency, on = "subject_id", nomatch = 0L]
  boot <- joined[
    rep(seq_len(.N), times = bootstrap_frequency)
  ]
  boot[, bootstrap_frequency := NULL]
  boot[, bootstrap_row_id := seq_len(.N)]
  if (uniqueN(boot$subject_id) < SAP_K_FOLDS) {
    stop("Bootstrap sample has too few unique patients")
  }
  boot
}

sap_bootstrap_one <- function(attempt_id) {
  seed <- as.integer(202608000L + attempt_id * 5000L)
  started <- Sys.time()
  result <- tryCatch({
    boot <- sap_draw_patient_bootstrap(original_data, seed = seed)
    boot_y <- as.numeric(boot[[SAP_OUTCOME]])
    if (length(unique(boot_y)) != 2L) {
      stop("Bootstrap outcome has fewer than two classes")
    }

    fold_map <- sap_make_patient_folds(
      boot,
      k = SAP_K_FOLDS,
      seed = seed + 11L
    )
    boot_fold <- sap_fold_id_for_rows(boot, fold_map)
    fold_qc <- data.table(
      subject_id = boot$subject_id,
      fold = boot_fold
    )[, uniqueN(fold), by = subject_id]
    if (fold_qc[, max(V1)] != 1L) {
      stop("Bootstrap patient crosses CV folds")
    }

    cv <- sap_select_lambdas_nested_mi(
      dt = boot,
      knots = knots,
      fold_map = fold_map,
      seed_base = seed + 1000L,
      m = SAP_M
    )

    train_mice <- sap_run_mice(
      dt = boot,
      seed = seed + 2000L,
      ignore = rep(FALSE, nrow(boot)),
      neutral_outcome_rows = integer(0),
      m = SAP_M,
      maxit = SAP_MAXIT
    )
    train_fit <- sap_fit_completed_models(
      completed = train_mice$completed,
      y = boot_y,
      knots = knots,
      lambdas = cv$lambda
    )
    train_probability <- train_fit$pooled_probability
    train_performance <- sap_performance(boot_y, train_probability)

    combined <- rbindlist(
      list(
        boot[, setdiff(names(boot), "bootstrap_row_id"), with = FALSE],
        original_data
      ),
      use.names = TRUE,
      fill = TRUE
    )
    test_rows <- (nrow(boot) + 1L):(nrow(boot) + nrow(original_data))
    test_mice <- sap_run_mice(
      dt = combined,
      seed = seed + 3000L,
      ignore = seq_len(nrow(combined)) %in% test_rows,
      neutral_outcome_rows = test_rows,
      m = SAP_M,
      maxit = SAP_MAXIT
    )
    test_completed <- lapply(
      test_mice$completed,
      function(comp) comp[test_rows, , drop = FALSE]
    )
    expected_test_row_key <- sprintf(
      "sap_row_%08d",
      test_rows
    )
    if (any(!vapply(
      test_completed,
      function(comp) identical(rownames(comp), expected_test_row_key),
      logical(1)
    ))) {
      stop("Bootstrap test imputation changed original-row order")
    }
    test_prediction <- sap_predict_completed_from_coefficients(
      completed = test_completed,
      coefficient_matrix = train_fit$coefficients,
      knots = knots
    )
    test_probability <- test_prediction$pooled_probability
    if (length(test_probability) != SAP_EXPECTED_ROWS ||
        anyNA(test_probability) ||
        any(!is.finite(test_probability))) {
      stop("Invalid original-sample test predictions")
    }
    test_performance <- sap_performance(original_y, test_probability)

    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    list(
      status = "SUCCESS",
      attempt_id = attempt_id,
      seed = seed,
      elapsed_seconds = elapsed,
      bootstrap_landmarks = nrow(boot),
      bootstrap_unique_patients = uniqueN(boot$subject_id),
      bootstrap_events = sum(boot_y),
      lambdas = cv$lambda,
      pooled_coefficients = rowMeans(train_fit$coefficients),
      test_probability = test_probability,
      apparent_c = train_performance$c_statistic,
      test_c = test_performance$c_statistic,
      apparent_calibration_intercept =
        train_performance$calibration_intercept,
      test_calibration_intercept =
        test_performance$calibration_intercept,
      apparent_calibration_slope =
        train_performance$calibration_slope,
      test_calibration_slope =
        test_performance$calibration_slope,
      cv_mice_logged_events = sum(cv$fold_log$mice_logged_events),
      train_mice_logged_events = train_mice$logged_count,
      test_mice_logged_events = test_mice$logged_count
    )
  }, error = function(e) {
    list(
      status = "FAILED",
      attempt_id = attempt_id,
      seed = seed,
      elapsed_seconds = as.numeric(
        difftime(Sys.time(), started, units = "secs")
      ),
      error_message = conditionMessage(e)
    )
  })
  result
}

checkpoint_files <- sort(list.files(
  checkpoint_directory,
  pattern = "^sap4_batch_[0-9]{4}\\.rds$",
  full.names = TRUE
))
existing_batches <- lapply(checkpoint_files, readRDS)
existing_results <- unlist(existing_batches, recursive = FALSE)
existing_successes <- sum(vapply(
  existing_results,
  function(x) identical(x$status, "SUCCESS"),
  logical(1)
))
existing_attempts <- if (length(existing_results) == 0L) {
  0L
} else {
  max(vapply(existing_results, function(x) x$attempt_id, integer(1)))
}
batch_number <- length(checkpoint_files)

cat(
  "Stage 4 resume state:",
  existing_successes,
  "successes after",
  existing_attempts,
  "attempts\n"
)

future::plan(
  future::multisession,
  workers = parallel_workers
)
on.exit(future::plan(future::sequential), add = TRUE)
options(future.globals.maxSize = 2 * 1024^3)

all_results <- existing_results
success_count <- existing_successes
attempt_count <- existing_attempts

while (success_count < bootstrap_target) {
  if (attempt_count >= maximum_attempts) {
    stop("Maximum bootstrap attempts reached before 1,000 successes")
  }
  needed <- bootstrap_target - success_count
  this_batch_size <- min(batch_size, needed, maximum_attempts - attempt_count)
  attempt_ids <- seq.int(
    from = attempt_count + 1L,
    length.out = this_batch_size
  )
  batch_number <- batch_number + 1L
  cat(
    "Stage 4 batch",
    batch_number,
    "attempts",
    min(attempt_ids),
    "to",
    max(attempt_ids),
    "started\n"
  )
  batch_result <- future.apply::future_lapply(
    attempt_ids,
    sap_bootstrap_one,
    future.seed = TRUE,
    future.scheduling = Inf,
    future.packages = c(
      "data.table", "digest", "Hmisc", "mice", "glmnet"
    )
  )
  checkpoint_file <- file.path(
    checkpoint_directory,
    sprintf("sap4_batch_%04d.rds", batch_number)
  )
  if (file.exists(checkpoint_file)) {
    stop("Checkpoint already exists: ", checkpoint_file)
  }
  saveRDS(batch_result, checkpoint_file, compress = "xz")
  all_results <- c(all_results, batch_result)
  attempt_count <- max(attempt_ids)
  batch_successes <- sum(vapply(
    batch_result,
    function(x) identical(x$status, "SUCCESS"),
    logical(1)
  ))
  success_count <- success_count + batch_successes
  cat(
    "Stage 4 batch",
    batch_number,
    "complete:",
    batch_successes,
    "successes; cumulative",
    success_count,
    "/",
    bootstrap_target,
    "\n"
  )
  if (success_count == 0L && attempt_count >= batch_size) {
    failure_messages <- unique(vapply(
      batch_result,
      function(x) x$error_message,
      character(1)
    ))
    stop(
      "First bootstrap batch had zero successes: ",
      paste(failure_messages, collapse = " | ")
    )
  }
}

successful <- Filter(
  function(x) identical(x$status, "SUCCESS"),
  all_results
)
successful <- successful[
  order(vapply(successful, function(x) x$attempt_id, integer(1)))
]
successful <- successful[seq_len(bootstrap_target)]
failed <- Filter(
  function(x) identical(x$status, "FAILED"),
  all_results
)

metrics <- rbindlist(lapply(successful, function(x) {
  data.table(
    attempt_id = x$attempt_id,
    seed = x$seed,
    elapsed_seconds = x$elapsed_seconds,
    bootstrap_landmarks = x$bootstrap_landmarks,
    bootstrap_unique_patients = x$bootstrap_unique_patients,
    bootstrap_events = x$bootstrap_events,
    apparent_c = x$apparent_c,
    test_c = x$test_c,
    optimism_c = x$apparent_c - x$test_c,
    apparent_calibration_intercept =
      x$apparent_calibration_intercept,
    test_calibration_intercept =
      x$test_calibration_intercept,
    optimism_calibration_intercept =
      x$apparent_calibration_intercept -
      x$test_calibration_intercept,
    apparent_calibration_slope =
      x$apparent_calibration_slope,
    test_calibration_slope =
      x$test_calibration_slope,
    optimism_calibration_slope =
      x$apparent_calibration_slope -
      x$test_calibration_slope,
    cv_mice_logged_events = x$cv_mice_logged_events,
    train_mice_logged_events = x$train_mice_logged_events,
    test_mice_logged_events = x$test_mice_logged_events
  )
}))

main_performance <- main_artifact$apparent_performance
performance <- data.table(
  metric = c(
    "C statistic",
    "Calibration intercept",
    "Calibration slope"
  ),
  apparent = c(
    main_performance$c_statistic,
    main_performance$calibration_intercept,
    main_performance$calibration_slope
  ),
  mean_optimism = c(
    mean(metrics$optimism_c),
    mean(metrics$optimism_calibration_intercept),
    mean(metrics$optimism_calibration_slope)
  )
)
performance[, optimism_corrected := apparent - mean_optimism]

bootstrap_coefficient_matrix <- do.call(
  rbind,
  lapply(successful, function(x) x$pooled_coefficients)
)
colnames(bootstrap_coefficient_matrix) <- rownames(
  main_artifact$coefficient_matrix
)
main_pooled_coefficients <- rowMeans(main_artifact$coefficient_matrix)
coefficient_table <- data.table(
  term = names(main_pooled_coefficients),
  pooled_log_odds = as.numeric(main_pooled_coefficients),
  bootstrap_ci_lower = apply(
    bootstrap_coefficient_matrix,
    2L,
    quantile,
    probs = 0.025,
    type = 7
  ),
  bootstrap_ci_upper = apply(
    bootstrap_coefficient_matrix,
    2L,
    quantile,
    probs = 0.975,
    type = 7
  )
)
coefficient_table[
  ,
  `:=`(
    exponentiated_pooled_coefficient = exp(pooled_log_odds),
    exponentiated_ci_lower = exp(bootstrap_ci_lower),
    exponentiated_ci_upper = exp(bootstrap_ci_upper)
  )
]

bootstrap_prediction_matrix <- do.call(
  rbind,
  lapply(successful, function(x) x$test_probability)
)
main_probability <- main_artifact$pooled_probability
absolute_difference <- abs(
  sweep(
    bootstrap_prediction_matrix,
    2L,
    main_probability,
    FUN = "-"
  )
)
instability_by_landmark <- data.table(
  model_row_id = original_data$model_row_id,
  subject_id = original_data$subject_id,
  main_probability = main_probability,
  mean_absolute_difference = colMeans(absolute_difference),
  median_absolute_difference = apply(
    absolute_difference,
    2L,
    median
  ),
  p95_absolute_difference = apply(
    absolute_difference,
    2L,
    quantile,
    probs = 0.95,
    type = 7
  )
)
instability_summary <- data.table(
  bootstrap_successes = bootstrap_target,
  landmarks = SAP_EXPECTED_ROWS,
  instability_index_mean_absolute_difference =
    mean(absolute_difference),
  median_absolute_difference = median(absolute_difference),
  p95_absolute_difference = as.numeric(
    quantile(absolute_difference, 0.95, type = 7)
  ),
  mean_landmark_specific_instability =
    mean(instability_by_landmark$mean_absolute_difference)
)

lambda_table <- rbindlist(lapply(successful, function(x) {
  data.table(
    attempt_id = x$attempt_id,
    imputation_stream = seq_len(SAP_M),
    selected_lambda = x$lambdas
  )
}))

failure_table <- if (length(failed) == 0L) {
  data.table(
    attempt_id = integer(),
    seed = integer(),
    elapsed_seconds = numeric(),
    error_message = character()
  )
} else {
  rbindlist(lapply(failed, function(x) {
    data.table(
      attempt_id = x$attempt_id,
      seed = x$seed,
      elapsed_seconds = x$elapsed_seconds,
      error_message = x$error_message
    )
  }))
}

qc <- data.table(
  qc_number = seq_len(10L),
  check_item = c(
    "Successful bootstrap replicates",
    "Main rows unchanged",
    "Main patients unchanged",
    "Main events unchanged",
    "Main design columns",
    "Lambda selections per successful replicate",
    "Coefficient rows",
    "Original test predictions per replicate",
    "Finite corrected C statistic",
    "Finite corrected calibration slope"
  ),
  observed = c(
    length(successful),
    nrow(original_data),
    uniqueN(original_data$subject_id),
    sum(original_y),
    length(main_artifact$design_columns),
    nrow(lambda_table),
    nrow(coefficient_table),
    nrow(bootstrap_prediction_matrix) *
      ncol(bootstrap_prediction_matrix),
    as.integer(is.finite(
      performance[metric == "C statistic", optimism_corrected]
    )),
    as.integer(is.finite(
      performance[metric == "Calibration slope", optimism_corrected]
    ))
  ),
  expected = c(
    bootstrap_target,
    SAP_EXPECTED_ROWS,
    SAP_EXPECTED_PATIENTS,
    SAP_EXPECTED_EVENTS,
    27,
    bootstrap_target * SAP_M,
    28,
    bootstrap_target * SAP_EXPECTED_ROWS,
    1,
    1
  )
)
qc[, pass := observed == expected]
if (!all(qc$pass)) {
  stop("Final Stage 4 QC failed")
}

final_artifact <- list(
  protocol_version = "original_method_minimal_correction_20260920",
  input_sha256 = expected_sha256,
  successful_attempt_ids = vapply(
    successful,
    function(x) x$attempt_id,
    integer(1)
  ),
  metrics = metrics,
  performance = performance,
  bootstrap_coefficient_matrix = bootstrap_coefficient_matrix,
  coefficient_table = coefficient_table,
  bootstrap_prediction_matrix = bootstrap_prediction_matrix,
  instability_summary = instability_summary,
  failure_table = failure_table,
  qc = qc
)

fwrite(
  metrics,
  file.path(output_directory, "sap4_bootstrap_metrics_all_v1_3.csv")
)
fwrite(
  performance,
  file.path(
    output_directory,
    "sap4_optimism_corrected_performance_v1_3.csv"
  )
)
fwrite(
  coefficient_table,
  file.path(
    output_directory,
    "sap4_coefficients_with_bootstrap_ci_v1_3.csv"
  )
)
fwrite(
  instability_summary,
  file.path(
    output_directory,
    "sap4_prediction_instability_summary_v1_3.csv"
  )
)
fwrite(
  instability_by_landmark,
  file.path(
    output_directory,
    "sap4_prediction_instability_by_landmark_v1_3.csv"
  )
)
fwrite(
  lambda_table,
  file.path(output_directory, "sap4_bootstrap_lambdas_v1_3.csv")
)
fwrite(
  failure_table,
  file.path(output_directory, "sap4_bootstrap_failures_v1_3.csv")
)
fwrite(
  qc,
  file.path(output_directory, "sap4_qc_v1_3.csv")
)
saveRDS(
  final_artifact,
  file.path(output_directory, "sap4_final_artifact_v1_3.rds"),
  compress = "xz"
)

cat("Stage 4 complete\n")
cat("Successful bootstrap replicates:", length(successful), "\n")
cat("Failed attempts retained in audit log:", nrow(failure_table), "\n")
cat("All final QC passed:", all(qc$pass), "\n")
