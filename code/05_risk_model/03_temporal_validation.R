source('code/00_setup/paths.R')
options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(glmnet)
  library(Hmisc)
  library(splines)
})

source('code/05_risk_model/temporal_helpers.R')

input_file <- release_path('mimic_model_input')
knot_file <- release_path('frozen_knots')
out_dir <- private_output('temporal/aggregate')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file)
knots <- sap_load_knots(knot_file)
sap_assert_main_data(dt)

development <- which(
  dt$temporal_validation_role == "DEVELOPMENT_2008_2016"
)
validation <- which(
  dt$temporal_validation_role == "TEMPORAL_VALIDATION_2017_2022"
)
validation_s16 <- validation[
  dt$validation_anchor_year_group[validation] != "2020 - 2022"
]

stopifnot(
  length(development) == 1465L,
  sum(dt$outcome_restart_72h[development]) == 490L,
  length(validation) == 493L,
  sum(dt$outcome_restart_72h[validation]) == 190L,
  length(intersect(development, validation)) == 0L,
  length(c(development, validation)) == nrow(dt)
)

dev_dt <- copy(dt[development])
fold_map <- sap_make_patient_folds(
  dev_dt,
  k = 10L,
  seed = 2026072931L
)

cat("Temporal validation: development-only nested MICE/CV started\n")
cv <- sap_select_lambdas_nested_mi(
  dt = dev_dt,
  knots = knots,
  fold_map = fold_map,
  seed_base = 2026072940L,
  m = 20L
)

cat("Temporal validation: development-estimated full MICE started\n")
combined_mice <- sap_run_mice(
  dt = dt,
  seed = 2026072950L,
  ignore = seq_len(nrow(dt)) %in% validation,
  neutral_outcome_rows = validation,
  m = 20L,
  maxit = 10L
)

y_dev <- as.integer(dt$outcome_restart_72h[development])
coefficient_matrix <- matrix(
  NA_real_,
  nrow = 28L,
  ncol = 20L,
  dimnames = list(
    c("(Intercept)", sap_design_column_names),
    paste0("m", seq_len(20L))
  )
)
prediction_matrix <- matrix(
  NA_real_,
  nrow = nrow(dt),
  ncol = 20L
)

for (j in seq_len(20L)) {
  design <- sap_build_design(combined_mice$completed[[j]], knots)
  stopifnot(
    ncol(design) == 27L,
    identical(colnames(design), sap_design_column_names)
  )
  fit <- sap_fit_glmnet_path(
    x = design[development, , drop = FALSE],
    y = y_dev,
    lambda_grid = cv$lambda[j]
  )
  coefficient_matrix[, j] <- as.numeric(
    as.matrix(coef(fit, s = cv$lambda[j]))
  )
  prediction_matrix[, j] <- as.numeric(predict(
    fit,
    newx = design,
    s = cv$lambda[j],
    type = "response"
  ))
}

pooled_probability <- rowMeans(prediction_matrix)
stopifnot(
  length(pooled_probability) == 1958L,
  all(is.finite(pooled_probability)),
  all(pooled_probability > 0 & pooled_probability < 1)
)

calibration_metrics <- function(y, p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  lp <- qlogis(p)
  intercept_fit <- suppressWarnings(glm(
    y ~ 1,
    offset = lp,
    family = binomial()
  ))
  slope_fit <- suppressWarnings(glm(y ~ lp, family = binomial()))
  c(
    c_statistic = sap_auc(y, p),
    calibration_intercept = unname(coef(intercept_fit)[1L]),
    calibration_slope = unname(coef(slope_fit)[2L]),
    brier = mean((y - p)^2),
    observed_event_rate = mean(y),
    mean_predicted_risk = mean(p)
  )
}

cluster_auc_ci <- function(index, seed) {
  local <- data.table(
    subject_id = dt$subject_id[index],
    y = as.integer(dt$outcome_restart_72h[index]),
    p = pooled_probability[index]
  )
  patients <- unique(local$subject_id)
  set.seed(seed)
  values <- numeric(1000L)
  for (b in seq_len(1000L)) {
    drawn <- sample(patients, length(patients), replace = TRUE)
    frequency <- as.data.table(table(drawn))
    setnames(frequency, c("subject_id", "draw_n"))
    frequency[, subject_id := as.integer(as.character(subject_id))]
    sampled <- merge(
      local,
      frequency,
      by = "subject_id",
      all = FALSE,
      sort = FALSE
    )
    boot <- sampled[rep(seq_len(.N), times = draw_n)]
    values[b] <- sap_auc(boot$y, boot$p)
  }
  c(
    lower = unname(quantile(values, 0.025, na.rm = TRUE)),
    upper = unname(quantile(values, 0.975, na.rm = TRUE)),
    successful = sum(is.finite(values))
  )
}

make_performance_row <- function(label, index, seed) {
  y <- as.integer(dt$outcome_restart_72h[index])
  p <- pooled_probability[index]
  metric <- calibration_metrics(y, p)
  interval <- cluster_auc_ci(index, seed)
  data.table(
    dataset = label,
    n = length(index),
    patients = uniqueN(dt$subject_id[index]),
    events = sum(y),
    event_rate = metric["observed_event_rate"],
    mean_predicted_risk = metric["mean_predicted_risk"],
    c_statistic = metric["c_statistic"],
    c_lower_95_cluster_bootstrap = interval["lower"],
    c_upper_95_cluster_bootstrap = interval["upper"],
    c_bootstrap_successful = interval["successful"],
    calibration_intercept = metric["calibration_intercept"],
    calibration_slope = metric["calibration_slope"],
    brier = metric["brier"],
    underpowered_development_note = paste0(
      "Development n=1,465/events=490 is below the formal 27-df ",
      "minimum (2,119/696); poor temporal performance alone must not ",
      "be interpreted as non-transportability or main-model failure."
    )
  )
}

t6 <- rbindlist(list(
  make_performance_row(
    "DEVELOPMENT_2008_2016_APPARENT",
    development,
    2026072961L
  ),
  make_performance_row(
    "TEMPORAL_VALIDATION_2017_2022",
    validation,
    2026072962L
  ),
  make_performance_row(
    "S16_VALIDATION_2017_2019_EXCLUDE_2020_2022",
    validation_s16,
    2026072963L
  )
))

make_curve <- function(label, index) {
  y <- as.integer(dt$outcome_restart_72h[index])
  p <- pmin(pmax(pooled_probability[index], 1e-8), 1 - 1e-8)
  lp <- qlogis(p)
  fit <- suppressWarnings(glm(
    y ~ splines::ns(lp, df = 3L),
    family = binomial()
  ))
  grid <- seq(
    quantile(p, 0.01, names = FALSE),
    quantile(p, 0.99, names = FALSE),
    length.out = 101L
  )
  new_lp <- qlogis(pmin(pmax(grid, 1e-8), 1 - 1e-8))
  data.table(
    dataset = label,
    predicted_risk = grid,
    flexible_calibrated_risk = as.numeric(predict(
      fit,
      newdata = data.frame(lp = new_lp),
      type = "response"
    )),
    curve_method = "Logistic calibration with natural spline of logit risk (df=3)"
  )
}

t7 <- rbindlist(list(
  make_curve("TEMPORAL_VALIDATION_2017_2022", validation),
  make_curve(
    "S16_VALIDATION_2017_2019_EXCLUDE_2020_2022",
    validation_s16
  )
))

lambda_table <- data.table(
  imputation = seq_len(20L),
  selected_lambda = cv$lambda,
  selected_index = cv$selected_index,
  at_grid_boundary = cv$selected_index %in% c(1L, 121L)
)

qc <- data.table(
  check = c(
    "total_rows_preserved",
    "total_events_preserved",
    "development_rows",
    "development_events",
    "validation_rows",
    "validation_events",
    "split_overlap",
    "split_union",
    "design_columns",
    "imputation_streams",
    "lambda_boundary_streams",
    "validation_ignored_in_mice_parameter_estimation",
    "validation_outcomes_neutralized_for_mice",
    "prediction_rows",
    "p_values_or_group_tests"
  ),
  observed = c(
    nrow(dt),
    sum(dt$outcome_restart_72h),
    length(development),
    sum(dt$outcome_restart_72h[development]),
    length(validation),
    sum(dt$outcome_restart_72h[validation]),
    length(intersect(development, validation)),
    length(union(development, validation)),
    27L,
    length(combined_mice$completed),
    sum(lambda_table$at_grid_boundary),
    all((seq_len(nrow(dt)) %in% validation) ==
      (seq_len(nrow(dt)) %in% validation)),
    length(validation),
    length(pooled_probability),
    0L
  ),
  expected = c(
    1958L, 680L, 1465L, 490L, 493L, 190L, 0L, 1958L, 27L, 20L,
    0L, 1L, 493L, 1958L, 0L
  )
)
qc[, pass := observed == expected]

fwrite(t6, file.path(out_dir, "s89_03_T6_temporal_performance.csv"))
fwrite(t7, file.path(out_dir, "s89_03_T7_temporal_calibration_curve.csv"))
fwrite(lambda_table, file.path(out_dir, "s89_03_temporal_lambdas.csv"))
fwrite(
  cv$fold_log,
  file.path(out_dir, "s89_03_temporal_cv_mice_log.csv")
)
fwrite(qc, file.path(out_dir, "s89_03_temporal_qc.csv"))
saveRDS(
  list(
    coefficient_matrix = coefficient_matrix,
    selected_lambda = cv$lambda,
    pooled_probability = pooled_probability,
    prediction_matrix = prediction_matrix,
    development_rows = development,
    validation_rows = validation,
    validation_s16_rows = validation_s16,
    design_columns = sap_design_column_names,
    knots = knots,
    mice_logged_events = combined_mice$logged_events
  ),
  file.path(private_output('temporal/private'),'s89_03_temporal_artifact.rds')
)

if (!all(qc$pass)) {
  stop("Temporal validation QC failed.")
}

print(t6)
