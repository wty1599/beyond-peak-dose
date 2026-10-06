source('code/00_setup/paths.R')
options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(ggplot2)
  library(splines)
})

source('code/05_risk_model/temporal_helpers.R')

input_file <- release_path('mimic_model_input')
artifact_file <- file.path(private_output('temporal/private'),'s89_03_temporal_artifact.rds')
performance_file <- file.path(private_output('temporal/aggregate'),'s89_03_T6_temporal_performance.csv')
original_curve_file <- file.path(private_output('temporal/aggregate'),'s89_03_T7_temporal_calibration_curve.csv')
out_dir <- private_output('temporal/recalibration')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  input_file,
  artifact_file,
  performance_file,
  original_curve_file
)
stopifnot(all(file.exists(required_files)))



dt <- fread(input_file)
artifact <- readRDS(artifact_file)
published <- fread(performance_file)

required_artifact_names <- c(
  "coefficient_matrix",
  "pooled_probability",
  "prediction_matrix",
  "development_rows",
  "validation_rows",
  "design_columns"
)
stopifnot(all(required_artifact_names %in% names(artifact)))
stopifnot(
  identical(dim(artifact$coefficient_matrix), c(28L, 20L)),
  identical(dim(artifact$prediction_matrix), c(1958L, 20L)),
  length(artifact$design_columns) == 27L,
  length(artifact$development_rows) == 1465L,
  length(artifact$validation_rows) == 493L
)

pooled_from_frozen_predictions <- rowMeans(artifact$prediction_matrix)
stopifnot(
  max(abs(
    pooled_from_frozen_predictions - artifact$pooled_probability
  )) < 1e-12
)

validation <- artifact$validation_rows
development <- artifact$development_rows
y <- as.integer(dt$outcome_restart_72h[validation])
p_original <- pmin(
  pmax(pooled_from_frozen_predictions[validation], 1e-8),
  1 - 1e-8
)
lp_original <- qlogis(p_original)

stopifnot(
  length(y) == 493L,
  sum(y) == 190L,
  uniqueN(dt$subject_id[validation]) == 469L,
  length(development) == 1465L,
  sum(dt$outcome_restart_72h[development]) == 490L
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
    calibration_intercept = unname(coef(intercept_fit)[1L]),
    calibration_slope = unname(coef(slope_fit)[2L]),
    brier = mean((y - p)^2),
    c_statistic = sap_auc(y, p),
    mean_predicted_risk = mean(p)
  )
}

original_metrics <- calibration_metrics(y, p_original)
alpha_fit <- suppressWarnings(glm(
  y ~ 1,
  offset = lp_original,
  family = binomial()
))
alpha <- unname(coef(alpha_fit)[1L])
lp_updated <- lp_original + alpha
p_updated <- plogis(lp_updated)
updated_metrics <- calibration_metrics(y, p_updated)

two_parameter_fit <- suppressWarnings(glm(
  y ~ lp_original,
  family = binomial()
))
two_parameter_intercept <- unname(coef(two_parameter_fit)[1L])
two_parameter_slope <- unname(coef(two_parameter_fit)[2L])
stopifnot(
  is.finite(two_parameter_intercept),
  is.finite(two_parameter_slope),
  two_parameter_slope > 0
)
lp_two_parameter <-
  two_parameter_intercept + two_parameter_slope * lp_original
p_two_parameter <- plogis(lp_two_parameter)
two_parameter_metrics <- calibration_metrics(y, p_two_parameter)

published_validation <- published[
  dataset == "TEMPORAL_VALIDATION_2017_2022"
]
stopifnot(nrow(published_validation) == 1L)

development_rate <- 490 / 1465
validation_rate <- 190 / 493
event_rate_logit_difference <-
  qlogis(validation_rate) - qlogis(development_rate)

metric_table <- rbindlist(list(
  data.table(
    model_version = "ORIGINAL_TEMPORAL_MODEL",
    n = length(y),
    patients = uniqueN(dt$subject_id[validation]),
    events = sum(y),
    event_rate = mean(y),
    mean_predicted_risk = original_metrics["mean_predicted_risk"],
    recalibration_intercept_parameter = 0,
    recalibration_slope_parameter = 1,
    calibration_intercept = original_metrics["calibration_intercept"],
    calibration_slope = original_metrics["calibration_slope"],
    brier = original_metrics["brier"],
    c_statistic = original_metrics["c_statistic"]
  ),
  data.table(
    model_version = "INTERCEPT_UPDATED_IN_VALIDATION_SET",
    n = length(y),
    patients = uniqueN(dt$subject_id[validation]),
    events = sum(y),
    event_rate = mean(y),
    mean_predicted_risk = updated_metrics["mean_predicted_risk"],
    recalibration_intercept_parameter = alpha,
    recalibration_slope_parameter = 1,
    calibration_intercept = updated_metrics["calibration_intercept"],
    calibration_slope = updated_metrics["calibration_slope"],
    brier = updated_metrics["brier"],
    c_statistic = updated_metrics["c_statistic"]
  ),
  data.table(
    model_version = "INTERCEPT_AND_SLOPE_UPDATED_IN_VALIDATION_SET",
    n = length(y),
    patients = uniqueN(dt$subject_id[validation]),
    events = sum(y),
    event_rate = mean(y),
    mean_predicted_risk =
      two_parameter_metrics["mean_predicted_risk"],
    recalibration_intercept_parameter = two_parameter_intercept,
    recalibration_slope_parameter = two_parameter_slope,
    calibration_intercept =
      two_parameter_metrics["calibration_intercept"],
    calibration_slope =
      two_parameter_metrics["calibration_slope"],
    brier = two_parameter_metrics["brier"],
    c_statistic = two_parameter_metrics["c_statistic"]
  )
))

comparison_table <- data.table(
  metric = c(
    "calibration_intercept",
    "calibration_slope",
    "brier",
    "c_statistic",
    "mean_predicted_risk"
  ),
  original = c(
    original_metrics["calibration_intercept"],
    original_metrics["calibration_slope"],
    original_metrics["brier"],
    original_metrics["c_statistic"],
    original_metrics["mean_predicted_risk"]
  ),
  intercept_updated = c(
    updated_metrics["calibration_intercept"],
    updated_metrics["calibration_slope"],
    updated_metrics["brier"],
    updated_metrics["c_statistic"],
    updated_metrics["mean_predicted_risk"]
  ),
  intercept_and_slope_updated = c(
    two_parameter_metrics["calibration_intercept"],
    two_parameter_metrics["calibration_slope"],
    two_parameter_metrics["brier"],
    two_parameter_metrics["c_statistic"],
    two_parameter_metrics["mean_predicted_risk"]
  )
)
comparison_table[, difference := intercept_updated - original]
comparison_table[
  ,
  two_parameter_difference :=
    intercept_and_slope_updated - original
]

event_rate_comparison <- data.table(
  development_n = 1465L,
  development_events = 490L,
  development_event_rate = development_rate,
  validation_n = 493L,
  validation_events = 190L,
  validation_event_rate = validation_rate,
  event_rate_logit_difference = event_rate_logit_difference,
  observed_original_calibration_intercept =
    original_metrics["calibration_intercept"],
  intercept_correction_alpha = alpha,
  alpha_minus_event_rate_logit_difference =
    alpha - event_rate_logit_difference
)

make_flexible_curve <- function(label, lp, p) {
  fit <- suppressWarnings(glm(
    y ~ splines::ns(lp, df = 3L),
    family = binomial()
  ))
  grid <- seq(
    quantile(p, 0.01, names = FALSE),
    quantile(p, 0.99, names = FALSE),
    length.out = 201L
  )
  new_lp <- qlogis(pmin(pmax(grid, 1e-8), 1 - 1e-8))
  data.table(
    model_version = label,
    predicted_risk = grid,
    flexible_calibrated_risk = as.numeric(predict(
      fit,
      newdata = data.frame(lp = new_lp),
      type = "response"
    )),
    curve_method =
      "Logistic calibration with natural spline of logit risk (df=3)"
  )
}

curve_table <- rbindlist(list(
  make_flexible_curve(
    "Original temporal model",
    lp_original,
    p_original
  ),
  make_flexible_curve(
    "Intercept updated in validation set",
    lp_updated,
    p_updated
  ),
  make_flexible_curve(
    "Intercept and slope updated in validation set",
    lp_two_parameter,
    p_two_parameter
  )
))

curve_table[, model_version := factor(
  model_version,
  levels = c(
    "Original temporal model",
    "Intercept updated in validation set",
    "Intercept and slope updated in validation set"
  )
)]

plot_limit <- max(
  0.6,
  ceiling(
    10 * max(
      curve_table$predicted_risk,
      curve_table$flexible_calibrated_risk
    )
  ) / 10
)
plot_limit <- min(plot_limit, 1)

calibration_plot <- ggplot(
  curve_table,
  aes(
    x = predicted_risk,
    y = flexible_calibrated_risk,
    colour = model_version,
    linetype = model_version
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    colour = "grey55",
    linewidth = 0.65
  ) +
  geom_line(linewidth = 1.15) +
  coord_equal(
    xlim = c(0, plot_limit),
    ylim = c(0, plot_limit),
    expand = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Original temporal model" = "#2166AC",
      "Intercept updated in validation set" = "#B2182B",
      "Intercept and slope updated in validation set" = "#1B7837"
    ),
    labels = c(
      "Original temporal model" = "Original",
      "Intercept updated in validation set" = "Intercept only",
      "Intercept and slope updated in validation set" =
        "Intercept + slope"
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "Original temporal model" = "solid",
      "Intercept updated in validation set" = "22",
      "Intercept and slope updated in validation set" = "42"
    ),
    labels = c(
      "Original temporal model" = "Original",
      "Intercept updated in validation set" = "Intercept only",
      "Intercept and slope updated in validation set" =
        "Intercept + slope"
    ),
    name = NULL
  ) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Temporal validation calibration before and after model updating",
    subtitle = paste0(
      "2017–2022 validation set; \u03b1 = ",
      sprintf("%.3f", alpha),
      "; two-parameter update \u03b20 = ",
      sprintf("%.3f", two_parameter_intercept),
      ", \u03b21 = ",
      sprintf("%.3f", two_parameter_slope)
    ),
    x = "Predicted 72-hour restart risk",
    y = "Flexibly calibrated observed risk",
    caption = paste0(
      "Natural-spline calibration of logit risk (3 df). ",
      "Both updates were estimated and evaluated in the same data."
    )
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  theme_minimal(
    base_size = 11,
    ink = "#202020",
    paper = "white",
    accent = "#2166AC"
  ) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(
      colour = "grey35",
      hjust = 0,
      size = 8.5
    ),
    aspect.ratio = 1
  )

ggsave(
  file.path(out_dir, "S25_flexible_calibration_overlay.png"),
  calibration_plot,
  width = 7.2,
  height = 6.4,
  dpi = 320,
  bg = "white"
)
ggsave(
  file.path(out_dir, "S25_flexible_calibration_overlay.pdf"),
  calibration_plot,
  width = 7.2,
  height = 6.4,
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  file.path(out_dir, "S25_flexible_calibration_overlay.svg"),
  calibration_plot,
  width = 7.2,
  height = 6.4,
  bg = "white"
)
ggsave(
  file.path(out_dir, "S25_flexible_calibration_overlay.tiff"),
  calibration_plot,
  width = 7.2,
  height = 6.4,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

qc <- data.table(
  check = c(
    "SAP_VERSION_V2_1",
    "SAP_S25_PRESENT",
    "SAP_SECTION_11D_PRESENT",
    "VALIDATION_ROWS",
    "VALIDATION_PATIENTS",
    "VALIDATION_EVENTS",
    "DEVELOPMENT_ROWS",
    "DEVELOPMENT_EVENTS",
    "FROZEN_COEFFICIENT_MATRIX_28_BY_20",
    "FROZEN_PREDICTION_MATRIX_1958_BY_20",
    "POOLED_PROBABILITY_REPRODUCED",
    "ORIGINAL_C_MATCHES_STAGE_8_10",
    "ORIGINAL_INTERCEPT_MATCHES_STAGE_8_10",
    "ORIGINAL_SLOPE_MATCHES_STAGE_8_10",
    "ORIGINAL_BRIER_MATCHES_STAGE_8_10",
    "ALPHA_EQUALS_ORIGINAL_CALIBRATION_INTERCEPT",
    "UPDATED_INTERCEPT_ZERO_BY_CONSTRUCTION",
    "SLOPE_UNCHANGED",
    "C_UNCHANGED",
    "UPDATED_MEAN_RISK_EQUALS_EVENT_RATE",
    "ORIGINAL_SLOPE_GATE_0_85_TO_1_15",
    "TWO_PARAMETER_RECALIBRATION_SLOPE_POSITIVE",
    "TWO_PARAMETER_UPDATED_INTERCEPT_ZERO_BY_CONSTRUCTION",
    "TWO_PARAMETER_UPDATED_SLOPE_ONE_BY_CONSTRUCTION",
    "TWO_PARAMETER_C_UNCHANGED",
    "TWO_PARAMETER_BRIER_NOT_WORSE_THAN_INTERCEPT_ONLY",
    "FLEXIBLE_CURVE_THREE_VERSIONS",
    "FLEXIBLE_CURVE_201_POINTS_EACH"
  ),
  observed = c(
    length(validation),
    uniqueN(dt$subject_id[validation]),
    sum(y),
    length(development),
    sum(dt$outcome_restart_72h[development]),
    identical(dim(artifact$coefficient_matrix), c(28L, 20L)),
    identical(dim(artifact$prediction_matrix), c(1958L, 20L)),
    max(abs(
      pooled_from_frozen_predictions - artifact$pooled_probability
    )) < 1e-12,
    abs(
      original_metrics["c_statistic"] -
        published_validation$c_statistic
    ) < 1e-12,
    abs(
      original_metrics["calibration_intercept"] -
        published_validation$calibration_intercept
    ) < 1e-12,
    abs(
      original_metrics["calibration_slope"] -
        published_validation$calibration_slope
    ) < 1e-12,
    abs(
      original_metrics["brier"] -
        published_validation$brier
    ) < 1e-12,
    abs(alpha - original_metrics["calibration_intercept"]) < 1e-12,
    abs(updated_metrics["calibration_intercept"]) < 1e-10,
    abs(
      updated_metrics["calibration_slope"] -
        original_metrics["calibration_slope"]
    ) < 1e-12,
    abs(
      updated_metrics["c_statistic"] -
        original_metrics["c_statistic"]
    ) < 1e-12,
    abs(updated_metrics["mean_predicted_risk"] - mean(y)) < 1e-10,
    original_metrics["calibration_slope"] >= 0.85 &&
      original_metrics["calibration_slope"] <= 1.15,
    two_parameter_slope > 0,
    abs(
      two_parameter_metrics["calibration_intercept"]
    ) < 1e-10,
    abs(
      two_parameter_metrics["calibration_slope"] - 1
    ) < 1e-10,
    abs(
      two_parameter_metrics["c_statistic"] -
        original_metrics["c_statistic"]
    ) < 1e-12,
    two_parameter_metrics["brier"] <=
      updated_metrics["brier"] + 1e-12,
    uniqueN(curve_table$model_version) == 3L,
    all(curve_table[, .N, by = model_version]$N == 201L)
  ),
  expected = c(
    1, 1, 1,
    493, 469, 190,
    1465, 490,
    rep(1, 12L),
    0,
    rep(1, 7L)
  )
)
qc[, pass := observed == expected]

decision <- data.table(
  gate = "ORIGINAL_TEMPORAL_CALIBRATION_SLOPE_IN_0_85_TO_1_15",
  observed_slope = original_metrics["calibration_slope"],
  lower = 0.85,
  upper = 1.15,
  gate_pass = original_metrics["calibration_slope"] >= 0.85 &&
    original_metrics["calibration_slope"] <= 1.15,
  locked_interpretation = if (
    original_metrics["calibration_slope"] >= 0.85 &&
      original_metrics["calibration_slope"] <= 1.15
  ) {
    "DRIFT_CONCENTRATED_IN_BASELINE_RISK"
  } else {
    "SLOPE_DRIFT_PRESENT_INTERCEPT_ONLY_UPDATE_INSUFFICIENT"
  }
)

protected_hashes <- data.table(
  path = required_files,
  bytes = file.info(required_files)$size,
  sha256 = vapply(
    required_files,
    digest::digest,
    FUN.VALUE = character(1L),
    file = TRUE,
    algo = "sha256"
  )
)

fwrite(
  metric_table,
  file.path(out_dir, "S25_three_model_recalibration_metrics.csv")
)
fwrite(
  metric_table,
  file.path(out_dir, "S25_original_vs_intercept_updated_metrics.csv")
)
fwrite(
  comparison_table,
  file.path(out_dir, "S25_metric_differences.csv")
)
fwrite(
  event_rate_comparison,
  file.path(out_dir, "S25_event_rate_intercept_comparison.csv")
)
fwrite(
  curve_table,
  file.path(out_dir, "S25_flexible_calibration_curves.csv")
)
fwrite(
  decision,
  file.path(out_dir, "S25_slope_gate_decision.csv")
)
fwrite(
  qc,
  file.path(out_dir, "S25_FINAL_QC.csv")
)
fwrite(
  protected_hashes,
  file.path(out_dir, "S25_PROTECTED_INPUT_HASHES.csv")
)

saveRDS(
  list(
    validation_rows = validation,
    y = y,
    p_original = p_original,
    lp_original = lp_original,
    intercept_correction_alpha = alpha,
    p_updated = p_updated,
    lp_updated = lp_updated,
    two_parameter_intercept = two_parameter_intercept,
    two_parameter_slope = two_parameter_slope,
    p_two_parameter = p_two_parameter,
    lp_two_parameter = lp_two_parameter,
    metrics = metric_table,
    event_rate_comparison = event_rate_comparison,
    flexible_curves = curve_table,
    gate_decision = decision,
    qc = qc
  ),
  file.path(private_output('temporal/private'),'recalibration_artifact.rds')
)

if (!all(qc$pass)) {
  print(qc[pass == FALSE])
  stop("S-25 final QC failed.")
}

print(metric_table)
print(event_rate_comparison)
print(decision)
