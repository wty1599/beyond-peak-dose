options(stringsAsFactors = FALSE, warn = 1)

SAP_OUTCOME <- "outcome_restart_72h"
SAP_ID_NAMES <- c(
  "model_row_id", "d2v2_event_id", "subject_id", "hadm_id", "stay_id", "t0"
)
SAP_PREDICTORS <- c(
  "predictor_age",
  "predictor_female",
  "predictor_episode_peak_nee",
  "predictor_last_positive_nee",
  "predictor_positive_nee_hours",
  "predictor_multivaso_peak",
  "predictor_sofa_24h",
  "predictor_mechanical_ventilation_t0",
  "predictor_rrt_t0",
  "predictor_mean_map_4h",
  "predictor_map_source_invasive",
  "predictor_mean_heart_rate_4h",
  "predictor_urine_output_24h_ml"
)
SAP_SPLINE_NAMES <- c(
  "predictor_age",
  "predictor_episode_peak_nee",
  "predictor_last_positive_nee",
  "predictor_positive_nee_hours",
  "predictor_sofa_24h",
  "predictor_mean_map_4h",
  "predictor_urine_output_24h_ml"
)
SAP_BINARY_NAMES <- c(
  "predictor_female",
  "predictor_multivaso_peak",
  "predictor_mechanical_ventilation_t0",
  "predictor_rrt_t0",
  "predictor_map_source_invasive"
)
SAP_LINEAR_NAMES <- "predictor_mean_heart_rate_4h"
SAP_IMPUTE_PMM <- c(
  "predictor_mean_map_4h",
  "predictor_mean_heart_rate_4h",
  "predictor_urine_output_24h_ml"
)
SAP_IMPUTE_LOGREG <- "predictor_map_source_invasive"
SAP_M <- 20L
SAP_MAXIT <- 10L
SAP_DONORS <- 5L
SAP_K_FOLDS <- 10L
SAP_LAMBDA_GRID <- 10^seq(2, -6, length.out = 121)
SAP_EXPECTED_ROWS <- 1958L
SAP_EXPECTED_PATIENTS <- 1856L
SAP_EXPECTED_EVENTS <- 680L

sap_design_column_names <- c(
  unlist(lapply(
    SAP_SPLINE_NAMES,
    function(v) paste0(v, "_rcs", 1:3)
  )),
  SAP_BINARY_NAMES,
  SAP_LINEAR_NAMES
)

sap_assert_main_data <- function(dt) {
  stopifnot(
    nrow(dt) == SAP_EXPECTED_ROWS,
    data.table::uniqueN(dt$model_row_id) == SAP_EXPECTED_ROWS,
    data.table::uniqueN(dt$subject_id) == SAP_EXPECTED_PATIENTS,
    sum(dt[[SAP_OUTCOME]]) == SAP_EXPECTED_EVENTS,
    all(dt[[SAP_OUTCOME]] %in% c(0L, 1L)),
    all(c(SAP_ID_NAMES, SAP_OUTCOME, SAP_PREDICTORS) %in% names(dt)),
    !("predictor_peak_lactate_12h" %in% names(dt))
  )
  invisible(TRUE)
}

sap_load_knots <- function(knot_file) {
  knot_table <- data.table::fread(knot_file)
  required <- c("predictor", "knot_index", "probability", "knot_value")
  stopifnot(all(required %in% names(knot_table)))
  knots <- split(knot_table$knot_value, knot_table$predictor)
  stopifnot(
    setequal(names(knots), SAP_SPLINE_NAMES),
    all(vapply(knots, length, integer(1)) == 4L),
    all(vapply(knots, function(x) length(unique(x)), integer(1)) == 4L)
  )
  knots[SAP_SPLINE_NAMES]
}

sap_build_design <- function(frame, knots) {
  pieces <- list()
  for (v in SAP_SPLINE_NAMES) {
    x <- as.numeric(frame[[v]])
    if (anyNA(x) || any(!is.finite(x))) {
      stop("Non-finite value in completed spline predictor: ", v)
    }
    basis <- Hmisc::rcspline.eval(
      x,
      knots = knots[[v]],
      inclx = TRUE,
      norm = 2
    )
    if (ncol(basis) != 3L) {
      stop("RCS basis for ", v, " produced ", ncol(basis), " columns")
    }
    colnames(basis) <- paste0(v, "_rcs", 1:3)
    pieces[[v]] <- basis
  }
  for (v in SAP_BINARY_NAMES) {
    x <- frame[[v]]
    if (is.factor(x)) {
      x <- as.character(x)
    }
    x <- as.numeric(x)
    if (anyNA(x) || any(!(x %in% c(0, 1)))) {
      stop("Invalid completed binary predictor: ", v)
    }
    pieces[[v]] <- matrix(x, ncol = 1L, dimnames = list(NULL, v))
  }
  for (v in SAP_LINEAR_NAMES) {
    x <- as.numeric(frame[[v]])
    if (anyNA(x) || any(!is.finite(x))) {
      stop("Non-finite completed linear predictor: ", v)
    }
    pieces[[v]] <- matrix(x, ncol = 1L, dimnames = list(NULL, v))
  }
  design <- do.call(cbind, pieces)
  design <- design[, sap_design_column_names, drop = FALSE]
  stopifnot(
    nrow(design) == nrow(frame),
    ncol(design) == 27L,
    identical(colnames(design), sap_design_column_names),
    !anyNA(design),
    all(is.finite(design))
  )
  storage.mode(design) <- "double"
  design
}

sap_make_patient_folds <- function(dt, k = SAP_K_FOLDS, seed = 1L) {
  patient <- dt[
    ,
    .(
      n_rows = .N,
      n_events = sum(get(SAP_OUTCOME)),
      any_event = as.integer(any(get(SAP_OUTCOME) == 1L))
    ),
    by = subject_id
  ]
  if (nrow(patient) < k) {
    stop("Fewer patients than CV folds")
  }
  set.seed(seed)
  patient[, random_tie := runif(.N)]
  data.table::setorder(patient, -n_events, -n_rows, -any_event, random_tie)
  fold_events <- numeric(k)
  fold_rows <- numeric(k)
  fold_patients <- integer(k)
  patient[, fold := NA_integer_]
  for (i in seq_len(nrow(patient))) {
    score <- data.table::data.table(
      fold = seq_len(k),
      events = fold_events,
      rows = fold_rows,
      patients = fold_patients
    )
    data.table::setorder(score, events, rows, patients, fold)
    selected <- score$fold[1L]
    patient$fold[i] <- selected
    fold_events[selected] <- fold_events[selected] + patient$n_events[i]
    fold_rows[selected] <- fold_rows[selected] + patient$n_rows[i]
    fold_patients[selected] <- fold_patients[selected] + 1L
  }
  patient[, random_tie := NULL]
  patient[]
}

sap_fold_id_for_rows <- function(dt, fold_map) {
  lookup <- fold_map$fold[match(dt$subject_id, fold_map$subject_id)]
  if (anyNA(lookup)) {
    stop("Patient fold map does not cover every row")
  }
  if (data.table::data.table(
    subject_id = dt$subject_id,
    fold = lookup
  )[, data.table::uniqueN(fold), by = subject_id][, max(V1)] != 1L) {
    stop("A patient crosses CV folds")
  }
  as.integer(lookup)
}

sap_prepare_mice_frame <- function(dt) {
  frame <- as.data.frame(dt[, c(SAP_OUTCOME, SAP_PREDICTORS), with = FALSE])
  rownames(frame) <- sprintf("sap_row_%08d", seq_len(nrow(dt)))
  frame[[SAP_OUTCOME]] <- as.numeric(frame[[SAP_OUTCOME]])
  frame[[SAP_IMPUTE_LOGREG]] <- factor(
    frame[[SAP_IMPUTE_LOGREG]],
    levels = c(0, 1)
  )
  frame
}

sap_run_mice <- function(
  dt,
  seed,
  ignore = rep(FALSE, nrow(dt)),
  neutral_outcome_rows = integer(0),
  m = SAP_M,
  maxit = SAP_MAXIT
) {
  if (length(ignore) != nrow(dt)) {
    stop("MICE ignore vector has wrong length")
  }
  frame <- sap_prepare_mice_frame(dt)
  expected_row_key <- rownames(frame)
  training_rows <- !ignore
  if (!any(training_rows)) {
    stop("MICE has no training rows")
  }
  if (length(neutral_outcome_rows) > 0L) {
    training_prevalence <- mean(frame[[SAP_OUTCOME]][training_rows])
    frame[[SAP_OUTCOME]][neutral_outcome_rows] <- training_prevalence
  }

  initial <- mice::mice(
    frame,
    m = 1L,
    maxit = 0L,
    printFlag = FALSE
  )
  method <- rep("", ncol(frame))
  names(method) <- names(frame)
  for (v in SAP_IMPUTE_PMM) {
    if (anyNA(frame[[v]])) {
      method[v] <- "pmm"
    }
  }
  if (anyNA(frame[[SAP_IMPUTE_LOGREG]])) {
    method[SAP_IMPUTE_LOGREG] <- "logreg"
  }

  predictor_matrix <- matrix(
    0L,
    nrow = ncol(frame),
    ncol = ncol(frame),
    dimnames = list(names(frame), names(frame))
  )
  targets <- names(method)[method != ""]
  for (target in targets) {
    predictor_matrix[target, setdiff(names(frame), target)] <- 1L
  }
  predictor_matrix[SAP_OUTCOME, ] <- 0L

  mids <- mice::mice(
    frame,
    m = as.integer(m),
    maxit = as.integer(maxit),
    method = method,
    predictorMatrix = predictor_matrix,
    ignore = as.logical(ignore),
    donors = SAP_DONORS,
    seed = as.integer(seed),
    printFlag = FALSE
  )

  completed <- lapply(seq_len(m), function(i) {
    comp <- mice::complete(mids, action = i)
    if (!identical(rownames(comp), expected_row_key)) {
      stop("MICE changed row order")
    }
    comp[[SAP_IMPUTE_LOGREG]] <- as.numeric(
      as.character(comp[[SAP_IMPUTE_LOGREG]])
    )
    if (nrow(comp) != nrow(dt)) {
      stop("MICE changed row count")
    }
    if (anyNA(comp[SAP_PREDICTORS])) {
      missing_after <- names(which(vapply(
        comp[SAP_PREDICTORS],
        anyNA,
        logical(1)
      )))
      stop(
        "MICE left missing predictors: ",
        paste(missing_after, collapse = ", ")
      )
    }
    attr(comp, "sap_row_key") <- expected_row_key
    comp
  })

  logged <- mids$loggedEvents
  logged_count <- if (is.null(logged)) 0L else nrow(logged)
  list(
    completed = completed,
    mids = mids,
    logged_count = logged_count,
    logged_events = logged,
    method = method,
    predictor_matrix = predictor_matrix
    ,
    row_key = expected_row_key
  )
}

sap_binomial_deviance <- function(y, p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  -2 * colSums(
    y * log(p) + (1 - y) * log(1 - p)
  )
}

sap_fit_glmnet_path <- function(x, y, lambda_grid = SAP_LAMBDA_GRID) {
  suppressWarnings(
    glmnet::glmnet(
      x = x,
      y = y,
      family = "binomial",
      alpha = 0,
      lambda = lambda_grid,
      standardize = TRUE,
      intercept = TRUE,
      thresh = 1e-8,
      maxit = 100000
    )
  )
}

sap_select_lambdas_nested_mi <- function(
  dt,
  knots,
  fold_map,
  seed_base,
  m = SAP_M
) {
  row_fold <- sap_fold_id_for_rows(dt, fold_map)
  y <- as.numeric(dt[[SAP_OUTCOME]])
  loss_sum <- matrix(
    0,
    nrow = m,
    ncol = length(SAP_LAMBDA_GRID)
  )
  validation_rows <- integer(m)
  mice_log <- vector("list", SAP_K_FOLDS)

  for (fold in seq_len(SAP_K_FOLDS)) {
    validation <- which(row_fold == fold)
    training <- which(row_fold != fold)
    if (length(validation) == 0L || length(unique(y[training])) != 2L) {
      stop("Invalid CV fold ", fold)
    }

    fold_mice <- sap_run_mice(
      dt = dt,
      seed = as.integer(seed_base + fold * 101L),
      ignore = seq_len(nrow(dt)) %in% validation,
      neutral_outcome_rows = validation,
      m = m,
      maxit = SAP_MAXIT
    )
    mice_log[[fold]] <- data.table::data.table(
      fold = fold,
      training_rows = length(training),
      validation_rows = length(validation),
      training_events = sum(y[training]),
      validation_events = sum(y[validation]),
      mice_logged_events = fold_mice$logged_count
    )

    for (j in seq_len(m)) {
      comp <- fold_mice$completed[[j]]
      x <- sap_build_design(comp, knots)
      fit <- sap_fit_glmnet_path(
        x = x[training, , drop = FALSE],
        y = y[training],
        lambda_grid = SAP_LAMBDA_GRID
      )
      prediction <- predict(
        fit,
        newx = x[validation, , drop = FALSE],
        s = SAP_LAMBDA_GRID,
        type = "response"
      )
      prediction <- as.matrix(prediction)
      if (nrow(prediction) != length(validation) ||
          ncol(prediction) != length(SAP_LAMBDA_GRID)) {
        stop("Unexpected CV prediction dimensions")
      }
      loss_sum[j, ] <- loss_sum[j, ] +
        sap_binomial_deviance(y[validation], prediction)
      validation_rows[j] <- validation_rows[j] + length(validation)
    }
  }

  if (any(validation_rows != nrow(dt))) {
    stop("CV did not evaluate every row exactly once per imputation stream")
  }

  mean_deviance <- loss_sum / validation_rows
  selected_index <- apply(mean_deviance, 1L, which.min)
  if (any(selected_index %in% c(1L, length(SAP_LAMBDA_GRID)))) {
    stop("Selected lambda lies on the prespecified grid boundary")
  }
  selected_lambda <- SAP_LAMBDA_GRID[selected_index]
  list(
    lambda = selected_lambda,
    selected_index = selected_index,
    mean_deviance = mean_deviance,
    fold_log = data.table::rbindlist(mice_log)
  )
}

sap_fit_completed_models <- function(completed, y, knots, lambdas) {
  m <- length(completed)
  if (length(lambdas) != m) {
    stop("One lambda is required per completed dataset")
  }
  coefficient_matrix <- matrix(
    NA_real_,
    nrow = length(sap_design_column_names) + 1L,
    ncol = m,
    dimnames = list(c("(Intercept)", sap_design_column_names), paste0("m", seq_len(m)))
  )
  prediction_matrix <- matrix(
    NA_real_,
    nrow = length(y),
    ncol = m
  )
  for (j in seq_len(m)) {
    x <- sap_build_design(completed[[j]], knots)
    fit <- sap_fit_glmnet_path(x, y, lambda_grid = lambdas[j])
    coefficient_matrix[, j] <- as.numeric(
      as.matrix(coef(fit, s = lambdas[j]))
    )
    prediction_matrix[, j] <- as.numeric(
      predict(fit, newx = x, s = lambdas[j], type = "response")
    )
  }
  list(
    coefficients = coefficient_matrix,
    predictions = prediction_matrix,
    pooled_probability = rowMeans(prediction_matrix)
  )
}

sap_predict_completed_from_coefficients <- function(
  completed,
  coefficient_matrix,
  knots
) {
  m <- length(completed)
  if (ncol(coefficient_matrix) != m) {
    stop("Coefficient/imputation count mismatch")
  }
  prediction_matrix <- matrix(
    NA_real_,
    nrow = nrow(completed[[1L]]),
    ncol = m
  )
  for (j in seq_len(m)) {
    x <- sap_build_design(completed[[j]], knots)
    beta <- coefficient_matrix[, j]
    prediction_matrix[, j] <- plogis(
      beta[1L] + as.numeric(x %*% beta[-1L])
    )
  }
  list(
    predictions = prediction_matrix,
    pooled_probability = rowMeans(prediction_matrix)
  )
}

sap_auc <- function(y, score) {
  y <- as.integer(y)
  if (length(unique(y)) != 2L || anyNA(score)) {
    return(NA_real_)
  }
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

sap_calibration <- function(y, probability) {
  p <- pmin(pmax(probability, 1e-8), 1 - 1e-8)
  lp <- qlogis(p)
  fit <- suppressWarnings(
    glm.fit(
      x = cbind("(Intercept)" = 1, "logit_probability" = lp),
      y = as.numeric(y),
      family = binomial()
    )
  )
  coefficients <- fit$coefficients
  data.table::data.table(
    calibration_intercept = as.numeric(coefficients[1L]),
    calibration_slope = as.numeric(coefficients[2L])
  )
}

sap_performance <- function(y, probability) {
  calibration <- sap_calibration(y, probability)
  data.table::data.table(
    c_statistic = sap_auc(y, probability),
    calibration_intercept = calibration$calibration_intercept,
    calibration_slope = calibration$calibration_slope
  )
}

sap_pool_coefficient_table <- function(coefficient_matrix) {
  data.table::data.table(
    term = rownames(coefficient_matrix),
    pooled_log_odds = rowMeans(coefficient_matrix),
    between_imputation_sd = apply(coefficient_matrix, 1L, stats::sd),
    minimum_across_imputations = apply(coefficient_matrix, 1L, min),
    maximum_across_imputations = apply(coefficient_matrix, 1L, max),
    exponentiated_pooled_coefficient = exp(rowMeans(coefficient_matrix))
  )
}

sap_mice_log_table <- function(mice_result, context) {
  if (is.null(mice_result$logged_events)) {
    return(data.table::data.table(
      context = context,
      iteration = integer(),
      imputation = integer(),
      dependent = character(),
      method = character(),
      message = character()
    ))
  }
  logged <- data.table::as.data.table(mice_result$logged_events)
  data.table::setnames(
    logged,
    old = intersect(
      c("it", "im", "dep", "meth", "out"),
      names(logged)
    ),
    new = c(
      it = "iteration",
      im = "imputation",
      dep = "dependent",
      meth = "method",
      out = "message"
    )[intersect(c("it", "im", "dep", "meth", "out"), names(logged))]
  )
  logged[, context := context]
  logged
}
