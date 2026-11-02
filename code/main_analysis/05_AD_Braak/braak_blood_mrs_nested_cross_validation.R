#!/usr/bin/env Rscript
# Evaluate binary brain-derived MRS models using repeated nested cross-validation.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--input_rds", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--endpoints", type = "character", default = "binary"),
  make_option("--top_grid", type = "character",
              default = "5,10,20,30,50,75,100,150,200,300,500,all"),
  make_option("--blend_grid", type = "character",
              default = "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"),
  make_option("--outer_folds", type = "integer", default = 5),
  make_option("--outer_repeats", type = "integer", default = 3),
  make_option("--inner_folds", type = "integer", default = 5),
  make_option("--min_features", type = "integer", default = 3),
  make_option("--bootstrap", type = "integer", default = 500),
  make_option("--cores", type = "integer", default = 1),
  make_option("--seed", type = "integer", default = 20260719)
)

opts <- parse_args(OptionParser(option_list = option_list))
for (name in c("input_rds", "output_dir")) {
  if (is.null(opts[[name]]) || !nzchar(opts[[name]])) stop("--", name, " is required.")
}
if (!file.exists(opts$input_rds)) stop("Input RDS not found: ", opts$input_rds)
set.seed(opts$seed)
setDTthreads(1L)

DIRS <- list(
  root = opts$output_dir,
  source = file.path(opts$output_dir, "source_data"),
  audit = file.path(opts$output_dir, "audit")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

ENDPOINTS <- trimws(strsplit(opts$endpoints, ",", fixed = TRUE)[[1]])
if (!identical(ENDPOINTS, "binary")) {
  stop("This analysis implements the binary advanced-Braak endpoint only.")
}
MODELS <- c("Baseline", "EWAS", "haEWAS", "Combined")
SCHEMES <- "Method_age_sex_EpiDISH"
MODEL_METHODS <- list(
  EWAS = "beta",
  haEWAS = c("chalm", "camda"),
  Combined = c("beta", "chalm", "camda")
)
BASE_COVARIATES <- c("blood_age", "blood_sex")
MAJOR_CELL_COVARIATES <- c(
  "blood_major_B", "blood_major_NK", "blood_major_CD4T",
  "blood_major_CD8T", "blood_major_Mono", "blood_major_Neutro",
  "blood_major_Eos"
)
BLOOD_COVARIATES <- c(BASE_COVARIATES, setdiff(MAJOR_CELL_COVARIATES, "blood_major_Eos"))
BLOOD_ASSOC_COVARIATES <- c(BLOOD_COVARIATES, "blood_slide")

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
normalize_id <- function(x) gsub("\\.", "-", trimws(as.character(x)))
parse_numeric_grid <- function(x, allow_all = FALSE) {
  tokens <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out <- suppressWarnings(as.numeric(tokens))
  if (allow_all) out[tolower(tokens) == "all"] <- 0
  sort(unique(out[is.finite(out)]))
}
top_grid <- as.integer(parse_numeric_grid(opts$top_grid, allow_all = TRUE))
blend_grid <- parse_numeric_grid(opts$blend_grid)
if (!length(top_grid) || !length(blend_grid)) {
  stop("A tuning grid is empty.")
}

bundle <- readRDS(opts$input_rds)
metadata <- as.data.table(bundle$feature_metadata)
samples <- as.data.table(bundle$sample_data)
blood_matrices <- bundle$method_matrices

required_metadata <- c(
  "endpoint", "cpg", "feature_id", "value_method", "brain_meta_weight", "brain_meta_p"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) stop("Feature metadata lacks: ", paste(missing_metadata, collapse = ", "))
required_samples <- c("blood_sample", "stage", "late_braak", "blood_slide",
                      MAJOR_CELL_COVARIATES, BLOOD_COVARIATES)
missing_samples <- setdiff(required_samples, names(samples))
if (length(missing_samples)) stop("Sample data lacks: ", paste(missing_samples, collapse = ", "))

samples[, blood_sample := normalize_id(blood_sample)]
samples[, stage := to_num(stage)]
samples[, late_braak := as.integer(to_num(late_braak))]
samples <- unique(samples[is.finite(stage) & late_braak %in% c(0L, 1L)], by = "blood_sample")
setorder(samples, blood_sample)

metadata[, `:=`(
  cpg = trimws(as.character(cpg)),
  value_method = tolower(trimws(as.character(value_method))),
  brain_meta_weight = to_num(brain_meta_weight),
  brain_meta_p = to_num(brain_meta_p)
)]
metadata <- metadata[
  endpoint == "binary" & nzchar(cpg) & is.finite(brain_meta_weight) &
    is.finite(brain_meta_p) & brain_meta_p >= 0 & brain_meta_p < 0.05
]
metadata <- unique(metadata, by = c("endpoint", "feature_id"))
metadata[, p_rank := frank(brain_meta_p, ties.method = "first"), by = "endpoint"]
if (!nrow(metadata)) stop("No brain-meta candidate feature remains.")
if (any(!metadata$value_method %in% c("beta", "chalm", "camda"))) {
  stop("The MRS input must contain only beta, CHALM and CAMDA features.")
}

align_method_matrix <- function(matrix, method, sample_ids, feature_ids) {
  colnames(matrix) <- normalize_id(colnames(matrix))
  keep_samples <- match(sample_ids, colnames(matrix))
  if (anyNA(keep_samples)) stop(method, " matrix lacks one or more blood samples.")
  matrix <- matrix[, keep_samples, drop = FALSE]
  rownames(matrix) <- paste(rownames(matrix), method, sep = "__")
  matrix <- matrix[intersect(feature_ids, rownames(matrix)), , drop = FALSE]
  t(matrix)
}

candidate_ids <- unique(metadata$feature_id)
blood_parts <- list()
matrix_methods <- intersect(c("beta", "chalm", "camda"), names(blood_matrices))
for (method in matrix_methods) {
  method_ids <- metadata[value_method == method, feature_id]
  blood_parts[[method]] <- align_method_matrix(
    blood_matrices[[method]], method, samples$blood_sample, method_ids
  )
}
blood_matrix <- do.call(cbind, blood_parts)
shared_features <- intersect(candidate_ids, colnames(blood_matrix))
blood_matrix <- blood_matrix[, shared_features, drop = FALSE]
metadata <- metadata[feature_id %in% shared_features]
rownames(blood_matrix) <- samples$blood_sample
if (!ncol(blood_matrix)) stop("No blood methylation feature remains.")

fwrite(
  metadata[, .(
    endpoint, cpg, feature_id, value_method, brain_meta_weight,
    brain_meta_p, p_rank
  )],
  file.path(DIRS$audit, "brain_meta_candidate_pool.csv")
)

mode_value <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("Unknown")
  names(sort(table(x), decreasing = TRUE))[[1]]
}

build_design <- function(train_idx, test_idx, covariates) {
  train_parts <- list()
  test_parts <- list()
  numeric_covariates <- setdiff(covariates, c("blood_sex", "brain_sex", "blood_slide", "brain_slide"))
  factor_covariates <- intersect(
    c("blood_sex", "brain_sex", "blood_slide", "brain_slide"), covariates
  )
  for (column in numeric_covariates) {
    if (!column %in% names(samples)) next
    train <- to_num(samples[[column]][train_idx])
    test <- to_num(samples[[column]][test_idx])
    fill <- median(train, na.rm = TRUE)
    if (!is.finite(fill)) next
    train[!is.finite(train)] <- fill
    test[!is.finite(test)] <- fill
    center <- mean(train)
    scale <- sd(train)
    if (!is.finite(scale) || scale <= 0) next
    train_parts[[column]] <- (train - center) / scale
    test_parts[[column]] <- (test - center) / scale
  }
  for (column in factor_covariates) {
    if (!column %in% names(samples)) next
    train <- trimws(as.character(samples[[column]][train_idx]))
    test <- trimws(as.character(samples[[column]][test_idx]))
    fill <- mode_value(train)
    train[is.na(train) | !nzchar(train)] <- fill
    test[is.na(test) | !nzchar(test)] <- fill
    levels <- names(sort(table(train), decreasing = TRUE))
    if (length(levels) <= 1L) next
    for (level in levels[-1L]) {
      name <- paste0(column, "_", make.names(level))
      train_parts[[name]] <- as.numeric(train == level)
      test_parts[[name]] <- as.numeric(test == level)
    }
  }
  train <- if (length(train_parts)) do.call(cbind, train_parts) else matrix(numeric(), length(train_idx), 0)
  test <- if (length(test_parts)) do.call(cbind, test_parts) else matrix(numeric(), length(test_idx), 0)
  list(
    train = cbind(intercept = rep(1, length(train_idx)), train),
    test = cbind(intercept = rep(1, length(test_idx)), test)
  )
}

independent_design <- function(x) {
  if (!ncol(x)) return(matrix(1, nrow(x), 1))
  q <- qr(x)
  x[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
}

make_folds <- function(outcome, k, seed) {
  k <- max(2L, min(as.integer(k), length(outcome) - 1L))
  set.seed(seed)
  folds <- integer(length(outcome))
  if (length(unique(outcome)) <= 10L) {
    strata <- outcome
  } else {
    breaks <- unique(quantile(
      outcome,
      probs = seq(0, 1, length.out = min(6L, length(outcome))),
      na.rm = TRUE
    ))
    strata <- if (length(breaks) >= 2L) {
      cut(outcome, breaks = breaks, include.lowest = TRUE)
    } else {
      rep("all", length(outcome))
    }
  }
  for (stratum in unique(strata)) {
    idx <- which(strata == stratum)
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

auc_rank <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]; score <- score[keep]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

spearman_value <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  if (sum(keep) < 4L || length(unique(y[keep])) < 2L || length(unique(score[keep])) < 2L) return(NA_real_)
  suppressWarnings(cor(y[keep], score[keep], method = "spearman"))
}

fit_outcome_model <- function(train_idx, test_idx, endpoint, score_train = NULL,
                              score_test = NULL, covariates = BLOOD_COVARIATES) {
  design <- build_design(train_idx, test_idx, covariates)
  train_x <- design$train
  test_x <- design$test
  if (!is.null(score_train)) {
    score_train <- as.matrix(score_train)
    score_test <- as.matrix(score_test)
    if (is.null(colnames(score_train))) {
      colnames(score_train) <- paste0("MRS", seq_len(ncol(score_train)))
    }
    colnames(score_test) <- colnames(score_train)
    if (!all(is.finite(score_train)) || !all(is.finite(score_test))) {
      return(rep(NA_real_, length(test_idx)))
    }
    colnames(score_train) <- paste0("MRS_", colnames(score_train))
    colnames(score_test) <- colnames(score_train)
    train_x <- cbind(train_x, score_train)
    test_x <- cbind(test_x, score_test)
  }
  train_x <- independent_design(train_x)
  test_x <- test_x[, colnames(train_x), drop = FALSE]
  y <- if (endpoint == "continuous") samples$stage[train_idx] else samples$late_braak[train_idx]
  if (endpoint == "continuous") {
    fit <- lm.fit(train_x, y)
    coefficients <- fit$coefficients
    coefficients[!is.finite(coefficients)] <- 0
    as.numeric(test_x %*% coefficients)
  } else {
    fit <- suppressWarnings(glm.fit(
      train_x, y, family = binomial(), control = glm.control(maxit = 100L)
    ))
    if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients))) {
      return(rep(NA_real_, length(test_idx)))
    }
    coefficients <- fit$coefficients
    as.numeric(plogis(test_x %*% coefficients))
  }
}

impute_matrix <- function(train, test) {
  fill <- vapply(seq_len(ncol(train)), function(j) {
    median(train[, j], na.rm = TRUE)
  }, numeric(1))
  names(fill) <- colnames(train)
  spread <- vapply(seq_len(ncol(train)), function(j) {
    sd(train[, j], na.rm = TRUE)
  }, numeric(1))
  names(spread) <- colnames(train)
  valid <- is.finite(fill) & is.finite(spread) & spread > 0
  fill[!is.finite(fill)] <- 0
  train_fill <- matrix(rep(fill, each = nrow(train)), nrow = nrow(train))
  test_fill <- matrix(rep(fill, each = nrow(test)), nrow = nrow(test))
  train_missing <- !is.finite(train)
  test_missing <- !is.finite(test)
  train[train_missing] <- train_fill[train_missing]
  test[test_missing] <- test_fill[test_missing]
  list(
    train = train[, valid, drop = FALSE],
    test = test[, valid, drop = FALSE],
    valid = names(fill)[valid]
  )
}

estimate_blood_effect <- function(x, y, design) {
  design <- independent_design(design)
  fit_one <- function(j) {
    feature_design <- cbind(design, methylation = x[, j])
    fit <- tryCatch(
      suppressWarnings(glm.fit(
        feature_design, y, family = binomial(),
        control = glm.control(maxit = 100L)
      )),
      error = function(e) NULL
    )
    if (is.null(fit) || !isTRUE(fit$converged) || fit$rank < ncol(feature_design)) {
      return(NA_real_)
    }
    effect <- unname(fit$coefficients[ncol(feature_design)])
    if (is.finite(effect)) effect else NA_real_
  }
  indices <- seq_len(ncol(x))
  effects <- if (.Platform$OS.type == "unix" && opts$cores > 1L) {
    parallel::mclapply(indices, fit_one, mc.cores = opts$cores)
  } else {
    lapply(indices, fit_one)
  }
  unlist(effects, use.names = FALSE)
}

prepare_fold <- function(train_idx, test_idx, endpoint_name, blood_covariates) {
  endpoint_meta <- metadata[endpoint == endpoint_name]
  feature_ids <- endpoint_meta$feature_id
  blood_raw_train <- blood_matrix[train_idx, feature_ids, drop = FALSE]
  blood_raw_test <- blood_matrix[test_idx, feature_ids, drop = FALSE]
  blood <- impute_matrix(blood_raw_train, blood_raw_test)
  if (!length(blood$valid)) return(NULL)

  blood_design <- build_design(train_idx, integer(), BLOOD_ASSOC_COVARIATES)$train
  disease_status <- samples$late_braak[train_idx]
  blood_effect <- estimate_blood_effect(
    blood$train, disease_status, blood_design
  )
  stats <- endpoint_meta[match(blood$valid, feature_id)]
  stats[, `:=`(
    blood_effect = blood_effect,
    direction_concordant = sign(brain_meta_weight) == sign(blood_effect)
  )]
  list(
    train_idx = train_idx, test_idx = test_idx,
    blood_train = blood$train, blood_test = blood$test,
    stats = stats, covariates = blood_covariates
  )
}

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}

finite_positive_fraction <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x > 0) else NA_real_
}

select_fold_features <- function(prepared, config, model) {
  methods <- MODEL_METHODS[[model]]
  top_k <- if ("top_k" %in% names(config)) {
    to_num(config$top_k[[1]])
  } else if ("requested_top_k" %in% names(config)) {
    to_num(config$requested_top_k[[1]])
  } else {
    NA_real_
  }
  if (length(top_k) != 1L || !is.finite(top_k) || top_k < 0) {
    stop("Invalid Top-N value in the selected configuration.")
  }
  top_k <- as.integer(top_k)
  selected <- prepared$stats[
    value_method %in% methods &
      is.finite(brain_meta_p) & brain_meta_p < 0.05 &
      is.finite(brain_meta_weight) & is.finite(blood_effect)
  ]
  setorder(selected, brain_meta_p, feature_id)
  selected <- unique(selected, by = "feature_id")
  if (top_k > 0L) {
    selected <- head(selected, top_k)
  }
  selected
}

predict_config <- function(prepared, config, model, endpoint) {
  selected <- select_fold_features(prepared, config, model)
  if (nrow(selected) < opts$min_features) return(NULL)
  brain_weight <- selected$brain_meta_weight
  blood_weight <- selected$blood_effect
  combined_weight <- config$blend_fraction * brain_weight +
    (1 - config$blend_fraction) * blood_weight
  if (any(!is.finite(combined_weight))) return(NULL)
  ids <- selected$feature_id
  score_train <- as.matrix(prepared$blood_train[, ids, drop = FALSE] %*% combined_weight)
  score_test <- as.matrix(prepared$blood_test[, ids, drop = FALSE] %*% combined_weight)
  score_name <- if (model == "EWAS") "EWAS_MRS" else "haEWAS_MRS"
  colnames(score_train) <- score_name
  colnames(score_test) <- score_name
  prediction <- fit_outcome_model(
    prepared$train_idx, prepared$test_idx, endpoint,
    score_train, score_test, prepared$covariates
  )
  selected[, `:=`(
    brain_weight = brain_weight,
    blood_weight = blood_weight,
    blended_weight = combined_weight
  )]
  list(
    prediction = prediction,
    selected = selected,
    score_components = colnames(score_train),
    score_train = score_train,
    score_test = score_test
  )
}

full_grid <- CJ(top_k = top_grid, blend_fraction = blend_grid)
full_grid[, `:=`(
  significance_rule = "P_0.05",
  direction_rule = "none"
)]
evaluated_grid <- copy(full_grid)
fwrite(full_grid, file.path(DIRS$audit, "complete_blended_weight_hyperparameter_space.csv"))
fwrite(evaluated_grid, file.path(DIRS$audit, "evaluated_blended_weight_configurations.csv"))
message(
  "[INFO] Blood samples: ", nrow(samples),
  "; P<0.05 binary brain-meta candidates: ", metadata[endpoint == "binary", uniqueN(feature_id)],
  "; evaluated configurations per score model: ", nrow(evaluated_grid),
  "; Top-N pooled across haEWAS measures; complete inner-CV grid"
)

evaluate_config <- function(config, prepared_folds, endpoint, model) {
  fold_rows <- lapply(seq_along(prepared_folds), function(i) {
    prepared <- prepared_folds[[i]]
    result <- predict_config(prepared, config, model, endpoint)
    if (is.null(result)) return(NULL)
    observed <- if (endpoint == "continuous") samples$stage[prepared$test_idx] else samples$late_braak[prepared$test_idx]
    metric <- if (endpoint == "continuous") {
      spearman_value(observed, result$prediction)
    } else {
      auc_rank(observed, result$prediction)
    }
    if (!is.finite(metric)) return(NULL)
    data.table(fold = i, metric = metric, n_features = nrow(result$selected))
  })
  fold_rows <- rbindlist(fold_rows)
  if (nrow(fold_rows) < max(2L, length(prepared_folds) - 1L)) return(NULL)
  data.table(
    significance_rule = config$significance_rule,
    direction_rule = config$direction_rule,
    requested_top_k = config$top_k,
    blend_fraction = config$blend_fraction,
    mean_metric = mean(fold_rows$metric),
    metric_se = if (nrow(fold_rows) > 1L) sd(fold_rows$metric) / sqrt(nrow(fold_rows)) else NA_real_,
    mean_features = mean(fold_rows$n_features),
    min_features = min(fold_rows$n_features),
    valid_inner_folds = nrow(fold_rows)
  )
}

select_tuned_config <- function(evaluated) {
  if (!nrow(evaluated)) return(NULL)
  evaluated[, lambda_distance_from_equal := abs(blend_fraction - 0.5)]
  setorder(evaluated, -mean_metric, mean_features, lambda_distance_from_equal,
           requested_top_k, blend_fraction)
  evaluated[1]
}

tune_score_model <- function(prepared_folds, endpoint, model) {
  rows <- split(evaluated_grid, seq_len(nrow(evaluated_grid)))
  worker <- function(config) evaluate_config(config, prepared_folds, endpoint, model)
  evaluated <- if (.Platform$OS.type == "unix" && opts$cores > 1L) {
    parallel::mclapply(rows, worker, mc.cores = opts$cores)
  } else {
    lapply(rows, worker)
  }
  evaluated <- rbindlist(evaluated, use.names = TRUE, fill = TRUE)
  if (nrow(evaluated)) evaluated[, search_stage := "inner_cv"]
  best <- select_tuned_config(evaluated)
  if (is.null(best)) return(NULL)
  list(best = best, evaluated = evaluated)
}

predict_component_joint <- function(prepared, base_config, haewas_config, endpoint,
                                    base_model = "EWAS") {
  base_result <- predict_config(prepared, base_config, base_model, endpoint)
  haewas <- predict_config(prepared, haewas_config, "haEWAS", endpoint)
  if (is.null(base_result) || is.null(haewas)) return(NULL)

  prefix_scores <- function(x, prefix) {
    x <- as.matrix(x)
    colnames(x) <- paste0(prefix, "_", colnames(x))
    x
  }
  score_train <- cbind(
    prefix_scores(base_result$score_train, base_model),
    prefix_scores(haewas$score_train, "haEWAS")
  )
  score_test <- cbind(
    prefix_scores(base_result$score_test, base_model),
    prefix_scores(haewas$score_test, "haEWAS")
  )
  prediction <- fit_outcome_model(
    prepared$train_idx, prepared$test_idx, endpoint,
    score_train, score_test, prepared$covariates
  )
  selected <- rbindlist(list(
    copy(base_result$selected)[, component_model := base_model],
    copy(haewas$selected)[, component_model := "haEWAS"]
  ), use.names = TRUE, fill = TRUE)
  list(
    prediction = prediction,
    selected = selected,
    score_components = colnames(score_train)
  )
}

evaluate_component_joint <- function(prepared_folds, base_config, haewas_config, endpoint,
                                     base_model = "EWAS") {
  rows <- lapply(seq_along(prepared_folds), function(i) {
    result <- predict_component_joint(
      prepared_folds[[i]], base_config, haewas_config, endpoint, base_model
    )
    if (is.null(result)) return(NULL)
    observed <- if (endpoint == "continuous") {
      samples$stage[prepared_folds[[i]]$test_idx]
    } else {
      samples$late_braak[prepared_folds[[i]]$test_idx]
    }
    metric <- if (endpoint == "continuous") {
      spearman_value(observed, result$prediction)
    } else {
      auc_rank(observed, result$prediction)
    }
    if (!is.finite(metric)) return(NULL)
    data.table(fold = i, metric = metric, n_features = nrow(result$selected))
  })
  rows <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (nrow(rows) < max(2L, length(prepared_folds) - 1L)) return(NULL)
  data.table(
    mean_metric = mean(rows$metric),
    metric_se = if (nrow(rows) > 1L) sd(rows$metric) / sqrt(nrow(rows)) else NA_real_,
    mean_features = mean(rows$n_features),
    valid_inner_folds = nrow(rows)
  )
}

bootstrap_ci <- function(y, prediction, endpoint) {
  keep <- is.finite(y) & is.finite(prediction)
  y <- y[keep]; prediction <- prediction[keep]
  draws <- replicate(opts$bootstrap, {
    if (endpoint == "binary") {
      idx <- c(
        sample(which(y == 0), sum(y == 0), replace = TRUE),
        sample(which(y == 1), sum(y == 1), replace = TRUE)
      )
      auc_rank(y[idx], prediction[idx])
    } else {
      idx <- sample(seq_along(y), length(y), replace = TRUE)
      spearman_value(y[idx], prediction[idx])
    }
  })
  draws <- draws[is.finite(draws)]
  if (!length(draws)) return(c(NA_real_, NA_real_))
  as.numeric(quantile(draws, c(0.025, 0.975), names = FALSE))
}

roc_points <- function(y, score) {
  keep <- is.finite(y) & is.finite(score)
  y <- y[keep]; score <- score[keep]
  thresholds <- c(Inf, sort(unique(score), decreasing = TRUE), -Inf)
  rbindlist(lapply(thresholds, function(threshold) {
    predicted <- score >= threshold
    data.table(
      threshold = threshold,
      fpr = sum(predicted & y == 0) / sum(y == 0),
      tpr = sum(predicted & y == 1) / sum(y == 1)
    )
  }))
}

outer_predictions <- list()
tuning_results <- list()
selected_feature_rows <- list()

for (endpoint_index in seq_along(ENDPOINTS)) {
  endpoint_name <- ENDPOINTS[[endpoint_index]]
  endpoint_seed_index <- match(endpoint_name, c("continuous", "binary"))
  outcome <- if (endpoint_name == "continuous") samples$stage else samples$late_braak
  for (repeat_index in seq_len(opts$outer_repeats)) {
    outer_foldid <- make_folds(
      outcome, opts$outer_folds,
      opts$seed + endpoint_seed_index * 100000L + repeat_index * 1000L
    )
    for (outer_fold in sort(unique(outer_foldid))) {
      outer_test <- which(outer_foldid == outer_fold)
      outer_train <- which(outer_foldid != outer_fold)
      inner_foldid <- make_folds(
        outcome[outer_train], opts$inner_folds,
        opts$seed + endpoint_seed_index * 100000L + repeat_index * 1000L + outer_fold
      )
      for (scheme in SCHEMES) {
        covariates <- BLOOD_COVARIATES
        selected_cells <- "EpiDISH_B_NK_CD4T_CD8T_Mono_Neutro;Eos_reference"

        message(
          "[OUTER] endpoint=", endpoint_name,
          "; repeat=", repeat_index, "/", opts$outer_repeats,
          "; fold=", outer_fold, "/", max(outer_foldid),
          "; covariates=", scheme,
          "; cells=", selected_cells
        )
        baseline <- fit_outcome_model(outer_train, outer_test, endpoint_name, covariates = covariates)
        baseline_key <- paste(endpoint_name, repeat_index, outer_fold, scheme, "Baseline", sep = "__")
        outer_predictions[[baseline_key]] <- data.table(
          endpoint = endpoint_name, outer_repeat = repeat_index, outer_fold = outer_fold,
          covariate_scheme = scheme, model = "Baseline",
          blood_sample = samples$blood_sample[outer_test], observed = outcome[outer_test],
          prediction = baseline, selected_cell_subset = selected_cells
        )

        inner_folds <- sort(unique(inner_foldid))
        prepared_folds <- lapply(inner_folds, function(inner_fold) {
          train_idx <- outer_train[inner_foldid != inner_fold]
          test_idx <- outer_train[inner_foldid == inner_fold]
          prepare_fold(train_idx, test_idx, endpoint_name, covariates)
        })
        if (any(vapply(prepared_folds, is.null, logical(1)))) next
        outer_prepared <- prepare_fold(outer_train, outer_test, endpoint_name, covariates)
        if (is.null(outer_prepared)) next



        models_to_tune <- c("EWAS", "haEWAS")
        best_configs <- list()
        final_results <- list()
        for (model_name in models_to_tune) {
          tuning_result <- tune_score_model(prepared_folds, endpoint_name, model_name)
          if (is.null(tuning_result)) next
          evaluated <- tuning_result$evaluated
          best <- tuning_result$best
          final_result <- predict_config(outer_prepared, best, model_name, endpoint_name)
          if (is.null(final_result)) next
          best_configs[[model_name]] <- best
          final_results[[model_name]] <- final_result
          message(
            "[SELECT] endpoint=", endpoint_name,
            "; model=", model_name,
            "; top=", ifelse(best$requested_top_k == 0L, "all", best$requested_top_k),
            "; brain_fraction=", best$blend_fraction,
            "; features=", nrow(final_result$selected),
            "; score_components=", paste(final_result$score_components, collapse = "+")
          )
          key <- paste(endpoint_name, repeat_index, outer_fold, scheme, model_name, sep = "__")
          outer_predictions[[key]] <- data.table(
            endpoint = endpoint_name, outer_repeat = repeat_index, outer_fold = outer_fold,
            covariate_scheme = scheme, model = model_name,
            blood_sample = samples$blood_sample[outer_test], observed = outcome[outer_test],
            prediction = final_result$prediction,
            significance_rule = best$significance_rule,
            direction_rule = best$direction_rule,
            requested_top_k = best$requested_top_k,
            blend_fraction = best$blend_fraction,
            n_features = nrow(final_result$selected),
            n_score_components = length(final_result$score_components),
            score_components = paste(final_result$score_components, collapse = ";"),
            selected_cell_subset = selected_cells
          )
          evaluated[, `:=`(
            endpoint = endpoint_name, outer_repeat = repeat_index, outer_fold = outer_fold,
            covariate_scheme = scheme, model = model_name,
            selected = significance_rule == best$significance_rule &
              direction_rule == best$direction_rule &
              requested_top_k == best$requested_top_k & blend_fraction == best$blend_fraction
          )]
          tuning_results[[key]] <- evaluated
          selected_features <- final_result$selected
          selected_features[, `:=`(
            endpoint = endpoint_name, outer_repeat = repeat_index, outer_fold = outer_fold,
            covariate_scheme = scheme, model = model_name,
            significance_rule = best$significance_rule,
            direction_rule = best$direction_rule,
            requested_top_k = best$requested_top_k,
            blend_fraction = best$blend_fraction,
            selected_cell_subset = selected_cells
          )]
          selected_feature_rows[[key]] <- selected_features
        }

        for (component_baseline in "EWAS") {
            if (!all(c(component_baseline, "haEWAS") %in% names(best_configs)) ||
                !all(c(component_baseline, "haEWAS") %in% names(final_results))) next
            base_config <- best_configs[[component_baseline]]
            ha_config <- best_configs$haEWAS
            joint_inner <- evaluate_component_joint(
              prepared_folds, base_config, ha_config, endpoint_name, component_baseline
            )
            joint_final <- predict_component_joint(
              outer_prepared, base_config, ha_config, endpoint_name, component_baseline
            )
            if (is.null(joint_inner) || is.null(joint_final)) next
            joint_model_name <- "Combined"
            combined_key <- paste(
              endpoint_name, repeat_index, outer_fold, scheme, joint_model_name, sep = "__"
            )
            outer_predictions[[combined_key]] <- data.table(
              endpoint = endpoint_name,
              outer_repeat = repeat_index,
              outer_fold = outer_fold,
              covariate_scheme = scheme,
              model = joint_model_name,
              blood_sample = samples$blood_sample[outer_test],
              observed = outcome[outer_test],
              prediction = joint_final$prediction,
              significance_rule = "component_specific",
              direction_rule = "component_specific",
              requested_top_k = NA_integer_,
              blend_fraction = NA_real_,
              ewas_blend_fraction = base_config$blend_fraction,
              haewas_blend_fraction = ha_config$blend_fraction,
              base_component_model = component_baseline,
              n_features = nrow(joint_final$selected),
              n_score_components = length(joint_final$score_components),
              score_components = paste(joint_final$score_components, collapse = ";"),
              selected_cell_subset = selected_cells
            )
            joint_audit <- copy(joint_inner)
            joint_audit[, `:=`(
              endpoint = endpoint_name,
              outer_repeat = repeat_index,
              outer_fold = outer_fold,
              covariate_scheme = scheme,
              model = joint_model_name,
              significance_rule = "component_specific",
              direction_rule = "component_specific",
              requested_top_k = NA_integer_,
              blend_fraction = NA_real_,
              ewas_blend_fraction = base_config$blend_fraction,
              haewas_blend_fraction = ha_config$blend_fraction,
              base_component_model = component_baseline,
              search_stage = "component_joint",
              selected = TRUE
            )]
            tuning_results[[combined_key]] <- joint_audit

            selected_features <- copy(joint_final$selected)
            selected_features[, `:=`(
              endpoint = endpoint_name,
              outer_repeat = repeat_index,
              outer_fold = outer_fold,
              covariate_scheme = scheme,
              model = joint_model_name,
              significance_rule = fifelse(
                component_model == component_baseline,
                base_config$significance_rule,
                ha_config$significance_rule
              ),
              direction_rule = fifelse(
                component_model == component_baseline,
                base_config$direction_rule,
                ha_config$direction_rule
              ),
              requested_top_k = fifelse(
                component_model == component_baseline,
                base_config$requested_top_k,
                ha_config$requested_top_k
              ),
              blend_fraction = fifelse(
                component_model == component_baseline,
                base_config$blend_fraction,
                ha_config$blend_fraction
              ),
              selected_cell_subset = selected_cells
            )]
            selected_feature_rows[[combined_key]] <- selected_features
            message(
              "[COMPONENT JOINT] endpoint=", endpoint_name,
              "; model=", joint_model_name,
              "; base=", component_baseline,
              "; base lambda=", sprintf("%.3f", base_config$blend_fraction),
              "; haEWAS lambda=", sprintf("%.3f", ha_config$blend_fraction),
              "; features=", nrow(joint_final$selected),
              "; score_components=", paste(joint_final$score_components, collapse = "+")
            )
        }

      }
    }
  }
}

raw_predictions <- rbindlist(outer_predictions, use.names = TRUE, fill = TRUE)
tuning <- rbindlist(tuning_results, use.names = TRUE, fill = TRUE)
selected_features <- rbindlist(selected_feature_rows, use.names = TRUE, fill = TRUE)
if (!nrow(raw_predictions)) stop("No outer-CV prediction was generated.")
prediction_defaults <- list(
  significance_rule = NA_character_, direction_rule = NA_character_,
  requested_top_k = NA_integer_,
  blend_fraction = NA_real_,
  ewas_blend_fraction = NA_real_, haewas_blend_fraction = NA_real_,
  selected_cell_subset = NA_character_, n_features = NA_integer_
)
for (column in names(prediction_defaults)) {
  if (!column %in% names(raw_predictions)) {
    set(raw_predictions, j = column, value = prediction_defaults[[column]])
  }
}

tuning_defaults <- list(
  endpoint = NA_character_, outer_repeat = NA_integer_, outer_fold = NA_integer_,
  covariate_scheme = NA_character_, model = NA_character_,
  blend_fraction = NA_real_, mean_metric = NA_real_, mean_features = NA_real_
)
for (column in names(tuning_defaults)) {
  if (!column %in% names(tuning)) {
    set(tuning, j = column, value = tuning_defaults[[column]])
  }
}

feature_defaults <- list(
  cpg = NA_character_, value_method = NA_character_,
  brain_weight = NA_real_, blood_weight = NA_real_,
  blended_weight = NA_real_
)
for (column in names(feature_defaults)) {
  if (!column %in% names(selected_features)) {
    set(selected_features, j = column, value = feature_defaults[[column]])
  }
}

predictions <- raw_predictions[, .(
  observed = observed[[1]],
  prediction = mean(prediction, na.rm = TRUE),
  prediction_sd = sd(prediction, na.rm = TRUE),
  n_outer_predictions = sum(is.finite(prediction))
), by = .(endpoint, covariate_scheme, model, blood_sample)]

performance <- predictions[, {
  keep <- is.finite(observed) & is.finite(prediction)
  ci <- bootstrap_ci(observed[keep], prediction[keep], endpoint[[1]])
  if (endpoint[[1]] == "continuous") {
    tss <- sum((observed[keep] - mean(observed[keep]))^2)
    list(
      n = sum(keep), spearman = spearman_value(observed, prediction),
      pearson = cor(observed[keep], prediction[keep]),
      r2_external = 1 - sum((observed[keep] - prediction[keep])^2) / tss,
      rmse = sqrt(mean((observed[keep] - prediction[keep])^2)),
      mae = mean(abs(observed[keep] - prediction[keep])), auc = NA_real_,
      primary_ci_low = ci[[1]], primary_ci_high = ci[[2]]
    )
  } else {
    list(
      n = sum(keep), spearman = NA_real_, pearson = NA_real_,
      r2_external = NA_real_, rmse = NA_real_, mae = NA_real_,
      auc = auc_rank(observed, prediction),
      primary_ci_low = ci[[1]], primary_ci_high = ci[[2]]
    )
  }
}, by = .(endpoint, covariate_scheme, model)]

metric_difference <- function(y, first, second, endpoint_name) {
  if (endpoint_name == "continuous") {
    spearman_value(y, first) - spearman_value(y, second)
  } else {
    auc_rank(y, first) - auc_rank(y, second)
  }
}

comparison_pairs <- data.table(
  first_model = c("EWAS", "haEWAS", "Combined", "haEWAS", "Combined", "Combined"),
  second_model = c("Baseline", "Baseline", "Baseline", "EWAS", "EWAS", "haEWAS")
)

paired_differences <- rbindlist(lapply(ENDPOINTS, function(endpoint_name) {
  rbindlist(lapply(SCHEMES, function(scheme) {
    wide <- dcast(
      predictions[endpoint == endpoint_name & covariate_scheme == scheme],
      blood_sample + observed ~ model,
      value.var = "prediction"
    )
    rbindlist(lapply(seq_len(nrow(comparison_pairs)), function(i) {
      first_model <- comparison_pairs$first_model[[i]]
      second_model <- comparison_pairs$second_model[[i]]
      if (!all(c(first_model, second_model) %in% names(wide))) return(NULL)
      y <- wide$observed
      first <- wide[[first_model]]
      second <- wide[[second_model]]
      keep <- is.finite(y) & is.finite(first) & is.finite(second)
      y <- y[keep]; first <- first[keep]; second <- second[keep]
      point <- metric_difference(y, first, second, endpoint_name)
      draws <- replicate(opts$bootstrap, {
        idx <- if (endpoint_name == "binary") {
          c(
            sample(which(y == 0), sum(y == 0), replace = TRUE),
            sample(which(y == 1), sum(y == 1), replace = TRUE)
          )
        } else {
          sample(seq_along(y), length(y), replace = TRUE)
        }
        metric_difference(y[idx], first[idx], second[idx], endpoint_name)
      })
      draws <- draws[is.finite(draws)]
      ci <- if (length(draws)) {
        quantile(draws, c(0.025, 0.975), names = FALSE)
      } else {
        c(NA_real_, NA_real_)
      }
      p_value <- if (length(draws)) {
        min(1, 2 * min(
          (sum(draws <= 0) + 1) / (length(draws) + 1),
          (sum(draws >= 0) + 1) / (length(draws) + 1)
        ))
      } else {
        NA_real_
      }
      data.table(
        endpoint = endpoint_name,
        covariate_scheme = scheme,
        first_model = first_model,
        second_model = second_model,
        contrast = paste0(first_model, " - ", second_model),
        metric = ifelse(endpoint_name == "continuous", "Spearman", "AUC"),
        difference = point,
        ci_low = ci[[1]],
        ci_high = ci[[2]],
        bootstrap_two_sided_p = p_value,
        bootstrap_reps = length(draws)
      )
    }), use.names = TRUE, fill = TRUE)
  }), use.names = TRUE, fill = TRUE)
}), use.names = TRUE, fill = TRUE)

roc_data <- predictions[endpoint == "binary", roc_points(observed, prediction),
                        by = .(endpoint, covariate_scheme, model)]

hyperparameter_frequency <- raw_predictions[model != "Baseline", .(
  outer_splits = uniqueN(paste(outer_repeat, outer_fold, sep = "__"))
), by = .(
  endpoint, covariate_scheme, model, significance_rule, direction_rule,
  requested_top_k, blend_fraction, selected_cell_subset
)]

lambda_profile <- tuning[
  model %in% c("EWAS", "haEWAS") &
    is.finite(blend_fraction) & is.finite(mean_metric)
]
if (nrow(lambda_profile)) {
  setorder(
    lambda_profile,
    endpoint, covariate_scheme, model, outer_repeat, outer_fold,
    blend_fraction, -mean_metric, mean_features
  )
  lambda_profile <- lambda_profile[, .SD[1], by = .(
    endpoint, covariate_scheme, model, outer_repeat, outer_fold, blend_fraction
  )]
  setnames(lambda_profile, "blend_fraction", "brain_weight_fraction_lambda")
}
lambda_summary <- if (nrow(lambda_profile)) {
  lambda_profile[, .(
    mean_inner_metric = mean(mean_metric),
    se_across_outer_splits = if (.N > 1L) sd(mean_metric) / sqrt(.N) else NA_real_,
    n_outer_splits = .N
  ), by = .(
    endpoint, covariate_scheme, model, brain_weight_fraction_lambda
  )]
} else data.table()

selected_lambdas <- unique(
  raw_predictions[
    model %in% c("EWAS", "haEWAS") & is.finite(blend_fraction),
    .(
      endpoint, covariate_scheme, model, outer_repeat, outer_fold,
      brain_weight_fraction_lambda = blend_fraction,
      blood_weight_fraction = 1 - blend_fraction,
      significance_rule, direction_rule,
      requested_top_k, n_features, selected_cell_subset
    )
  ],
  by = c("endpoint", "covariate_scheme", "model", "outer_repeat", "outer_fold")
)
lambda_totals <- selected_lambdas[, .(
  total_outer_splits = uniqueN(paste(outer_repeat, outer_fold, sep = "__"))
), by = .(endpoint, covariate_scheme, model)]
lambda_frequency <- selected_lambdas[, .(
  selected_outer_splits = .N
), by = .(
  endpoint, covariate_scheme, model,
  brain_weight_fraction_lambda, blood_weight_fraction
)]
lambda_frequency <- merge(
  lambda_frequency, lambda_totals,
  by = c("endpoint", "covariate_scheme", "model"), all.x = TRUE
)
lambda_frequency[, selected_fraction := selected_outer_splits / total_outer_splits]
setorder(
  lambda_frequency,
  endpoint, covariate_scheme, model, -selected_outer_splits,
  brain_weight_fraction_lambda
)

feature_stability <- if (nrow(selected_features)) {
  denominator <- raw_predictions[model != "Baseline", .(
    outer_splits = uniqueN(paste(outer_repeat, outer_fold, sep = "__"))
  ), by = .(endpoint, covariate_scheme, model)]
  out <- selected_features[, .(
    selected_splits = uniqueN(paste(outer_repeat, outer_fold, sep = "__")),
    median_brain_weight = finite_median(brain_weight),
    median_blood_weight = finite_median(blood_weight),
    median_blended_weight = finite_median(blended_weight),
    positive_blended_fraction = finite_positive_fraction(blended_weight)
  ), by = .(endpoint, covariate_scheme, model, feature_id, cpg, value_method)]
  out <- merge(out, denominator, by = c("endpoint", "covariate_scheme", "model"), all.x = TRUE)
  out[, selection_frequency := selected_splits / outer_splits]
  setorder(out, endpoint, covariate_scheme, model, -selection_frequency, feature_id)
  out
} else data.table()

cell_frequency <- raw_predictions[model == "Baseline", .(
  outer_splits = uniqueN(paste(outer_repeat, outer_fold, sep = "__"))
), by = .(endpoint, covariate_scheme, selected_cell_subset)]

fwrite(raw_predictions, file.path(DIRS$source, "blended_mrs_outer_predictions_all_repeats.csv"))
fwrite(predictions, file.path(DIRS$source, "blended_mrs_outer_predictions_aggregated.csv"))
fwrite(performance, file.path(DIRS$source, "blended_mrs_repeated_nested_cv_performance.csv"))
fwrite(paired_differences, file.path(DIRS$source, "blended_mrs_paired_performance_differences.csv"))
fwrite(roc_data, file.path(DIRS$source, "blended_mrs_roc_points.csv"))
fwrite(tuning, file.path(DIRS$audit, "blended_mrs_inner_cv_tuning_results.csv"))
fwrite(hyperparameter_frequency, file.path(DIRS$audit, "blended_mrs_selected_hyperparameter_frequency.csv"))
fwrite(lambda_profile, file.path(DIRS$audit, "brain_blood_lambda_inner_cv_profile.csv"))
fwrite(lambda_summary, file.path(DIRS$audit, "brain_blood_lambda_inner_cv_summary.csv"))
fwrite(selected_lambdas, file.path(DIRS$audit, "brain_blood_lambda_selected_by_outer_split.csv"))
fwrite(lambda_frequency, file.path(DIRS$audit, "brain_blood_lambda_selection_frequency.csv"))
fwrite(selected_features, file.path(DIRS$audit, "blended_mrs_selected_features_all_outer_splits.csv"))
fwrite(feature_stability, file.path(DIRS$audit, "blended_mrs_feature_stability.csv"))
fwrite(cell_frequency, file.path(DIRS$audit, "blood_cell_covariate_scheme_frequency.csv"))

manifest <- data.table(
  field = c("n_blood", "n_advanced", "n_nonadvanced", "candidate_p",
            "blood_effect_model", "weight_formula", "mrs_formula", "haEWAS_score",
            "combined_score", "outer_cv", "inner_folds", "top_n_grid",
            "lambda_grid", "baseline_covariates", "cell_reference"),
  value = c(
    as.character(nrow(samples)),
    as.character(sum(samples$late_braak == 1L)),
    as.character(sum(samples$late_braak == 0L)),
    "0.05",
    "status~M_value+age+sex+EpiDISH+slide",
    "lambda*beta_brain+(1-lambda)*beta_blood",
    "sum(M_value*weight)",
    "pooled_CHALM_CAMDA",
    "EWAS_MRS+haEWAS_MRS",
    paste0(opts$outer_repeats, "x", opts$outer_folds),
    as.character(opts$inner_folds),
    paste(ifelse(top_grid == 0L, "all", top_grid), collapse = ","),
    paste(blend_grid, collapse = ","),
    paste(c("blood_age", "blood_sex", setdiff(MAJOR_CELL_COVARIATES, "blood_major_Eos")), collapse = ";"),
    "blood_major_Eos"
  )
)
fwrite(manifest, file.path(DIRS$root, "analysis_manifest.csv"))
runtime_warnings <- warnings()
if (!is.null(runtime_warnings)) {
  writeLines(
    capture.output(print(runtime_warnings)),
    file.path(DIRS$audit, "runtime_warnings.txt")
  )
}
writeLines(
  paste("completed", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  file.path(DIRS$root, "brain_blood_blended_mrs.done")
)

message("[DONE] Brain-blood blended MRS nested validation: ", opts$output_dir)
print(performance)
print(paired_differences)
