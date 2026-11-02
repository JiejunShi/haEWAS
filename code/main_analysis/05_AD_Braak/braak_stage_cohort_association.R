#!/usr/bin/env Rscript
# Run single-CpG Braak association models, using binary disease status as the logistic outcome.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

options(stringsAsFactors = FALSE)

option_list <- list(
  make_option("--matrix", type = "character"),
  make_option("--phenotype", type = "character"),
  make_option("--output", type = "character"),
  make_option("--sample-column", type = "character", default = "Sample"),
  make_option("--stage-column", type = "character"),
  make_option("--covariates", type = "character"),
  make_option("--feature-column", type = "character", default = "feature_id"),
  make_option("--outcome-mode", type = "character", default = "advanced_binary"),
  make_option("--advanced-cutoff", type = "double", default = 5),
  make_option("--matrix-scale", type = "character", default = "proportion"),
  make_option("--epsilon", type = "double", default = 1e-6),
  make_option("--cores", type = "integer", default = 1L)
)
opts <- parse_args(OptionParser(option_list = option_list))
if (is.null(opts$matrix) || is.null(opts$phenotype) || is.null(opts$output) ||
    is.null(opts$stage_column) || is.null(opts$covariates) || opts$cores < 1L) {
  stop("--matrix, --phenotype, --stage-column, --covariates, and --output are required.")
}
if (!opts$outcome_mode %in% c("continuous", "advanced_binary")) {
  stop("--outcome-mode must be continuous or advanced_binary.")
}
if (!opts$matrix_scale %in% c("proportion", "m_value")) {
  stop("--matrix-scale must be proportion or m_value.")
}

read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(as.data.table(readRDS(path)))
  fread(path, sep = if (ext == "tsv") "\t" else ",", data.table = TRUE,
        check.names = FALSE, showProgress = FALSE)
}
normalize_id <- function(x) gsub("\\.", "-", trimws(as.character(x)))

if (tolower(tools::file_ext(opts$matrix)) == "rds") {
  methylation <- as.matrix(readRDS(opts$matrix))
} else {
  matrix_input <- fread(opts$matrix, data.table = FALSE, check.names = FALSE,
                        showProgress = FALSE)
  if (!opts$feature_column %in% names(matrix_input)) stop("The matrix lacks --feature-column.")
  rownames(matrix_input) <- as.character(matrix_input[[opts$feature_column]])
  methylation <- as.matrix(matrix_input[, setdiff(names(matrix_input), opts$feature_column), drop = FALSE])
}
storage.mode(methylation) <- "double"
if (is.null(rownames(methylation)) || is.null(colnames(methylation))) {
  stop("The matrix must be oriented as features by samples with dimnames.")
}
colnames(methylation) <- normalize_id(colnames(methylation))
if (anyDuplicated(rownames(methylation)) || anyDuplicated(colnames(methylation))) {
  stop("Duplicate feature or sample IDs are not allowed.")
}

phenotype <- read_table(opts$phenotype)
required <- c(opts$sample_column, opts$stage_column, trimws(strsplit(opts$covariates, ",", fixed = TRUE)[[1L]]))
if (!all(required %in% names(phenotype))) stop("A requested phenotype or covariate column is missing.")
phenotype[, (opts$sample_column) := normalize_id(get(opts$sample_column))]
if (anyDuplicated(phenotype[[opts$sample_column]])) stop("Duplicate phenotype sample IDs.")
covariates <- trimws(strsplit(opts$covariates, ",", fixed = TRUE)[[1L]])
covariates <- covariates[nzchar(covariates)]
if (!length(covariates) || opts$stage_column %in% covariates) {
  stop("List at least one covariate and do not include the stage column.")
}
common <- colnames(methylation)[colnames(methylation) %in% phenotype[[opts$sample_column]]]
if (length(common) < 3L) stop("Fewer than three samples overlap the matrix and phenotype.")
phenotype <- phenotype[match(common, get(opts$sample_column))]
stage <- suppressWarnings(as.numeric(as.character(phenotype[[opts$stage_column]])))
trait <- if (opts$outcome_mode == "continuous") stage else {
  ifelse(is.finite(stage), as.integer(stage >= opts$advanced_cutoff), NA_integer_)
}
if (length(unique(stats::na.omit(trait))) < 2L) stop("The outcome has fewer than two observed values.")

covariate_data <- as.data.frame(phenotype[, ..covariates])
for (name in names(covariate_data)) {
  if (is.character(covariate_data[[name]])) covariate_data[[name]] <- factor(covariate_data[[name]])
  if (length(unique(stats::na.omit(covariate_data[[name]]))) <= 1L) covariate_data[[name]] <- NULL
}
covariates <- names(covariate_data)
if (!length(covariates)) stop("No varying covariates remain after checking the cohort data.")
keep <- is.finite(trait) & complete.cases(covariate_data)
if (sum(keep) < 3L || length(unique(trait[keep])) < 2L) {
  stop("Insufficient complete samples or outcome variation after covariate filtering.")
}
common <- common[keep]
trait <- trait[keep]
covariate_data <- droplevels(covariate_data[keep, , drop = FALSE])
model_data <- data.frame(trait = trait, covariate_data, check.names = FALSE)
continuous_design <- model.matrix(reformulate(c("trait", covariates)), data = model_data)
trait_index <- match("trait", colnames(continuous_design))
if (is.na(trait_index)) stop("The trait term is not estimable in the design matrix.")
binary_design <- model.matrix(reformulate(covariates), data = model_data)

if (opts$matrix_scale == "proportion") {
  finite <- methylation[is.finite(methylation)]
  if (!length(finite)) stop("The methylation matrix contains no finite values.")
  methylation[is.finite(methylation)] <- pmin(
    pmax(methylation[is.finite(methylation)], opts$epsilon), 1 - opts$epsilon
  )
  methylation <- log2(methylation / (1 - methylation))
}
methylation <- methylation[, common, drop = FALSE]

fit_feature <- function(i) {
  response <- methylation[i, ]
  if (opts$outcome_mode == "advanced_binary") {
    valid <- is.finite(trait) & is.finite(response)
    n <- sum(valid)
    empty <- list(effect = NA_real_, se = NA_real_, p = NA_real_, n = n)
    if (n < 3L || length(unique(trait[valid])) < 2L ||
        length(unique(response[valid])) < 2L) return(empty)
    x <- cbind(binary_design[valid, , drop = FALSE], methylation = response[valid])
    fit <- tryCatch(
      suppressWarnings(glm.fit(x, trait[valid], family = binomial(),
                               control = glm.control(maxit = 100L))),
      error = function(e) NULL
    )
    feature_index <- ncol(x)
    if (is.null(fit) || !isTRUE(fit$converged) || fit$rank < ncol(x)) return(empty)
    estimate <- unname(fit$coefficients[feature_index])
    probability <- pmin(pmax(fit$fitted.values, 1e-8), 1 - 1e-8)
    information <- crossprod(x, x * (probability * (1 - probability)))
    inverse <- tryCatch(solve(information), error = function(e) NULL)
    if (is.null(inverse)) return(empty)
    se <- sqrt(inverse[feature_index, feature_index])
    if (!is.finite(estimate) || !is.finite(se) || se <= 0) return(empty)
    return(list(effect = estimate, se = se,
                p = 2 * pnorm(-abs(estimate / se)), n = n))
  }

  valid <- is.finite(response)
  n <- sum(valid)
  empty <- list(effect = NA_real_, se = NA_real_, p = NA_real_, n = n)
  if (n < 3L || length(unique(response[valid])) < 2L) return(empty)
  x <- continuous_design[valid, , drop = FALSE]
  y <- response[valid]
  fit <- tryCatch(lm.fit(x, y), error = function(e) NULL)
  if (is.null(fit) || fit$rank < 1L) return(empty)
  retained <- fit$qr$pivot[seq_len(fit$rank)]
  position <- match(trait_index, retained)
  df <- n - fit$rank
  if (is.na(position) || df < 1L) return(empty)
  reduced <- x[, retained, drop = FALSE]
  inverse <- tryCatch(solve(crossprod(reduced)), error = function(e) NULL)
  if (is.null(inverse)) return(empty)
  estimate <- unname(fit$coefficients[trait_index])
  se <- sqrt(sum(fit$residuals^2) / df * inverse[position, position])
  if (!is.finite(estimate) || !is.finite(se) || se <= 0) return(empty)
  list(effect = estimate, se = se, p = 2 * pt(-abs(estimate / se), df = df), n = n)
}

indices <- seq_len(nrow(methylation))
groups <- split(indices, rep(seq_len(min(opts$cores, length(indices))), length.out = length(indices)))
scan_group <- function(index) rbindlist(lapply(index, function(i) {
  as.data.table(c(list(feature_id = rownames(methylation)[i]), fit_feature(i)))
}), fill = TRUE)
if (.Platform$OS.type != "windows" && opts$cores > 1L) {
  result <- rbindlist(parallel::mclapply(groups, scan_group, mc.cores = opts$cores), fill = TRUE)
} else {
  result <- rbindlist(lapply(groups, scan_group), fill = TRUE)
}
result[, fdr := p.adjust(p, method = "BH")]
setorder(result, p)
dir.create(dirname(opts$output), recursive = TRUE, showWarnings = FALSE)
fwrite(result, opts$output, na = "NA")
