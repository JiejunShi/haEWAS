# Run CellDMC and export phenotype-specific EWAS and haEWAS results.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(EpiDISH)
})

option_list <- list(
  make_option("--beta-matrix", type = "character"),
  make_option("--chalm-matrix", type = "character"),
  make_option("--camda-matrix", type = "character"),
  make_option("--phenotype", type = "character"),
  make_option("--cell-fractions", type = "character"),
  make_option("--output-dir", type = "character"),
  make_option("--outcome-column", type = "character"),
  make_option("--case", type = "character"),
  make_option("--control", type = "character")
)
opts <- parse_args(OptionParser(option_list = option_list))
required <- c("beta_matrix", "chalm_matrix", "camda_matrix", "phenotype",
              "cell_fractions", "output_dir", "outcome_column", "case",
              "control")
if (any(vapply(required, function(name) is.null(opts[[name]]) || !nzchar(opts[[name]]), logical(1)))) {
  stop("Three matrices, phenotype, cell fractions, outcome, case/control, and output are required.")
}

sample_column <- "Sample"
cell_columns <- c("Hep", "Chol", "EC", "Kup", "Lym")
covariates <- c("Age", "Sex", paste0("PC_resid", 1:5))
sex_column <- "Sex"
fdr_threshold <- 0.05

normalize_id <- function(x) gsub("\\.", "-", trimws(as.character(x)))

read_matrix <- function(path) {
  if (tolower(tools::file_ext(path)) == "rds") {
    x <- as.matrix(readRDS(path))
  } else {
    input <- fread(path, data.table = FALSE, check.names = FALSE,
                   showProgress = FALSE)
    if (ncol(input) < 2L) stop("A methylation matrix needs feature and sample columns: ", path)
    feature <- trimws(as.character(input[[1L]]))
    x <- as.matrix(input[, -1L, drop = FALSE])
    rownames(x) <- feature
  }
  storage.mode(x) <- "double"
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("Methylation matrices must have feature and sample dimnames: ", path)
  }
  colnames(x) <- normalize_id(colnames(x))
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop("Duplicate feature or sample IDs in ", path)
  }
  finite <- x[is.finite(x)]
  if (!length(finite) || any(finite < -1e-8 | finite > 1 + 1e-8)) {
    stop("CellDMC input must be an imputed proportion-scale matrix in [0, 1]: ", path)
  }
  if (any(!is.finite(x))) stop("Methylation matrices must be imputed before CellDMC: ", path)
  pmin(pmax(x, 0), 1)
}

matrices <- list(
  EWAS = read_matrix(opts$beta_matrix),
  CHALM = read_matrix(opts$chalm_matrix),
  CAMDA = read_matrix(opts$camda_matrix)
)

read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(as.data.table(readRDS(path)))
  fread(path, sep = if (ext == "tsv") "\t" else ",", data.table = TRUE,
        check.names = FALSE, showProgress = FALSE)
}

phenotype <- read_table(opts$phenotype)
fractions <- read_table(opts$cell_fractions)
if (!all(c(sample_column, opts$outcome_column, covariates) %in% names(phenotype))) {
  stop("Phenotype or covariate columns are missing.")
}
if (!all(c(sample_column, cell_columns) %in% names(fractions))) {
  stop("Required liver cell-fraction columns are missing.")
}
phenotype[, (sample_column) := normalize_id(get(sample_column))]
fractions[, (sample_column) := normalize_id(get(sample_column))]
if (anyDuplicated(phenotype[[sample_column]]) || anyDuplicated(fractions[[sample_column]])) {
  stop("Duplicate sample IDs in phenotype or cell fractions.")
}

samples <- colnames(matrices$EWAS)
for (method in c("CHALM", "CAMDA")) {
  samples <- samples[samples %in% colnames(matrices[[method]])]
}
samples <- samples[samples %in% phenotype[[sample_column]] &
                     samples %in% fractions[[sample_column]]]
if (length(samples) < 3L) stop("Fewer than three samples overlap all inputs.")
phenotype <- phenotype[match(samples, get(sample_column))]
fractions <- fractions[match(samples, get(sample_column))]
if (!identical(phenotype[[sample_column]], fractions[[sample_column]])) {
  stop("Phenotype and cell-fraction sample order differs.")
}

outcome <- trimws(as.character(phenotype[[opts$outcome_column]]))
cell_data <- as.data.frame(fractions[, ..cell_columns])
for (name in cell_columns) cell_data[[name]] <- suppressWarnings(as.numeric(as.character(cell_data[[name]])))
covariate_data <- as.data.frame(phenotype[, ..covariates])
for (name in covariates) {
  if (name == sex_column) {
    covariate_data[[name]] <- factor(trimws(as.character(covariate_data[[name]])))
  } else {
    covariate_data[[name]] <- suppressWarnings(as.numeric(as.character(covariate_data[[name]])))
  }
}
keep <- outcome %in% c(opts$case, opts$control) &
  complete.cases(cell_data) & complete.cases(covariate_data)
if (sum(keep) < 3L) stop("Fewer than three complete case/control samples remain.")
phenotype <- phenotype[keep]
fractions <- fractions[keep]
cell_data <- cell_data[keep, , drop = FALSE]
covariate_data <- covariate_data[keep, , drop = FALSE]
outcome <- outcome[keep]
samples <- samples[keep]
if (!all(c(opts$case, opts$control) %in% outcome)) stop("Both case and control samples must remain after filtering.")
if (any(vapply(covariate_data, function(x) length(unique(x)) < 2L, logical(1)))) {
  stop("Every requested covariate must vary in the analyzed samples.")
}

cell_matrix <- as.matrix(cell_data)
storage.mode(cell_matrix) <- "double"
rownames(cell_matrix) <- samples
cell_sums <- rowSums(cell_matrix)
if (any(!is.finite(cell_sums)) || any(abs(cell_sums - 1) > 0.05)) {
  stop("The five liver cell fractions must sum to approximately one per sample.")
}
if (any(vapply(as.data.frame(cell_matrix), function(x) var(x) <= 0, logical(1)))) {
  stop("Each of the five liver cell fractions must vary across samples.")
}

for (method in names(matrices)) matrices[[method]] <- matrices[[method]][, samples, drop = FALSE]
case_status <- as.integer(outcome == opts$case)
covariate_design <- model.matrix(reformulate(covariates), data = covariate_data)
if (nlevels(factor(case_status)) != 2L) stop("The encoded outcome must contain 0 controls and 1 cases.")

run_celldmc <- function(method) {
  fit <- EpiDISH::CellDMC(
    beta.m = matrices[[method]], pheno.v = case_status, frac.m = cell_matrix,
    adjPMethod = "none", adjPThresh = 0, cov.mod = covariate_design,
    sort = FALSE, mc.cores = 1L
  )
  rbindlist(lapply(names(fit$coe), function(cell) {
    result <- as.data.table(fit$coe[[cell]], keep.rownames = "feature_id")
    result[, `:=`(CellType = cell, Method = method)]
    result
  }), use.names = TRUE)
}

ewas <- run_celldmc("EWAS")
chalm <- run_celldmc("CHALM")
camda <- run_celldmc("CAMDA")
ewas[, adjP := p.adjust(p, method = "BH"), by = "CellType"]
ewas[, FDR_family := "EWAS"]
haewas <- rbindlist(list(chalm, camda), use.names = TRUE)
haewas[, adjP := p.adjust(p, method = "BH"), by = "CellType"]
haewas[, FDR_family := "haEWAS_joint"]

phenotype_name <- paste(opts$outcome_column, opts$case, "vs", opts$control, sep = "_")
phenotype_name <- gsub("[^A-Za-z0-9_.-]+", "_", phenotype_name)
analysis_dir <- file.path(opts$output_dir, phenotype_name)
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

ewas_significant <- ewas[adjP < fdr_threshold]
setorder(ewas_significant, CellType, adjP, feature_id)
ha_significant_keys <- unique(haewas[adjP < fdr_threshold, .(feature_id, CellType)])
haewas_significant <- data.table(
  feature_id = character(), CellType = character(), SignificantMethods = character(),
  CHALM_Estimate = numeric(), CHALM_SE = numeric(), CHALM_t = numeric(),
  CHALM_P = numeric(), CHALM_FDR = numeric(),
  CAMDA_Estimate = numeric(), CAMDA_SE = numeric(), CAMDA_t = numeric(),
  CAMDA_P = numeric(), CAMDA_FDR = numeric()
)
if (nrow(ha_significant_keys)) {
  ha_significant <- haewas[ha_significant_keys, on = .(feature_id, CellType)][, .(
    SignificantMethods = paste(Method[adjP < fdr_threshold], collapse = ";"),
    CHALM_Estimate = Estimate[Method == "CHALM"][1L],
    CHALM_SE = SE[Method == "CHALM"][1L],
    CHALM_t = t[Method == "CHALM"][1L],
    CHALM_P = p[Method == "CHALM"][1L],
    CHALM_FDR = adjP[Method == "CHALM"][1L],
    CAMDA_Estimate = Estimate[Method == "CAMDA"][1L],
    CAMDA_SE = SE[Method == "CAMDA"][1L],
    CAMDA_t = t[Method == "CAMDA"][1L],
    CAMDA_P = p[Method == "CAMDA"][1L],
    CAMDA_FDR = adjP[Method == "CAMDA"][1L]
  ), by = .(feature_id, CellType)]
}
setorder(haewas_significant, CellType, feature_id)

fwrite(ewas, file.path(analysis_dir, "cellDMC_EWAS_all.csv"), na = "NA")
fwrite(chalm, file.path(analysis_dir, "cellDMC_haEWAS_CHALM_all.csv"), na = "NA")
fwrite(camda, file.path(analysis_dir, "cellDMC_haEWAS_CAMDA_all.csv"), na = "NA")
fwrite(ewas_significant, file.path(analysis_dir, "cellDMC_significant_EWAS.csv"), na = "NA")
fwrite(haewas_significant, file.path(analysis_dir, "cellDMC_significant_haEWAS.csv"), na = "NA")

audit <- data.table(
  phenotype = phenotype_name,
  outcome_column = opts$outcome_column,
  case = opts$case,
  control = opts$control,
  n_samples = length(samples),
  n_cases = sum(case_status == 1L),
  n_controls = sum(case_status == 0L),
  n_features_EWAS = nrow(matrices$EWAS),
  n_features_CHALM = nrow(matrices$CHALM),
  n_features_CAMDA = nrow(matrices$CAMDA),
  cell_types = paste(cell_columns, collapse = ","),
  covariates = paste(covariates, collapse = ","),
  fdr_threshold = fdr_threshold,
  haewas_fdr = "one BH across CHALM and CAMDA per cell type"
)
fwrite(audit, file.path(analysis_dir, "cellDMC_audit.csv"), na = "NA")
