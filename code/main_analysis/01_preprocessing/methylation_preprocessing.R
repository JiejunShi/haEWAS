# Process beta matrices and write batch-adjusted M- and beta-values.

suppressPackageStartupMessages({
  library(optparse)
  library(ChAMP)
  library(minfi)
  library(impute)
  library(data.table)
  library(sva)
})

option_list <- list(
  make_option("--gse", type = "character", default = NULL),
  make_option("--matrix", type = "character", default = NULL),
  make_option("--idat-root", type = "character", default = Sys.getenv("IDAT_ROOT", "")),
  make_option("--metadata", type = "character", default = NULL),
  make_option("--sample-column", type = "character", default = "Sample"),
  make_option("--batch-column", type = "character", default = NULL),
  make_option("--keep-columns", type = "character", default = ""),
  make_option("--output-dir", type = "character", default = file.path(tempdir(), "methylation_preprocessing")),
  make_option("--prefix", type = "character", default = NULL),
  make_option("--epsilon", type = "double", default = 1e-6)
)
opts <- parse_args(OptionParser(option_list = option_list))
if (xor(is.null(opts$gse), is.null(opts$matrix)) == FALSE) stop("Provide either --gse or --matrix")
if (!is.finite(opts$epsilon) || opts$epsilon <= 0 || opts$epsilon >= 0.5) stop("Invalid epsilon")

normalize_id <- function(x) gsub("\\.", "-", trimws(as.character(x)))

read_beta_matrix <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    beta <- as.matrix(readRDS(path))
  } else {
    input <- fread(path, sep = if (ext == "tsv") "\t" else ",",
                   data.table = FALSE, check.names = FALSE)
    if (ncol(input) < 2L) stop("The beta matrix must contain probe and sample columns")
    rownames(input) <- as.character(input[[1L]])
    beta <- as.matrix(data.frame(lapply(input[, -1L, drop = FALSE], as.numeric),
                                 check.names = FALSE))
    rownames(beta) <- rownames(input)
  }
  storage.mode(beta) <- "double"
  if (is.null(rownames(beta)) || is.null(colnames(beta)) || anyDuplicated(rownames(beta)) ||
      anyDuplicated(colnames(beta))) stop("The beta matrix needs unique probe and sample IDs")
  colnames(beta) <- normalize_id(colnames(beta))
  if (anyDuplicated(colnames(beta))) stop("Duplicate sample IDs after normalization")
  if (any(is.finite(beta) & (beta < 0 | beta > 1))) stop("Beta values must be in [0, 1]")
  beta
}

read_idat_beta <- function(gse_id) {
  idat_root <- opts$idat_root
  if (!nzchar(idat_root)) stop("Set IDAT_ROOT or provide --idat-root")
  base_dir <- if (dir.exists(file.path(idat_root, "All_GEO"))) file.path(idat_root, "All_GEO") else idat_root
  gse_dir <- file.path(base_dir, gse_id)
  if (!dir.exists(gse_dir)) stop("IDAT directory not found for ", gse_id)
  folders <- file.path(gse_dir, c("idat_files", "Idat", "IDAT", "idat"))
  folders <- folders[dir.exists(folders)]
  idat_dir <- if (length(folders)) folders[1L] else gse_dir
  if (!length(list.files(idat_dir, pattern = "idat(\\.gz)?$", recursive = TRUE, ignore.case = TRUE))) {
    stop("No IDAT files found for ", gse_id)
  }

  rg <- read.metharray.exp(idat_dir, recursive = TRUE, force = TRUE)
  sample_names <- sampleNames(rg)
  sample_ids <- normalize_id(sub("_.*$", "", basename(sample_names)))
  if (anyDuplicated(sample_ids)) stop("Duplicated sample IDs in ", gse_id)
  array_name <- as.character(minfi::annotation(rg)[["array"]])
  if (grepl("EPICv2", array_name, ignore.case = TRUE)) stop("ChAMP does not support EPIC v2 in this workflow")
  array_type <- if (grepl("450k", array_name, ignore.case = TRUE)) "450K" else
    if (grepl("EPIC", array_name, ignore.case = TRUE)) "EPIC" else NA_character_
  if (is.na(array_type)) stop("Unsupported array type: ", array_name)

  det_p <- detectionP(rg)
  beta <- getBeta(preprocessRaw(rg), type = "Illumina", offset = 100)
  colnames(beta) <- sample_ids
  colnames(det_p) <- sample_ids
  pd <- data.frame(Sample_Name = sample_ids, Basename = sample_names, row.names = sample_ids)
  filtered <- champ.filter(
    beta = beta, pd = pd, detP = det_p,
    filterDetP = TRUE, detPcut = 0.01,
    filterNoCG = TRUE, filterSNPs = TRUE, filterMultiHit = TRUE, filterXY = TRUE,
    filterBeads = FALSE, ProbeCutoff = 0.20, SampleCutoff = 1,
    autoimpute = FALSE, fixOutlier = FALSE, arraytype = array_type
  )
  beta <- if (is.list(filtered)) filtered$beta else filtered
  if (is.null(beta) || !nrow(beta) || !ncol(beta)) stop("ChAMP filtering returned no data")
  if (!all(colnames(beta) %in% colnames(det_p))) stop("ChAMP returned unexpected sample IDs")
  list(beta = beta, array_type = array_type)
}

if (is.null(opts$matrix)) {
  idat <- read_idat_beta(opts$gse)
  beta <- idat$beta
  prefix <- if (is.null(opts$prefix)) opts$gse else opts$prefix
  array_type <- idat$array_type
} else {
  if (!file.exists(opts$matrix)) stop("Beta matrix not found")
  beta <- read_beta_matrix(opts$matrix)
  prefix <- if (is.null(opts$prefix)) tools::file_path_sans_ext(basename(opts$matrix)) else opts$prefix
  array_type <- NA_character_
}

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
missing_fraction <- colMeans(is.na(beta))
keep_samples <- is.finite(missing_fraction) & missing_fraction <= 0.20
beta <- beta[, keep_samples, drop = FALSE]
message("Retained samples: ", ncol(beta))
if (ncol(beta) < 2L) stop("Fewer than two samples remain after the 20% missing-value filter")
beta <- beta[rowSums(!is.na(beta)) > 0L, , drop = FALSE]
if (!nrow(beta)) stop("No probes remain after removing all-missing rows")
if (any(!is.na(beta) & !is.finite(beta))) stop("The beta matrix contains non-finite values")
if (any(beta < 0 | beta > 1, na.rm = TRUE)) stop("Beta values must be in [0, 1]")
if (anyNA(beta)) {
  set.seed(1)
  beta <- impute.knn(as.matrix(beta), k = 5)$data
}
if (anyNA(beta)) stop("kNN imputation left missing values")

if (!is.na(array_type)) {
  normalized <- champ.norm(beta = beta, method = "BMIQ", arraytype = array_type,
                           cores = 2, resultsDir = opts$output_dir)
  beta <- if (is.list(normalized)) normalized$beta else normalized
  if (is.null(beta) || !nrow(beta) || !ncol(beta) || anyNA(beta)) stop("BMIQ returned invalid data")
  if (any(!is.finite(beta)) || any(beta < 0 | beta > 1)) stop("BMIQ returned values outside [0, 1]")
}

beta_for_m <- pmin(pmax(beta, opts$epsilon), 1 - opts$epsilon)
mvalue <- log2(beta_for_m / (1 - beta_for_m))
if (!is.null(opts$batch_column)) {
  if (is.null(opts$metadata) || !file.exists(opts$metadata)) stop("A metadata file is required for ComBat")
  keep_columns <- trimws(strsplit(opts$keep_columns, ",", fixed = TRUE)[[1L]])
  keep_columns <- keep_columns[nzchar(keep_columns)]
  if (!length(keep_columns)) stop("Specify biological covariates with --keep-columns")
  if (opts$batch_column %in% keep_columns) stop("Do not include the batch column in --keep-columns")

  metadata <- fread(opts$metadata, sep = if (tolower(tools::file_ext(opts$metadata)) == "tsv") "\t" else ",",
                    data.table = FALSE, check.names = FALSE)
  required <- c(opts$sample_column, opts$batch_column, keep_columns)
  if (!all(required %in% names(metadata))) stop("A requested metadata column is missing")
  metadata[[opts$sample_column]] <- normalize_id(metadata[[opts$sample_column]])
  if (anyDuplicated(metadata[[opts$sample_column]])) stop("Duplicate sample IDs in metadata")
  metadata <- metadata[match(colnames(beta_for_m), metadata[[opts$sample_column]]), , drop = FALSE]
  if (anyNA(metadata[[opts$sample_column]])) stop("Metadata does not cover every beta sample")
  if (anyNA(metadata[, c(opts$batch_column, keep_columns), drop = FALSE])) stop("Batch/covariate metadata contains missing values")

  batch <- droplevels(factor(metadata[[opts$batch_column]]))
  if (nlevels(batch) > 1L) {
    if (any(table(batch) < 2L)) stop("Each batch needs at least two samples for ComBat")
    covariates <- metadata[, keep_columns, drop = FALSE]
    covariates[] <- lapply(covariates, function(x) if (is.character(x)) factor(x) else x)
    mod <- model.matrix(reformulate(keep_columns), data = covariates)
    mvalue <- ComBat(dat = mvalue, batch = batch, mod = mod,
                     par.prior = TRUE, prior.plots = FALSE)
  } else {
    message("Only one batch level; writing unadjusted M-values")
  }
} else {
  message("No batch column supplied; writing unadjusted M-values")
}

mvalue_table <- as.data.table(mvalue)
setnames(mvalue_table, colnames(mvalue))
mvalue_table[, Probe := rownames(mvalue)]
setcolorder(mvalue_table, c("Probe", colnames(mvalue)))
mvalue_output <- file.path(opts$output_dir, paste0(prefix, "_mvalue.tsv"))
fwrite(mvalue_table, mvalue_output, sep = "\t", quote = FALSE)

beta_from_mvalue <- stats::plogis(mvalue * log(2))
beta_from_mvalue_table <- as.data.table(beta_from_mvalue)
setnames(beta_from_mvalue_table, colnames(beta))
beta_from_mvalue_table[, Probe := rownames(beta_from_mvalue)]
setcolorder(beta_from_mvalue_table, c("Probe", colnames(beta)))
assay_output <- file.path(opts$output_dir, paste0(prefix, "_assay.tsv"))
fwrite(beta_from_mvalue_table, assay_output, sep = "\t", quote = FALSE)

message("Processed M-value matrix: ", mvalue_output)
message("Inverse-transformed beta assay: ", assay_output)
