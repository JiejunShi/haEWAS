#!/usr/bin/env Rscript
# Fit method-specific eQTM associations from matched expression and methylation data.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(limma)
})

option_list <- list(
  make_option("--clinical-file", type = "character"),
  make_option("--sample-covariates-file", type = "character"),
  make_option("--expression-matrix", type = "character"),
  make_option("--beta-matrix", type = "character"),
  make_option("--chalm-matrix", type = "character"),
  make_option("--camda-matrix", type = "character"),
  make_option("--output-dir", type = "character"),
  make_option("--pair-file", type = "character", default = ""),
  make_option("--clinical-id-column", type = "character", default = "sample_id"),
  make_option("--covariate-id-column", type = "character", default = "sample_id"),
  make_option("--expression-id-column", type = "character", default = "expression_sample_id"),
  make_option("--methylation-id-column", type = "character", default = "methylation_sample_id"),
  make_option("--age-column", type = "character", default = "age"),
  make_option("--sex-column", type = "character", default = "sex"),
  make_option("--trait-column", type = "character", default = "trait"),
  make_option("--cell-column", type = "character", default = "cell_proportion"),
  make_option("--batch-column", type = "character", default = "batch"),
  make_option("--array-column", type = "character", default = "array_id"),
  make_option("--marker-genes", type = "character"),
  make_option("--pair-window", type = "integer", default = 250000L),
  make_option("--epsilon", type = "double", default = 1e-6),
  make_option("--fdr-cutoff", type = "double", default = 0.05)
)

opts <- parse_args(OptionParser(option_list = option_list))
required <- c(
  "clinical_file", "sample_covariates_file", "expression_matrix", "beta_matrix",
  "chalm_matrix", "camda_matrix", "output_dir", "marker_genes"
)
missing <- required[vapply(required, function(x) is.null(opts[[x]]) || !nzchar(opts[[x]]), logical(1))]
if (length(missing)) stop("Missing required options: ", paste(missing, collapse = ", "))
if (opts$pair_window < 0L || opts$epsilon <= 0 || opts$epsilon >= 0.5 ||
    opts$fdr_cutoff <= 0 || opts$fdr_cutoff >= 1) {
  stop("Invalid pair window, epsilon or FDR cutoff.")
}

input_paths <- unlist(opts[c(
  "clinical_file", "sample_covariates_file", "expression_matrix",
  "beta_matrix", "chalm_matrix", "camda_matrix"
)], use.names = FALSE)
missing_files <- input_paths[!file.exists(input_paths)]
if (length(missing_files)) stop("Input file not found: ", paste(missing_files, collapse = ", "))

output_dir <- opts$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
pair_file <- if (nzchar(opts$pair_file)) opts$pair_file else file.path(output_dir, "cpg_gene_pairs.csv")

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
normalize_id <- function(x) gsub("\\.", "-", trimws(as.character(x)))

read_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(as.data.table(readRDS(path)))
  fread(path, sep = if (ext == "tsv") "\t" else ",", data.table = TRUE,
        check.names = FALSE, showProgress = FALSE)
}

read_matrix <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") {
    matrix <- as.matrix(readRDS(path))
  } else {
    input <- fread(path, sep = if (ext == "tsv") "\t" else ",",
                   data.table = FALSE, check.names = FALSE, showProgress = FALSE)
    if (ncol(input) < 2L) stop("Matrix must have an ID column and sample columns: ", path)
    rownames(input) <- trimws(as.character(input[[1L]]))
    matrix <- as.matrix(input[, -1L, drop = FALSE])
  }
  storage.mode(matrix) <- "double"
  if (is.null(rownames(matrix)) || is.null(colnames(matrix))) stop("Matrix lacks dimnames: ", path)
  rownames(matrix) <- trimws(as.character(rownames(matrix)))
  colnames(matrix) <- normalize_id(colnames(matrix))
  if (anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix))) {
    stop("Duplicate feature or sample IDs in matrix: ", path)
  }
  matrix
}

read_feature_ids <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(rownames(readRDS(path)))
  ids <- fread(path, sep = if (ext == "tsv") "\t" else ",", select = 1L,
               data.table = FALSE, showProgress = FALSE)[[1L]]
  trimws(as.character(ids))
}

required_columns <- function(x, columns, label) {
  absent <- setdiff(columns, names(x))
  if (length(absent)) stop(label, " lacks columns: ", paste(absent, collapse = ", "))
}

clinical_raw <- read_table(opts$clinical_file)
covariates_raw <- read_table(opts$sample_covariates_file)
clinical_columns <- c(
  opts$clinical_id_column, opts$expression_id_column, opts$methylation_id_column,
  opts$age_column, opts$sex_column, opts$trait_column
)
covariate_columns <- c(
  opts$covariate_id_column, opts$cell_column, opts$batch_column, opts$array_column
)
required_columns(clinical_raw, clinical_columns, "Clinical table")
required_columns(covariates_raw, covariate_columns, "Sample-covariate table")

clinical <- data.table(
  subject_id = normalize_id(clinical_raw[[opts$clinical_id_column]]),
  expression_id = normalize_id(clinical_raw[[opts$expression_id_column]]),
  methylation_id = normalize_id(clinical_raw[[opts$methylation_id_column]]),
  age = to_num(clinical_raw[[opts$age_column]]),
  sex = as.factor(clinical_raw[[opts$sex_column]]),
  trait = to_num(clinical_raw[[opts$trait_column]])
)
covariates <- data.table(
  subject_id = normalize_id(covariates_raw[[opts$covariate_id_column]]),
  cell_proportion = to_num(covariates_raw[[opts$cell_column]]),
  batch = as.factor(covariates_raw[[opts$batch_column]]),
  array_id = as.factor(covariates_raw[[opts$array_column]])
)
if (anyDuplicated(clinical$subject_id) || anyDuplicated(covariates$subject_id)) {
  stop("Sample IDs must be unique in both input metadata tables.")
}
pheno <- merge(clinical, covariates, by = "subject_id", all = FALSE, sort = FALSE)
if (!nrow(pheno)) stop("The two metadata tables have no matching subject IDs.")
if (anyNA(pheno[, .(age, sex, cell_proportion, batch, array_id)])) {
  stop("Required model covariates contain missing values; resolve them before analysis.")
}

rna_matrix <- read_matrix(opts$expression_matrix)
rownames(rna_matrix) <- sub("\\..*$", "", rownames(rna_matrix))
if (anyDuplicated(rownames(rna_matrix))) stop("Duplicate expression gene IDs after version removal.")
marker_genes <- trimws(strsplit(opts$marker_genes, ",", fixed = TRUE)[[1L]])
marker_genes <- marker_genes[nzchar(marker_genes)]
if (!length(marker_genes) || any(!marker_genes %in% rownames(rna_matrix))) {
  stop("One or more --marker-genes are absent from the expression matrix.")
}

make_pairs <- function(beta_path, rna_ids, pair_window, path) {
  if (nzchar(opts$pair_file)) {
    if (!file.exists(path)) stop("--pair-file not found: ", path)
    pairs <- fread(path, data.table = TRUE, check.names = FALSE)
    required_columns(pairs, c("cpg", "gene"), "Pair table")
    if (!"gene_symbol" %in% names(pairs)) pairs[, gene_symbol := gene]
    return(unique(pairs[, .(cpg = as.character(cpg), gene = as.character(gene),
                            gene_symbol = as.character(gene_symbol))]))
  }
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(GenomeInfoDb)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(EnsDb.Hsapiens.v75)
  })
  cpg_ids <- trimws(as.character(read_feature_ids(beta_path)))
  gene_ids <- unique(sub("\\..*$", "", as.character(rna_ids)))
  data("Locations", package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  cpg_locations <- Locations[intersect(cpg_ids, rownames(Locations)), , drop = FALSE]
  cpg_ranges <- GRanges(
    seqnames = cpg_locations$chr,
    ranges = IRanges(start = cpg_locations$pos, end = cpg_locations$pos),
    cpg_id = rownames(cpg_locations)
  )
  gene_ranges <- genes(EnsDb.Hsapiens.v75)
  seqlevelsStyle(gene_ranges) <- "UCSC"
  gene_ranges <- gene_ranges[names(gene_ranges) %in% gene_ids]
  gene_ranges <- keepSeqlevels(
    gene_ranges, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse"
  )
  gene_id <- names(gene_ranges)
  gene_symbol <- if ("gene_name" %in% names(mcols(gene_ranges))) {
    as.character(mcols(gene_ranges)$gene_name)
  } else {
    gene_id
  }
  missing_symbol <- is.na(gene_symbol) | !nzchar(gene_symbol)
  gene_symbol[missing_symbol] <- gene_id[missing_symbol]
  start(gene_ranges) <- pmax(1L, start(gene_ranges) - pair_window)
  end(gene_ranges) <- end(gene_ranges) + pair_window
  hits <- findOverlaps(cpg_ranges, gene_ranges)
  pairs <- unique(data.table(
    cpg = mcols(cpg_ranges)$cpg_id[queryHits(hits)],
    gene = gene_id[subjectHits(hits)],
    gene_symbol = gene_symbol[subjectHits(hits)]
  ))
  fwrite(pairs, path)
  pairs
}

pairs <- make_pairs(opts$beta_matrix, rownames(rna_matrix), opts$pair_window, pair_file)
if (!nrow(pairs)) stop("No CpG-gene pairs were found.")
paired_cpgs <- unique(pairs$cpg)
paired_genes <- unique(pairs$gene)
marker_columns <- paste0("marker_", seq_along(marker_genes))
marker_expression <- t(log2(rna_matrix[marker_genes, , drop = FALSE] + 1))
colnames(marker_expression) <- marker_columns
rna_log <- log2(rna_matrix + 1)
results_all <- list()
audit <- list()

method_files <- list(
  BETA = opts$beta_matrix,
  CHALM = opts$chalm_matrix,
  CAMDA = opts$camda_matrix
)
for (method in names(method_files)) {
  message("[RUN] ", method)
  methylation <- read_matrix(method_files[[method]])
  samples <- pheno[
    expression_id %in% colnames(rna_matrix) & methylation_id %in% colnames(methylation)
  ]
  if (nrow(samples) < 4L) stop(method, ": fewer than four matched samples.")
  rna <- rna_log[, samples$expression_id, drop = FALSE]
  meth <- methylation[, samples$methylation_id, drop = FALSE]
  colnames(rna) <- samples$subject_id
  colnames(meth) <- samples$subject_id

  expression_design_data <- data.frame(
    age = samples$age, sex = samples$sex,
    marker_expression[
      match(samples$expression_id, rownames(marker_expression)), , drop = FALSE
    ],
    check.names = FALSE
  )
  expression_design <- model.matrix(
    reformulate(c("age", "sex", marker_columns)), data = expression_design_data
  )
  expression_residuals <- residuals(lmFit(rna, expression_design), rna)

  meth[!is.na(meth) & meth < opts$epsilon] <- opts$epsilon
  meth[!is.na(meth) & meth > (1 - opts$epsilon)] <- 1 - opts$epsilon
  m_values <- log2(meth / (1 - meth))
  methylation_design <- model.matrix(
    ~ cell_proportion + batch + array_id + age + sex, data = samples
  )
  methylation_residuals <- residuals(lmFit(m_values, methylation_design), m_values)

  pairs_valid <- pairs[
    cpg %in% rownames(methylation_residuals) & gene %in% rownames(expression_residuals)
  ]
  if (!nrow(pairs_valid)) {
    warning(method, ": no CpG-gene pairs overlap the matrices.")
    next
  }
  valid_samples <- is.finite(samples$trait)
  if (sum(valid_samples) < 4L) stop(method, ": fewer than four samples have a valid trait value.")
  y <- expression_residuals[pairs_valid$gene, valid_samples, drop = FALSE]
  x <- methylation_residuals[pairs_valid$cpg, valid_samples, drop = FALSE]
  trait_design <- model.matrix(~ trait, data = samples[valid_samples])
  y <- residuals(lmFit(y, trait_design), y)
  x <- residuals(lmFit(x, trait_design), x)

  y_centered <- y - rowMeans(y, na.rm = TRUE)
  x_centered <- x - rowMeans(x, na.rm = TRUE)
  n_valid <- rowSums(!is.na(y_centered) & !is.na(x_centered))
  sd_y <- sqrt(rowSums(y_centered^2, na.rm = TRUE) / (n_valid - 1))
  sd_x <- sqrt(rowSums(x_centered^2, na.rm = TRUE) / (n_valid - 1))
  cov_xy <- rowSums(x_centered * y_centered, na.rm = TRUE) / (n_valid - 1)
  valid <- is.finite(sd_y) & is.finite(sd_x) & sd_y > 0 & sd_x > 0 & n_valid > 3L

  correlation <- rep(NA_real_, nrow(pairs_valid))
  estimate <- se <- pvalue <- rep(NA_real_, nrow(pairs_valid))
  correlation[valid] <- cov_xy[valid] / (sd_y[valid] * sd_x[valid])
  correlation <- pmin(pmax(correlation, -0.999999), 0.999999)
  df <- n_valid - 3L
  statistic <- correlation * sqrt(df / (1 - correlation^2))
  pvalue[valid] <- 2 * pt(-abs(statistic[valid]), df = df[valid])
  estimate[valid] <- correlation[valid] * (sd_y[valid] / sd_x[valid])
  se[valid] <- (sd_y[valid] / sd_x[valid]) *
    sqrt((1 - correlation[valid]^2) / df[valid])

  result <- data.table(
    Method = method,
    cpg = pairs_valid$cpg,
    gene = pairs_valid$gene,
    gene_symbol = pairs_valid$gene_symbol,
    estimate = estimate,
    se = se,
    pvalue = pvalue
  )
  result <- result[is.finite(pvalue)]
  result[, fdr := p.adjust(pvalue, method = "BH")]
  setorder(result, pvalue)
  fwrite(result, file.path(output_dir, paste0("eqtm_results_", tolower(method), ".csv")))
  results_all[[method]] <- result
  audit[[method]] <- data.table(
    Method = method,
    matched_samples = nrow(samples),
    trait_complete_samples = sum(valid_samples),
    tested_pairs = nrow(pairs_valid),
    valid_results = nrow(result)
  )
}

if (!length(results_all)) stop("No method produced eQTM results.")
fwrite(rbindlist(results_all, use.names = TRUE, fill = TRUE),
       file.path(output_dir, "eqtm_results_merged_all.csv"))
fwrite(rbindlist(audit, use.names = TRUE, fill = TRUE),
       file.path(output_dir, "eqtm_analysis_audit.csv"))
message("[DONE] eQTM results written to: ", output_dir)
