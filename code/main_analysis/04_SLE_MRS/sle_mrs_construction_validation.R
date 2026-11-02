# Train and evaluate MRS prediction models.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(EpiDISH)
  library(pROC)
})

option_list <- list(
  make_option("--meta_dir", type = "character", help = "Two-cohort meta-analysis output directory"),
  make_option("--output_dir", type = "character", help = "Model and external-validation output directory"),
  make_option("--pheno_file", type = "character", help = "Merged GEO phenotype CSV"),
  make_option("--beta_dir", type = "character", help = "Beta matrix directory"),
  make_option("--chalm_dir", type = "character", help = "CHALM matrix directory"),
  make_option("--camda_dir", type = "character", help = "CAMDA matrix directory"),
  make_option("--train_cohort", type = "character", default = "GSE59250"),
  make_option("--tuning_cohorts", type = "character", default = "GSE59250,GSE82218"),
  make_option("--test_cohorts", type = "character",
              default = "GSE110607,GSE118144,GSE96879"),
  make_option("--phenotype", type = "character"),
  make_option("--control_label", type = "character"),
  make_option("--tissue", type = "character"),
  make_option("--fdr_cutoff", type = "numeric", default = 0.05),
  make_option("--candidate_p_cutoff", type = "numeric", default = 0.05),
  make_option("--top_min", type = "integer", default = 5L),
  make_option("--top_max", type = "integer", default = 150L),
  make_option("--top_step", type = "integer", default = 5L),
  make_option("--bootstrap_reps", type = "integer", default = 500L),
  make_option("--bootstrap_seed", type = "integer", default = 20260715L),
  make_option("--epsilon", type = "numeric", default = 1e-6)
)

opts <- parse_args(OptionParser(option_list = option_list))
required_opts <- c("meta_dir", "output_dir", "pheno_file", "beta_dir", "chalm_dir", "camda_dir",
                   "train_cohort", "tuning_cohorts", "test_cohorts",
                   "phenotype", "control_label", "tissue")
missing_opts <- required_opts[vapply(required_opts, function(x) is.null(opts[[x]]) || !nzchar(opts[[x]]), logical(1))]
if (length(missing_opts) > 0L) stop("Missing required option(s): ", paste(missing_opts, collapse = ", "))
if (!is.finite(opts$fdr_cutoff) || opts$fdr_cutoff <= 0 || opts$fdr_cutoff >= 1) stop("--fdr_cutoff must be in (0, 1).")
if (!is.finite(opts$candidate_p_cutoff) || opts$candidate_p_cutoff <= 0 || opts$candidate_p_cutoff >= 1) {
  stop("--candidate_p_cutoff must be in (0, 1).")
}
if (opts$top_min < 1L || opts$top_max < opts$top_min || opts$top_step < 1L) {
  stop("Invalid Top-N range.")
}
if (opts$bootstrap_reps < 1L) stop("--bootstrap_reps must be positive.")

TRAIN_COHORT <- opts$train_cohort
TUNING_COHORTS <- unique(trimws(unlist(strsplit(opts$tuning_cohorts, ",", fixed = TRUE))))
TUNING_COHORTS <- TUNING_COHORTS[nzchar(TUNING_COHORTS)]
if (!TRAIN_COHORT %in% TUNING_COHORTS) TUNING_COHORTS <- c(TRAIN_COHORT, TUNING_COHORTS)
TEST_COHORTS <- unique(trimws(unlist(strsplit(opts$test_cohorts, ",", fixed = TRUE))))
TEST_COHORTS <- TEST_COHORTS[nzchar(TEST_COHORTS)]
if (length(intersect(TUNING_COHORTS, TEST_COHORTS)) > 0L) {
  stop("Tuning and external-test cohorts must be disjoint: ", paste(intersect(TUNING_COHORTS, TEST_COHORTS), collapse = ", "))
}
if (TRAIN_COHORT != "GSE59250" ||
    !setequal(TUNING_COHORTS, c("GSE59250", "GSE82218")) ||
    !setequal(TEST_COHORTS, c("GSE110607", "GSE118144", "GSE96879"))) {
  stop("Cohort assignments do not match the SLE MRS method.")
}
SCORE_COHORTS <- unique(c(TUNING_COHORTS, TEST_COHORTS))
METHODS <- c("beta", "chalm", "camda")
TOP_N <- seq.int(opts$top_min, opts$top_max, by = opts$top_step)
EXPECTED_COUNTS <- list(
  GSE59250 = c(228L, 206L),
  GSE82218 = c(30L, 25L),
  GSE110607 = c(52L, 52L),
  GSE118144 = c(80L, 65L),
  GSE96879 = c(57L, 33L)
)

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
message("[INFO] Training cohort: ", TRAIN_COHORT)
message("[INFO] Development tuning cohorts: ", paste(TUNING_COHORTS, collapse = ", "))
message("[INFO] External cohorts: ", paste(TEST_COHORTS, collapse = ", "))
message("[INFO] Meta FDR cutoff: ", opts$fdr_cutoff)
message("[INFO] Candidate meta P cutoff: ", opts$candidate_p_cutoff)
message("[INFO] Top-N grid: ", paste(TOP_N, collapse = ","))
message("[INFO] Top-N is selected in development data by bootstrap one-SE rule; external labels are not used.")

first_existing <- function(nms, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0L) return(hit[[1]])
  if (required) {
    stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(nms, collapse = ", "))
  }
  NA_character_
}

normalize_sample <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("_.*$", "", x)
  x <- sub("^[^.]+\\.(GSM[0-9]+.*)$", "\\1", x)
  gsub("\\.", "-", x)
}

normalize_sex <- function(x) {
  low <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(low))
  out[low %in% c("f", "female", "woman", "0")] <- 0
  out[low %in% c("m", "male", "man", "1")] <- 1
  out
}

safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))

to_m_value <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat[is.infinite(mat)] <- NA_real_
  mat[is.finite(mat) & mat < 0] <- 0
  mat[is.finite(mat) & mat > 1] <- 1
  mat[is.finite(mat) & mat < opts$epsilon] <- opts$epsilon
  mat[is.finite(mat) & mat > 1 - opts$epsilon] <- 1 - opts$epsilon
  log2(mat / (1 - mat))
}

find_matrix_file <- function(cohort, method, base_dir) {
  exact <- switch(
    method,
    beta = c(sprintf("%s_assay.tsv", cohort), sprintf("%s_assay.tsv.gz", cohort)),
    chalm = c(sprintf("%s.chalm.filtered.tsv", cohort), sprintf("%s.chalm.filtered.tsv.gz", cohort)),
    camda = c(sprintf("%s.camda.filtered.tsv", cohort), sprintf("%s.camda.filtered.tsv.gz", cohort))
  )
  hits <- file.path(base_dir, exact)
  hits <- hits[file.exists(hits) & file.info(hits)$size > 0]
  if (length(hits) > 0L) return(hits[[1]])
  stop("Matrix file not found for ", cohort, " / ", method, " in ", base_dir)
}

read_matrix <- function(path, requested_samples, requested_cpgs = NULL) {
  sep <- if (grepl("\\.(tsv|txt)(\\.gz)?$", path, ignore.case = TRUE)) "\t" else ","
  header <- names(fread(path, sep = sep, nrows = 0L, check.names = FALSE))
  if (length(header) < 2L) stop("Matrix has fewer than two columns: ", path)
  normalized <- normalize_sample(header[-1])
  sample_map <- setNames(header[-1], normalized)
  keep_samples <- intersect(requested_samples, names(sample_map))
  if (length(keep_samples) == 0L) stop("No phenotype samples matched matrix: ", path)

  selected <- c(header[[1]], unname(sample_map[keep_samples]))
  dt <- fread(path, sep = sep, select = selected, check.names = FALSE, data.table = TRUE)
  cpg <- trimws(as.character(dt[[1]]))
  keep <- !is.na(cpg) & nzchar(cpg)
  if (!is.null(requested_cpgs)) keep <- keep & cpg %in% requested_cpgs
  dt <- dt[keep]
  cpg <- cpg[keep]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- cpg
  colnames(mat) <- keep_samples
  mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  mat
}

find_meta_file <- function(method) {
  source_method <- if (method == "beta") "m" else method
  current <- file.path(
    opts$meta_dir,
    paste0("meta_effectsize_", source_method, "_", opts$tissue, ".csv")
  )
  if (file.exists(current) && file.info(current)$size > 0) return(current)
  patterns <- c(
    sprintf("^meta_effectsize_.*_%s_raw_filtered\\.csv$", method),
    sprintf("^meta_effectsize_.*_%s_raw\\.csv$", method),
    sprintf("^meta_effectsize_.*_%s_%s\\.csv$", source_method, opts$tissue)
  )
  for (pattern in patterns) {
    hits <- list.files(opts$meta_dir, pattern = pattern, full.names = TRUE)
    hits <- hits[file.info(hits)$size > 0]
    if (length(hits) == 1L) return(hits[[1]])
    if (length(hits) > 1L) stop("Multiple meta files matched for ", method, ": ", paste(hits, collapse = "; "))
  }
  stop("Two-cohort meta file not found for method: ", method)
}

read_meta <- function(method) {
  path <- find_meta_file(method)
  dt <- fread(path, check.names = FALSE, data.table = TRUE)
  if ("cohorts" %in% names(dt)) {
    source_cohorts <- unique(trimws(unlist(strsplit(as.character(dt$cohorts), ",", fixed = TRUE))))
    source_cohorts <- source_cohorts[nzchar(source_cohorts)]
    if (!setequal(source_cohorts, c("GSE59250", "GSE82218"))) {
      stop("Meta file is not based on the two discovery cohorts: ", path)
    }
  } else {
    message("[WARN] Meta file does not record cohort provenance: ", basename(path))
  }
  cpg_col <- first_existing(names(dt), c("cpg", "CpG", "Name", "probe"), "CpG column")
  effect_col <- first_existing(
    names(dt), c("meta_Estimate", "meta_EffectSize", "Estimate_fixed", "Estimate"),
    "meta effect column"
  )
  p_col <- first_existing(names(dt), c("pval_final", "pval", "pValue", "p_value", "P.Value"), "meta P-value column")
  fdr_col <- first_existing(names(dt), c("fdr_all", "fdr", "fdr_U", "FDR", "adj.P.Val"),
                            "meta FDR column", required = FALSE)
  source_fdr <- if (!is.na(fdr_col)) safe_numeric(dt[[fdr_col]]) else rep(NA_real_, nrow(dt))
  out <- data.table(
    cpg = trimws(as.character(dt[[cpg_col]])),
    value_method = method,
    meta_weight = safe_numeric(dt[[effect_col]]),
    meta_p = safe_numeric(dt[[p_col]]),
    meta_fdr_method_specific = source_fdr,
    effect_column = effect_col,
    fdr_column = fdr_col,
    source_file = path
  )
  out <- out[!is.na(cpg) & nzchar(cpg) & is.finite(meta_weight) &
               is.finite(meta_p) & meta_p >= 0 & meta_p <= 1]
  out[, candidate := meta_p < opts$candidate_p_cutoff]
  setorder(out, meta_p, cpg)
  out <- unique(out, by = "cpg")
  message(
    "[META] ", method, ": tested=", nrow(out),
    ", P<", opts$candidate_p_cutoff, "=", sum(out$candidate)
  )
  out
}

meta_all <- setNames(lapply(METHODS, read_meta), METHODS)

beta_fdr <- meta_all$beta$meta_fdr_method_specific
missing_beta_fdr <- !is.finite(beta_fdr)
if (any(missing_beta_fdr)) {
  beta_fdr[missing_beta_fdr] <- p.adjust(meta_all$beta$meta_p, method = "BH")[missing_beta_fdr]
}
meta_all$beta[, `:=`(
  meta_fdr = beta_fdr,
  fdr_family = "EWAS_method_specific",
  significant = is.finite(meta_fdr) & meta_fdr < opts$fdr_cutoff
)]

ha_pooled <- rbindlist(lapply(c("chalm", "camda"), function(method) {
  x <- copy(meta_all[[method]][, .(cpg, meta_p)])
  x[, value_method := method]
  x
}), use.names = TRUE)
ha_pooled[, pool_id := paste(value_method, cpg, sep = "\r")]
if (anyDuplicated(ha_pooled$pool_id)) stop("Duplicate CpG-method tests in CHALM/CAMDA meta results")
ha_pooled[, meta_fdr_pooled := p.adjust(meta_p, method = "BH")]
ha_pooled[, significant_pooled := meta_fdr_pooled < opts$fdr_cutoff]
fwrite(ha_pooled[, .(value_method, cpg, meta_p, meta_fdr_pooled, significant_pooled)],
       file.path(opts$output_dir, "haewas_pooled_bh.csv"))

for (method in c("chalm", "camda")) {
  x <- meta_all[[method]]
  pool_index <- match(paste(method, x$cpg, sep = "\r"), ha_pooled$pool_id)
  if (anyNA(pool_index)) stop("Could not match pooled haEWAS FDR rows for ", method)
  x[, meta_fdr_pooled := ha_pooled$meta_fdr_pooled[pool_index]]
  x[, `:=`(
    meta_fdr = meta_fdr_pooled,
    fdr_family = "pooled_CHALM_CAMDA",
    significant = is.finite(meta_fdr_pooled) & meta_fdr_pooled < opts$fdr_cutoff
  )]
  meta_all[[method]] <- x
}

message("[POOLED BH] CHALM and CAMDA combined tests=", nrow(ha_pooled),
        "; significant CHALM=", sum(meta_all$chalm$significant),
        "; significant CAMDA=", sum(meta_all$camda$significant))
meta_sig <- lapply(meta_all, function(x) x[significant == TRUE])
meta_candidates <- lapply(meta_all, function(x) x[candidate == TRUE])
if (nrow(meta_candidates$beta) < opts$top_min) stop("Too few beta candidates pass the meta P-value cutoff.")
if (nrow(meta_candidates$chalm) + nrow(meta_candidates$camda) < opts$top_min) {
  stop("Too few CHALM/CAMDA candidates pass the meta P-value cutoff.")
}

sig_long <- rbindlist(meta_sig, fill = TRUE)
sig_long[, feature_group := ifelse(value_method == "beta", "EWAS", "haEWAS")]
sig_long[, feature_id := paste(cpg, toupper(value_method), sep = "-")]
sig_long[, ha_source := fifelse(
  feature_group == "haEWAS" & cpg %in% intersect(meta_sig$chalm$cpg, meta_sig$camda$cpg),
  "CHALM+CAMDA", fifelse(feature_group == "haEWAS", toupper(value_method), NA_character_)
)]
fwrite(sig_long, file.path(opts$output_dir, "meta_significant_features_long.csv"))

candidate_long <- rbindlist(meta_candidates, fill = TRUE)
candidate_long[, feature_group := ifelse(value_method == "beta", "EWAS", "haEWAS")]
candidate_long[, feature_id := paste(cpg, toupper(value_method), sep = "-")]
setorder(candidate_long, feature_group, meta_p, feature_id)
fwrite(candidate_long, file.path(opts$output_dir, "meta_candidates_long.csv"))

union_cpgs <- sort(unique(sig_long$cpg))
weight_long <- rbindlist(lapply(METHODS, function(method) {
  x <- meta_all[[method]][cpg %in% union_cpgs]
  x[, significant_in_method := significant]
  x
}), fill = TRUE)
weight_wide <- dcast(
  weight_long,
  cpg ~ value_method,
  value.var = c("meta_weight", "meta_p", "meta_fdr", "significant_in_method")
)
fwrite(weight_long, file.path(opts$output_dir, "meta_effects_long.csv"))
fwrite(weight_wide, file.path(opts$output_dir, "meta_effects_wide.csv"))

pheno_all <- fread(opts$pheno_file, check.names = FALSE, data.table = TRUE)
required_pheno <- c("GSE_ID", "GSM_ID", "Disease_State2", "Tissue2")
missing_pheno <- setdiff(required_pheno, names(pheno_all))
if (length(missing_pheno) > 0L) stop("Phenotype file lacks: ", paste(missing_pheno, collapse = ", "))
age_col <- first_existing(names(pheno_all), c("Age", "age", "age_years", "age_at_sampling"), "age column")
sex_col <- first_existing(names(pheno_all), c("Sex", "sex", "gender", "Gender"), "sex column")

read_pheno <- function(cohort) {
  dt <- pheno_all[
    trimws(as.character(GSE_ID)) == cohort &
      tolower(trimws(as.character(Tissue2))) == tolower(opts$tissue) &
      trimws(as.character(Disease_State2)) %in% c(opts$phenotype, opts$control_label)
  ]
  dt[, GSM_ID := normalize_sample(GSM_ID)]
  dt[, outcome := as.integer(trimws(as.character(Disease_State2)) == opts$phenotype)]
  dt[, Age_observed := if (!is.na(age_col)) safe_numeric(get(age_col)) else NA_real_]
  dt[, Age_model := Age_observed]
  dt[, Sex_model := if (!is.na(sex_col)) normalize_sex(get(sex_col)) else NA_real_]
  dt <- unique(dt[, .(GSM_ID, GSE_ID, Disease_State2, outcome, Age_observed, Age_model, Sex_model)], by = "GSM_ID")
  if (nrow(dt) < 10L || uniqueN(dt$outcome) != 2L) stop("Cohort lacks usable case/control blood samples: ", cohort)
  observed_counts <- c(sum(dt$outcome == 1L), sum(dt$outcome == 0L))
  if (!identical(as.integer(observed_counts), EXPECTED_COUNTS[[cohort]])) {
    stop(cohort, ": observed case/control counts ",
         paste(observed_counts, collapse = "/"), " do not match Method counts ",
         paste(EXPECTED_COUNTS[[cohort]], collapse = "/"))
  }
  if (any(!is.finite(dt$Age_model)) || any(!is.finite(dt$Sex_model))) {
    stop(cohort, ": missing age or sex values; imputation is not performed.")
  }
  message("[PHENO] ", cohort, ": N=", nrow(dt), ", case=", sum(dt$outcome == 1L), ", control=", sum(dt$outcome == 0L))
  dt
}

estimate_cells <- function(beta, cohort) {
  data("centDHSbloodDMC.m", package = "EpiDISH", envir = environment())
  ref <- get("centDHSbloodDMC.m", envir = environment())
  common <- intersect(rownames(beta), rownames(ref))
  if (length(common) < 50L) stop("Too few EpiDISH reference CpGs for ", cohort, ": ", length(common))
  beta_sub <- beta[common, , drop = FALSE]
  if (any(!is.finite(beta_sub))) stop(cohort, ": beta matrix has missing values before EpiDISH estimation.")
  beta_sub <- pmin(pmax(beta_sub, 0), 1)
  est <- EpiDISH::epidish(beta.m = beta_sub, ref.m = ref[common, , drop = FALSE], method = "RPC")$estF
  cells <- as.data.table(est, keep.rownames = "GSM_ID")
  cells[, GSM_ID := normalize_sample(GSM_ID)]
  old <- setdiff(names(cells), "GSM_ID")
  setnames(cells, old, paste0("cell_", make.names(old)))
  fwrite(cells, file.path(opts$output_dir, sprintf("%s_epidish_cell_proportions.csv", cohort)))
  cells
}

candidate_by_method <- lapply(meta_candidates, function(x) unique(x$cpg))
base_dirs <- list(beta = opts$beta_dir, chalm = opts$chalm_dir, camda = opts$camda_dir)

load_cohort <- function(cohort) {
  pheno <- read_pheno(cohort)
  beta_path <- find_matrix_file(cohort, "beta", opts$beta_dir)
  message("[MATRIX] ", cohort, " / beta: ", beta_path)
  beta_full <- read_matrix(beta_path, pheno$GSM_ID)
  pheno[, Age_model := Age_observed]
  pheno[, Age_source := fifelse(is.finite(Age_observed), "observed", "missing")]
  fwrite(
    pheno[, .(GSM_ID, GSE_ID, Disease_State2, outcome, Age_observed, Age_model, Age_source)],
    file.path(opts$output_dir, sprintf("%s_age_used.csv", cohort))
  )
  message(
    "[AGE] ", cohort, ": observed=", sum(pheno$Age_source == "observed"),
    ", missing=", sum(pheno$Age_source == "missing")
  )
  cells <- estimate_cells(beta_full, cohort)
  beta_target <- beta_full[rownames(beta_full) %in% candidate_by_method$beta, , drop = FALSE]
  rm(beta_full)
  gc()

  matrices <- list(beta = to_m_value(beta_target))
  for (method in c("chalm", "camda")) {
    if (length(candidate_by_method[[method]]) == 0L) {
      matrices[[method]] <- matrix(
        numeric(), nrow = 0L, ncol = length(pheno$GSM_ID),
        dimnames = list(character(), pheno$GSM_ID)
      )
    } else {
      path <- find_matrix_file(cohort, method, base_dirs[[method]])
      message("[MATRIX] ", cohort, " / ", method, ": ", path)
      matrices[[method]] <- to_m_value(read_matrix(path, pheno$GSM_ID, candidate_by_method[[method]]))
    }
  }

  sample_sets <- c(list(pheno$GSM_ID, cells$GSM_ID), lapply(matrices, colnames))
  samples <- Reduce(intersect, sample_sets)
  if (length(samples) < 10L) stop("Too few samples shared across phenotype, cells, and matrices for ", cohort)
  pheno <- pheno[match(samples, GSM_ID)]
  observed_counts <- c(sum(pheno$outcome == 1L), sum(pheno$outcome == 0L))
  if (!identical(as.integer(observed_counts), EXPECTED_COUNTS[[cohort]])) {
    stop(cohort, ": matrix matching changed the Method case/control counts to ",
         paste(observed_counts, collapse = "/"))
  }
  cells <- cells[match(samples, GSM_ID)]
  matrices <- lapply(matrices, function(x) x[, samples, drop = FALSE])
  pheno <- merge(pheno, cells, by = "GSM_ID", all.x = TRUE, sort = FALSE)
  pheno <- pheno[match(samples, GSM_ID)]

  data.table(
    cohort = cohort,
    method = METHODS,
    requested_features = vapply(candidate_by_method, length, integer(1)),
    observed_features = vapply(matrices, nrow, integer(1)),
    samples = length(samples)
  ) -> audit
  list(pheno = pheno, matrices = matrices, audit = audit)
}

cohort_data <- setNames(lapply(SCORE_COHORTS, load_cohort), SCORE_COHORTS)
fwrite(rbindlist(lapply(cohort_data, `[[`, "audit")), file.path(opts$output_dir, "matrix_availability_audit.csv"))

common_features <- lapply(METHODS, function(method) {
  Reduce(intersect, lapply(TUNING_COHORTS, function(cohort) {
    mat <- cohort_data[[cohort]]$matrices[[method]]
    rownames(mat)[rowSums(!is.finite(mat)) == 0L]
  }))
})
names(common_features) <- METHODS
all_cohort_features <- lapply(METHODS, function(method) {
  Reduce(intersect, lapply(SCORE_COHORTS, function(cohort) {
    mat <- cohort_data[[cohort]]$matrices[[method]]
    rownames(mat)[rowSums(!is.finite(mat)) == 0L]
  }))
})
names(all_cohort_features) <- METHODS

usable_features <- rbindlist(lapply(METHODS, function(method) {
  meta_candidates[[method]][
    cpg %in% common_features[[method]] & cpg %in% candidate_by_method[[method]]
  ]
}), fill = TRUE)
usable_features[, feature_group := ifelse(value_method == "beta", "EWAS", "haEWAS")]
usable_features[, feature_id := paste(cpg, toupper(value_method), sep = "-")]
setorder(usable_features, feature_group, meta_p, feature_id)
all_cohort_usable_features <- rbindlist(lapply(METHODS, function(method) {
  meta_candidates[[method]][
    cpg %in% all_cohort_features[[method]] & cpg %in% candidate_by_method[[method]]
  ]
}), fill = TRUE)
all_cohort_usable_features[, feature_group := ifelse(value_method == "beta", "EWAS", "haEWAS")]
all_cohort_usable_features[, feature_id := paste(cpg, toupper(value_method), sep = "-")]
setorder(all_cohort_usable_features, feature_group, meta_p, feature_id)
fwrite(all_cohort_usable_features,
       file.path(opts$output_dir, "mrs_features_available_in_all_cohorts.csv"))
fwrite(usable_features,
       file.path(opts$output_dir, "mrs_features_available_in_discovery.csv"))

ranked_ewas_features <- usable_features[feature_group == "EWAS"]
ranked_ha_features <- usable_features[feature_group == "haEWAS"]
if (nrow(ranked_ewas_features) < opts$top_min || nrow(ranked_ha_features) < opts$top_min) {
  stop("Too few EWAS or haEWAS features remain after cross-cohort availability QC.")
}
message("[FEATURE QC] Ranked EWAS candidate features: ", nrow(ranked_ewas_features))
message("[FEATURE QC] Ranked haEWAS candidate features: ", nrow(ranked_ha_features),
        " (CHALM=", sum(ranked_ha_features$value_method == "chalm"),
        ", CAMDA=", sum(ranked_ha_features$value_method == "camda"), ")")

calculate_mrs <- function(cohort, features) {
  samples <- cohort_data[[cohort]]$pheno$GSM_ID
  score <- setNames(rep(0, length(samples)), samples)
  for (method in unique(features$value_method)) {
    sub <- features[value_method == method]
    mat <- cohort_data[[cohort]]$matrices[[method]][sub$cpg, samples, drop = FALSE]
    weights <- setNames(sub$meta_weight, sub$cpg)
    score <- score + as.numeric(t(mat) %*% weights[rownames(mat)])
  }
  score
}

signed_auc <- function(score, outcome) {
  keep <- is.finite(score) & !is.na(outcome)
  score <- score[keep]
  outcome <- outcome[keep]
  n_case <- sum(outcome == 1L)
  n_control <- sum(outcome == 0L)
  if (n_case == 0L || n_control == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[outcome == 1L]) - n_case * (n_case + 1) / 2) / (n_case * n_control)
}

stratified_bootstrap_index <- function(outcome) {
  case_index <- which(outcome == 1L)
  control_index <- which(outcome == 0L)
  c(
    sample(case_index, length(case_index), replace = TRUE),
    sample(control_index, length(control_index), replace = TRUE)
  )
}

ranked_groups <- list(EWAS = ranked_ewas_features, haEWAS = ranked_ha_features)
top_n_grid <- list(
  EWAS = TOP_N[TOP_N <= nrow(ranked_ewas_features)],
  haEWAS = TOP_N[TOP_N <= nrow(ranked_ha_features)]
)
if (!length(top_n_grid$EWAS)) stop("No Top-N values are available for EWAS.")
if (!length(top_n_grid$haEWAS)) stop("No Top-N values are available for haEWAS.")
development_scores <- list()
tuning_rows <- list()

for (group_i in names(ranked_groups)) {
  development_scores[[group_i]] <- list()
  for (top_n_i in top_n_grid[[group_i]]) {
    selected_i <- head(ranked_groups[[group_i]], top_n_i)
    development_scores[[group_i]][[as.character(top_n_i)]] <- setNames(
      lapply(TUNING_COHORTS, function(cohort_i) calculate_mrs(cohort_i, selected_i)),
      TUNING_COHORTS
    )
    fold_auc <- vapply(TUNING_COHORTS, function(cohort_i) {
      ph <- cohort_data[[cohort_i]]$pheno
      score <- development_scores[[group_i]][[as.character(top_n_i)]][[cohort_i]]
      signed_auc(unname(score[ph$GSM_ID]), ph$outcome)
    }, numeric(1))
    tuning_row <- data.table(
      feature_group = group_i,
      top_n = top_n_i,
      mean_signed_auc = mean(fold_auc),
      min_signed_auc = min(fold_auc),
      max_signed_auc = max(fold_auc)
    )
    for (cohort_i in names(fold_auc)) {
      set(tuning_row, j = paste0("auc_", make.names(cohort_i)), value = unname(fold_auc[[cohort_i]]))
    }
    tuning_rows[[paste(group_i, top_n_i, sep = "__")]] <- tuning_row
  }
}
tuning_dt <- rbindlist(tuning_rows)

set.seed(opts$bootstrap_seed)
bootstrap_rows <- vector("list", opts$bootstrap_reps)
for (bootstrap_id in seq_len(opts$bootstrap_reps)) {
  sampled_index <- setNames(lapply(TUNING_COHORTS, function(cohort_i) {
    stratified_bootstrap_index(cohort_data[[cohort_i]]$pheno$outcome)
  }), TUNING_COHORTS)
  rows_i <- list()
  for (group_i in names(ranked_groups)) {
    for (top_n_i in top_n_grid[[group_i]]) {
      fold_auc <- vapply(TUNING_COHORTS, function(cohort_i) {
        ph <- cohort_data[[cohort_i]]$pheno
        idx <- sampled_index[[cohort_i]]
        score <- development_scores[[group_i]][[as.character(top_n_i)]][[cohort_i]]
        signed_auc(unname(score[ph$GSM_ID][idx]), ph$outcome[idx])
      }, numeric(1))
      rows_i[[paste(group_i, top_n_i, sep = "__")]] <- data.table(
        bootstrap_id = bootstrap_id,
        feature_group = group_i,
        top_n = top_n_i,
        mean_signed_auc = mean(fold_auc)
      )
    }
  }
  bootstrap_rows[[bootstrap_id]] <- rbindlist(rows_i)
}
bootstrap_dt <- rbindlist(bootstrap_rows)
bootstrap_summary <- bootstrap_dt[, .(
  bootstrap_mean_auc = mean(mean_signed_auc),
  bootstrap_sd_auc = sd(mean_signed_auc),
  bootstrap_ci_low = quantile(mean_signed_auc, 0.025, names = FALSE),
  bootstrap_ci_high = quantile(mean_signed_auc, 0.975, names = FALSE)
), by = .(feature_group, top_n)]
tuning_dt <- merge(tuning_dt, bootstrap_summary, by = c("feature_group", "top_n"), all.x = TRUE)
setorder(tuning_dt, feature_group, top_n)

select_one_se <- function(group_i) {
  dt <- copy(tuning_dt[feature_group == group_i])
  best_observed <- dt[order(-mean_signed_auc, top_n)][1]
  best_bootstrap <- dt[order(-bootstrap_mean_auc, top_n)][1]
  threshold <- best_bootstrap$bootstrap_mean_auc - best_bootstrap$bootstrap_sd_auc
  selected <- dt[bootstrap_mean_auc >= threshold][order(top_n)][1]
  selected[, `:=`(
    best_observed_top_n = best_observed$top_n,
    best_observed_mean_auc = best_observed$mean_signed_auc,
    best_bootstrap_top_n = best_bootstrap$top_n,
    best_bootstrap_mean_auc = best_bootstrap$bootstrap_mean_auc,
    one_se_threshold = threshold,
    selection_rule = "smallest Top-N within one bootstrap-estimated SE of maximum equal-cohort mean bootstrap AUC"
  )]
  selected
}

selected_n <- rbindlist(lapply(names(ranked_groups), select_one_se))
setorder(selected_n, feature_group)
fwrite(tuning_dt, file.path(opts$output_dir, "mrs_topn_tuning_summary.csv"))
fwrite(bootstrap_dt, file.path(opts$output_dir, "mrs_topn_bootstrap_replicates.csv"))
fwrite(selected_n, file.path(opts$output_dir, "mrs_selected_topn.csv"))

best_ewas_n <- selected_n[feature_group == "EWAS", top_n][[1]]
best_haewas_n <- selected_n[feature_group == "haEWAS", top_n][[1]]
ewas_features <- copy(head(ranked_ewas_features, best_ewas_n))
ha_features <- copy(head(ranked_ha_features, best_haewas_n))
ewas_features[, `:=`(selected_top_n = best_ewas_n, selected_rank = seq_len(.N))]
ha_features[, `:=`(selected_top_n = best_haewas_n, selected_rank = seq_len(.N))]
for (cohort in TEST_COHORTS) {
  required_features <- list(
    beta = ewas_features$cpg,
    chalm = ha_features[value_method == "chalm", cpg],
    camda = ha_features[value_method == "camda", cpg]
  )
  for (method in names(required_features)) {
    missing <- setdiff(required_features[[method]], rownames(cohort_data[[cohort]]$matrices[[method]]))
    if (length(missing)) {
      stop(cohort, " lacks ", length(missing), " frozen ", method,
           " features; external features cannot be used to alter discovery selection.")
    }
    frozen_matrix <- cohort_data[[cohort]]$matrices[[method]][required_features[[method]], , drop = FALSE]
    if (any(!is.finite(frozen_matrix))) {
      stop(cohort, " has missing values among frozen ", method,
           " features; imputation is not performed.")
    }
  }
}
frozen_features <- rbindlist(list(ewas_features, ha_features), fill = TRUE)
fwrite(frozen_features, file.path(opts$output_dir, "mrs_frozen_features.csv"))

message("[TOP-N] Selected EWAS N=", best_ewas_n, "; haEWAS N=", best_haewas_n)
message("[TOP-N] Frozen features are used unchanged in all external cohorts.")

score_tables <- lapply(SCORE_COHORTS, function(cohort_i) {
  ph <- copy(cohort_data[[cohort_i]]$pheno)
  ewas_score <- calculate_mrs(cohort_i, ewas_features)
  haewas_score <- calculate_mrs(cohort_i, ha_features)
  ph[, cohort := cohort_i]
  ph[, EWAS_MRS := ewas_score[GSM_ID]]
  ph[, haEWAS_MRS := haewas_score[GSM_ID]]
  ph
})
names(score_tables) <- SCORE_COHORTS

all_scores <- rbindlist(score_tables, fill = TRUE)
fwrite(all_scores, file.path(opts$output_dir, "mrs_sample_scores.csv"))

train_cells <- score_tables[[TRAIN_COHORT]]
cell_cols <- grep("^cell_", names(train_cells), value = TRUE)
cell_cols <- cell_cols[vapply(cell_cols, function(x) {
  vals <- safe_numeric(train_cells[[x]])
  sum(is.finite(vals)) >= 10L && is.finite(sd(vals, na.rm = TRUE)) &&
    sd(vals, na.rm = TRUE) > 0
}, logical(1))]
if (length(cell_cols) < 2L) stop("Training cohort lacks variable EpiDISH cell proportions.")
reference_cell <- if ("cell_Neutro" %in% cell_cols) "cell_Neutro" else if (length(cell_cols) > 0L) tail(cell_cols, 1L) else NA_character_
model_cell_cols <- setdiff(cell_cols, reference_cell)
message("[COVARIATES] EpiDISH cells retained: ", paste(model_cell_cols, collapse = ", "))
message("[COVARIATES] Compositional reference dropped: ", reference_cell)

train <- copy(score_tables[[TRAIN_COHORT]])
train_age_sd <- sd(train$Age_model, na.rm = TRUE)
if (any(!is.finite(train$Age_model)) || !is.finite(train_age_sd) || train_age_sd <= 0) {
  stop("Training cohort does not have enough observed age values for adjustment.")
}
if (any(!is.finite(train$Sex_model)) ||
    length(unique(na.omit(train$Sex_model))) < 2L) {
  stop("Training cohort does not have enough observed sex values for adjustment.")
}
base_covariates <- c("Age_model", "Sex_model", model_cell_cols)

covariate_audit <- rbindlist(lapply(SCORE_COHORTS, function(cohort) {
  dt <- score_tables[[cohort]]
  missing_columns <- setdiff(base_covariates, names(dt))
  if (length(missing_columns)) {
    stop(cohort, " lacks model covariates: ", paste(missing_columns, collapse = ", "))
  }
  for (column in base_covariates) {
    if (any(!is.finite(safe_numeric(dt[[column]])))) {
      stop(cohort, " has missing ", column, " values; imputation is not performed.")
    }
  }
  data.table(
    cohort = cohort,
    variable = base_covariates,
    n_missing = 0L,
    handling = "not_imputed"
  )
}))
fwrite(covariate_audit, file.path(opts$output_dir, "mrs_covariate_imputation_audit.csv"))
train <- score_tables[[TRAIN_COHORT]]

rhs0 <- paste(base_covariates, collapse = " + ")
formulas <- setNames(
  list(
    as.formula(paste("outcome ~ EWAS_MRS +", rhs0)),
    as.formula(paste("outcome ~ haEWAS_MRS +", rhs0)),
    as.formula(paste("outcome ~ EWAS_MRS + haEWAS_MRS +", rhs0))
  ),
  c("EWAS", "haEWAS", "EWAS+haEWAS")
)
models <- lapply(formulas, function(f) glm(f, data = train, family = binomial()))

coefficient_tables <- rbindlist(lapply(names(models), function(model_name) {
  sm <- summary(models[[model_name]])$coefficients
  data.table(
    model = model_name,
    term = rownames(sm),
    estimate = sm[, "Estimate"],
    std_error = sm[, "Std. Error"],
    z = sm[, "z value"],
    p_value = sm[, "Pr(>|z|)"],
    converged = models[[model_name]]$converged
  )
}))
fwrite(coefficient_tables, file.path(opts$output_dir, "mrs_model_coefficients.csv"))

auc_ci <- function(roc_obj) {
  ci <- suppressMessages(as.numeric(pROC::ci.auc(roc_obj, method = "delong")))
  c(low = ci[[1]], auc = ci[[2]], high = ci[[3]])
}

performance <- list()
roc_rows <- list()
prediction_rows <- list()
delong_rows <- list()

for (cohort in SCORE_COHORTS) {
  dt <- copy(score_tables[[cohort]])
  cohort_rocs <- list()
  for (model_name in names(models)) {
    prob <- as.numeric(predict(models[[model_name]], newdata = dt, type = "response"))
    roc_obj <- pROC::roc(dt$outcome, prob, levels = c(0, 1), direction = "<", quiet = TRUE)
    cohort_rocs[[model_name]] <- roc_obj
    ci <- auc_ci(roc_obj)
    performance[[length(performance) + 1L]] <- data.table(
      training_cohort = TRAIN_COHORT,
      evaluation_cohort = cohort,
      evaluation_type = fifelse(
        cohort == TRAIN_COHORT,
        "development_training",
        fifelse(cohort %in% TUNING_COHORTS, "development_tuning", "external")
      ),
      model = model_name,
      N = nrow(dt),
      N_case = sum(dt$outcome == 1L),
      N_control = sum(dt$outcome == 0L),
      AUC = ci[["auc"]],
      AUC_CI_low = ci[["low"]],
      AUC_CI_high = ci[["high"]]
    )
    prediction_rows[[length(prediction_rows) + 1L]] <- data.table(
      cohort = cohort, GSM_ID = dt$GSM_ID, outcome = dt$outcome,
      model = model_name, predicted_probability = prob
    )
    roc_rows[[length(roc_rows) + 1L]] <- data.table(
      cohort = cohort, model = model_name,
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      threshold = roc_obj$thresholds
    )
  }
  test <- tryCatch(
    pROC::roc.test(cohort_rocs$EWAS, cohort_rocs$haEWAS, paired = TRUE, method = "delong"),
    error = function(e) NULL
  )
  delong_rows[[length(delong_rows) + 1L]] <- data.table(
    cohort = cohort,
    EWAS_AUC = as.numeric(pROC::auc(cohort_rocs$EWAS)),
    haEWAS_AUC = as.numeric(pROC::auc(cohort_rocs$haEWAS)),
    delta_haEWAS_minus_EWAS = as.numeric(pROC::auc(cohort_rocs$haEWAS) - pROC::auc(cohort_rocs$EWAS)),
    paired_DeLong_p = if (!is.null(test)) test$p.value else NA_real_
  )
}

performance_dt <- rbindlist(performance)
predictions_dt <- rbindlist(prediction_rows)
roc_dt <- rbindlist(roc_rows)
delong_dt <- rbindlist(delong_rows)

cohen_d <- function(case_score, control_score) {
  case_score <- case_score[is.finite(case_score)]
  control_score <- control_score[is.finite(control_score)]
  n_case <- length(case_score)
  n_control <- length(control_score)
  if (n_case < 2L || n_control < 2L) return(NA_real_)
  pooled_sd <- sqrt(
    ((n_case - 1L) * var(case_score) + (n_control - 1L) * var(control_score)) /
      (n_case + n_control - 2L)
  )
  if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NA_real_)
  (mean(case_score) - mean(control_score)) / pooled_sd
}

external_mrs_tests <- list()
observed_deltas <- setNames(numeric(length(TEST_COHORTS)), TEST_COHORTS)
for (cohort in TEST_COHORTS) {
  dt <- score_tables[[cohort]]
  case <- dt$outcome == 1L
  control <- dt$outcome == 0L
  ewas_score <- dt$EWAS_MRS
  haewas_score <- dt$haEWAS_MRS
  for (model_name in c("EWAS", "haEWAS")) {
    score <- if (model_name == "EWAS") ewas_score else haewas_score
    case_score <- score[case & is.finite(score)]
    control_score <- score[control & is.finite(score)]
    test <- wilcox.test(case_score, control_score, exact = FALSE)
    external_mrs_tests[[length(external_mrs_tests) + 1L]] <- data.table(
      cohort = cohort,
      MRS = model_name,
      N_case = length(case_score),
      N_control = length(control_score),
      mean_case = mean(case_score),
      mean_control = mean(control_score),
      median_case = median(case_score),
      median_control = median(control_score),
      wilcoxon_p = test$p.value,
      cohens_d = cohen_d(case_score, control_score)
    )
  }
  observed_deltas[[cohort]] <-
    external_mrs_tests[[length(external_mrs_tests)]]$cohens_d -
    external_mrs_tests[[length(external_mrs_tests) - 1L]]$cohens_d
}
external_mrs_tests_dt <- rbindlist(external_mrs_tests)

set.seed(opts$bootstrap_seed + 1L)
delta_boot <- matrix(
  NA_real_,
  nrow = opts$bootstrap_reps,
  ncol = length(TEST_COHORTS),
  dimnames = list(NULL, TEST_COHORTS)
)
for (bootstrap_id in seq_len(opts$bootstrap_reps)) {
  for (cohort in TEST_COHORTS) {
    dt <- score_tables[[cohort]]
    case_index <- which(dt$outcome == 1L)
    control_index <- which(dt$outcome == 0L)
    idx <- c(
      sample(case_index, length(case_index), replace = TRUE),
      sample(control_index, length(control_index), replace = TRUE)
    )
    outcome <- dt$outcome[idx]
    ewas_score <- dt$EWAS_MRS[idx]
    haewas_score <- dt$haEWAS_MRS[idx]
    delta_boot[bootstrap_id, cohort] <-
      cohen_d(haewas_score[outcome == 1L], haewas_score[outcome == 0L]) -
      cohen_d(ewas_score[outcome == 1L], ewas_score[outcome == 0L])
  }
}
if (any(!is.finite(observed_deltas)) || any(!is.finite(delta_boot))) {
  stop("Could not calculate all external Cohen's d differences for the bootstrap.")
}
delta_boot_mean <- rowMeans(delta_boot, na.rm = FALSE)
delta_ci_rows <- lapply(TEST_COHORTS, function(cohort) {
  x <- delta_boot[, cohort]
  data.table(
    cohort = cohort,
    comparison = "haEWAS_minus_EWAS",
    observed_delta_cohens_d = observed_deltas[[cohort]],
    ci_low = quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    ci_high = quantile(x, 0.975, na.rm = TRUE, names = FALSE),
    bootstrap_replicates = opts$bootstrap_reps
  )
})
delta_ci_rows[[length(delta_ci_rows) + 1L]] <- data.table(
  cohort = "equally_weighted_mean",
  comparison = "haEWAS_minus_EWAS",
  observed_delta_cohens_d = mean(observed_deltas),
  ci_low = quantile(delta_boot_mean, 0.025, na.rm = TRUE, names = FALSE),
  ci_high = quantile(delta_boot_mean, 0.975, na.rm = TRUE, names = FALSE),
  bootstrap_replicates = opts$bootstrap_reps
)
external_cohen_d_delta <- rbindlist(delta_ci_rows)
external_cohen_d_bootstrap <- as.data.table(delta_boot)
external_cohen_d_bootstrap[, bootstrap_id := seq_len(.N)]
setcolorder(external_cohen_d_bootstrap, c("bootstrap_id", TEST_COHORTS))
external_cohen_d_bootstrap[, equally_weighted_mean := delta_boot_mean]

fwrite(performance_dt, file.path(opts$output_dir, "mrs_model_performance.csv"))
fwrite(
  performance_dt[evaluation_type == "external"],
  file.path(opts$output_dir, "mrs_external_model_performance.csv")
)
fwrite(predictions_dt, file.path(opts$output_dir, "mrs_model_predictions.csv"))
fwrite(roc_dt, file.path(opts$output_dir, "mrs_roc_curve_points.csv"))
fwrite(delong_dt, file.path(opts$output_dir, "mrs_auc_comparison.csv"))
fwrite(
  delong_dt[cohort %in% TEST_COHORTS],
  file.path(opts$output_dir, "mrs_external_auc_comparison.csv")
)
fwrite(external_mrs_tests_dt,
       file.path(opts$output_dir, "mrs_external_mrs_distribution_tests.csv"))
fwrite(external_cohen_d_delta,
       file.path(opts$output_dir, "mrs_external_cohen_d_difference.csv"))
fwrite(external_cohen_d_bootstrap,
       file.path(opts$output_dir, "mrs_external_cohen_d_difference_bootstrap.csv"))

summary_dt <- data.table(
  item = c(
    "beta_meta_significant", "chalm_meta_significant", "camda_meta_significant",
    "union_significant_cpgs", "beta_meta_p005_candidates", "chalm_meta_p005_candidates",
    "camda_meta_p005_candidates", "optimized_EWAS_features", "optimized_haEWAS_features",
    "optimized_haEWAS_unique_cpgs"
  ),
  value = c(
    nrow(meta_sig$beta), nrow(meta_sig$chalm), nrow(meta_sig$camda),
    length(union_cpgs), nrow(meta_candidates$beta), nrow(meta_candidates$chalm),
    nrow(meta_candidates$camda), nrow(ewas_features), nrow(ha_features), uniqueN(ha_features$cpg)
  )
)
fwrite(summary_dt, file.path(opts$output_dir, "mrs_three_model_summary.csv"))

done_file <- file.path(opts$output_dir, "mrs_three_models.complete")
writeLines(c(
  paste0("completed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("train_cohort=", TRAIN_COHORT),
  paste0("tuning_cohorts=", paste(TUNING_COHORTS, collapse = ",")),
  paste0("test_cohorts=", paste(TEST_COHORTS, collapse = ",")),
  paste0("fdr_cutoff=", opts$fdr_cutoff),
  paste0("candidate_p_cutoff=", opts$candidate_p_cutoff),
  paste0("top_n_grid=", paste(sort(unique(c(top_n_grid$EWAS, top_n_grid$haEWAS))), collapse = ",")),
  paste0("ewas_top_n_grid=", paste(top_n_grid$EWAS, collapse = ",")),
  paste0("haewas_top_n_grid=", paste(top_n_grid$haEWAS, collapse = ",")),
  paste0("selected_ewas_top_n=", best_ewas_n),
  paste0("selected_haewas_top_n=", best_haewas_n),
  "selection_rule=smallest Top-N within one bootstrap-estimated SE of maximum equal-cohort mean bootstrap AUC",
  "external_rule=frozen features, weights, and model; no external-label orientation",
  paste0("ewas_features=", nrow(ewas_features)),
  paste0("haewas_features=", nrow(ha_features))
), done_file)

message("[DONE] Three-model MRS workflow completed.")
message("[DONE] Output: ", opts$output_dir)
print(performance_dt)
print(delong_dt)
