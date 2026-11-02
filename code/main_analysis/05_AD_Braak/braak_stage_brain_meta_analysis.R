# Meta-analyze brain Braak associations and export significant loci.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(meta)
})

option_list <- list(
  make_option("--manifest", type = "character"),
  make_option("--output-dir", type = "character"),
  make_option("--smoking-file", type = "character"),
  make_option("--feature-column", type = "character", default = "feature_id"),
  make_option("--effect-column", type = "character", default = "effect"),
  make_option("--se-column", type = "character", default = "se")
)
opts <- parse_args(OptionParser(option_list = option_list))
if (is.null(opts$manifest) || is.null(opts$output_dir) || is.null(opts$smoking_file)) {
  stop("--manifest, --output-dir, and --smoking-file are required.")
}

manifest <- fread(opts$manifest, data.table = TRUE, check.names = FALSE,
                  showProgress = FALSE)
if (!all(c("cohort", "method", "path") %in% names(manifest))) {
  stop("The manifest must contain cohort, method, and path columns.")
}
manifest[, `:=`(
  cohort = trimws(as.character(cohort)),
  method = tolower(trimws(as.character(method))),
  path = trimws(as.character(path))
)]
methods <- c("beta", "chalm", "camda")
if (anyNA(manifest$cohort) || anyNA(manifest$method) || anyNA(manifest$path) ||
    any(!nzchar(manifest$cohort)) || any(!nzchar(manifest$path)) ||
    any(!manifest$method %chin% methods) || anyDuplicated(manifest[, .(cohort, method)])) {
  stop("The manifest contains invalid or duplicate cohort-method rows.")
}
cohorts <- sort(unique(manifest$cohort))
if (length(cohorts) != 3L || any(grepl("london", cohorts, ignore.case = TRUE))) {
  stop("The manifest must contain the three brain discovery cohorts and exclude London PFC.")
}
for (method in methods) {
  rows <- manifest[method == ..method]
  if (nrow(rows) != 3L || !setequal(rows$cohort, cohorts)) {
    stop("Each methylation method must include the same three brain cohorts: ", method)
  }
}
if (any(!file.exists(manifest$path))) stop("An input association file is missing.")

smoking <- fread(opts$smoking_file, data.table = TRUE, check.names = FALSE,
                 showProgress = FALSE)
probe_col <- intersect(c("Probe ID", "Probe.ID", "cpg", "CpG", "Name"), names(smoking))[1L]
p_col <- intersect(c("P-value", "P.value", "P_value", "p_value", "P"), names(smoking))[1L]
if (is.na(probe_col) || is.na(p_col)) stop("The smoking file must contain probe IDs and P-values.")
smoking_p <- suppressWarnings(as.numeric(smoking[[p_col]]))
smoking_probes <- unique(trimws(as.character(
  smoking[[probe_col]][is.finite(smoking_p) & smoking_p < 1e-7]
)))
smoking_probes <- smoking_probes[!is.na(smoking_probes) & nzchar(smoking_probes)]

read_effects <- function(row) {
  path <- row$path[[1L]]
  header <- names(fread(path, nrows = 0L, data.table = FALSE, check.names = FALSE))
  required <- c(opts$feature_column, opts$effect_column, opts$se_column)
  if (!all(required %in% header)) stop("An input association table lacks effect columns: ", path)
  dt <- fread(path, select = required, data.table = TRUE, check.names = FALSE,
              showProgress = FALSE)
  setnames(dt, required, c("feature_id", "effect", "se"))
  dt[, `:=`(
    feature_id = trimws(as.character(feature_id)),
    effect = suppressWarnings(as.numeric(effect)),
    se = suppressWarnings(as.numeric(se)),
    cohort = row$cohort[[1L]],
    method = row$method[[1L]]
  )]
  dt <- dt[!is.na(feature_id) & nzchar(feature_id) & is.finite(effect) &
             is.finite(se) & se > 0]
  if (anyDuplicated(dt$feature_id)) stop("Duplicate feature IDs in ", path)
  dt
}

all_effects <- rbindlist(lapply(seq_len(nrow(manifest)), function(i) {
  read_effects(manifest[i])
}), use.names = TRUE)

combine_feature <- function(feature, data_type) {
  dt <- all_effects[method == data_type & feature_id == feature]
  if (nrow(dt) < 2L) return(NULL)
  fit <- tryCatch(
    meta::metagen(TE = dt$effect, seTE = dt$se, studlab = dt$cohort),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  p_heterogeneity <- suppressWarnings(as.numeric(fit$pval.Q))
  if (!is.finite(p_heterogeneity)) return(NULL)
  use_fixed <- p_heterogeneity > 0.05
  list(
    feature_id = feature,
    effect = as.numeric(if (use_fixed) fit$TE.fixed else fit$TE.random),
    se = as.numeric(if (use_fixed) fit$seTE.fixed else fit$seTE.random),
    p = as.numeric(if (use_fixed) fit$pval.fixed else fit$pval.random),
    fdr = NA_real_,
    fdr_method_specific = NA_real_,
    n_cohorts = nrow(dt),
    cohorts = paste(dt$cohort, collapse = ";"),
    effect_fixed = as.numeric(fit$TE.fixed),
    se_fixed = as.numeric(fit$seTE.fixed),
    p_fixed = as.numeric(fit$pval.fixed),
    effect_random = as.numeric(fit$TE.random),
    se_random = as.numeric(fit$seTE.random),
    p_random = as.numeric(fit$pval.random),
    p_heterogeneity = p_heterogeneity,
    I2 = as.numeric(fit$I2),
    model = if (use_fixed) "fixed" else "random",
    method = data_type
  )
}

meta_results <- setNames(lapply(methods, function(data_type) {
  features <- unique(all_effects[method == data_type, feature_id])
  out <- rbindlist(lapply(features, combine_feature, data_type = data_type), fill = TRUE)
  if (!nrow(out)) stop("No valid meta-analysis rows for ", data_type)
  out
}), methods)

meta_results$beta[, `:=`(
  fdr_method_specific = p.adjust(p, method = "BH"),
  fdr_family = "EWAS"
)]
meta_results$beta[, fdr := fdr_method_specific]

haewas <- rbindlist(meta_results[c("chalm", "camda")], use.names = TRUE)
haewas[, `:=`(
  fdr_method_specific = p.adjust(p, method = "BH"),
  fdr_family = "haEWAS_joint"
)]
haewas[, fdr := fdr_method_specific]
meta_results$chalm <- haewas[method == "chalm"]
meta_results$camda <- haewas[method == "camda"]

dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
smoking_audit <- vector("list", length(methods))
for (i in seq_along(methods)) {
  data_type <- methods[[i]]
  result <- meta_results[[data_type]]
  n_removed <- sum(result$feature_id %chin% smoking_probes)
  meta_results[[data_type]] <- result[!feature_id %chin% smoking_probes]
  smoking_audit[[i]] <- data.table(
    method = data_type,
    tested_before_smoking_filter = nrow(result),
    smoking_probes_removed = n_removed,
    tested_after_smoking_filter = nrow(meta_results[[data_type]])
  )
}
smoking_audit <- rbindlist(smoking_audit)
fwrite(smoking_audit, file.path(opts$output_dir, "braak_smoking_filter_summary.csv"))

for (method in methods) {
  result <- meta_results[[method]]
  setorder(result, p)
  fwrite(result, file.path(opts$output_dir, paste0("braak_meta_", method, ".csv")))
}

ewas_significant <- meta_results$beta[fdr < 0.05]
setorder(ewas_significant, fdr, p)
fwrite(ewas_significant, file.path(opts$output_dir, "braak_significant_EWAS.csv"))

haewas_all <- rbindlist(meta_results[c("chalm", "camda")], use.names = TRUE)
significant_cpgs <- unique(haewas_all[fdr < 0.05, feature_id])
haewas_significant <- data.table(
  feature_id = character(), SignificantMethods = character(),
  CHALM_effect = numeric(), CHALM_P = numeric(), CHALM_FDR = numeric(),
  CAMDA_effect = numeric(), CAMDA_P = numeric(), CAMDA_FDR = numeric()
)
if (length(significant_cpgs)) {
  haewas_significant <- haewas_all[feature_id %chin% significant_cpgs, .(
    SignificantMethods = paste(toupper(method[fdr < 0.05]), collapse = ";"),
    CHALM_effect = effect[method == "chalm"][1L],
    CHALM_P = p[method == "chalm"][1L],
    CHALM_FDR = fdr[method == "chalm"][1L],
    CAMDA_effect = effect[method == "camda"][1L],
    CAMDA_P = p[method == "camda"][1L],
    CAMDA_FDR = fdr[method == "camda"][1L]
  ), by = "feature_id"]
}
setorder(haewas_significant, feature_id)
fwrite(haewas_significant, file.path(opts$output_dir, "braak_significant_haEWAS.csv"))

cat("[DONE] EWAS significant CpGs: ", nrow(ewas_significant), "\n", sep = "")
cat("[DONE] haEWAS significant unique CpGs after pooled BH: ",
    nrow(haewas_significant), "\n", sep = "")
cat("[DONE] Removed smoking-associated probes with P < 1e-7.\n")
