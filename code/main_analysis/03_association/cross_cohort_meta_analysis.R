# Meta-analyze cohort-level EWAS and haEWAS results.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(meta)
})

option_list <- list(
  make_option(c("--task_table"), type = "character"),
  make_option(c("--task_index"), type = "integer"),
  make_option(c("--method"), type = "character", help = "m, chalm, camda or haewas"),
  make_option(c("--smoking_file"), type = "character", default = NULL),
  make_option(c("--smoking_related"), action = "store_true", default = FALSE),
  make_option(c("--cores"), type = "integer", default = as.integer(Sys.getenv("META_CPUS", "1"))),
  make_option(c("--resume"), action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$task_table) || is.null(opts$task_index) || is.null(opts$method)) {
  stop("Required: --task_table --task_index --method")
}

method <- tolower(opts$method)
if (!method %in% c("m", "chalm", "camda", "haewas")) {
  stop("--method must be m, chalm, camda or haewas")
}
methods <- if (method == "m") "m" else c("chalm", "camda")
meta_cores <- if (is.finite(opts$cores) && opts$cores > 0L) opts$cores else 1L

tasks <- fread(opts$task_table, data.table = TRUE, check.names = FALSE)
if (opts$task_index < 1L || opts$task_index > nrow(tasks)) stop("Invalid task index")
task <- tasks[opts$task_index]
required <- c("phenotype_label", "tissue", "gse_list", "single_output_dir",
              "single_output_label", "output_dir", "output_label")
if (length(setdiff(required, names(task)))) stop("Meta task table is missing required columns")

phenotype_label <- as.character(task$phenotype_label[[1L]])
tissue <- as.character(task$tissue[[1L]])
gse_list <- unique(trimws(strsplit(as.character(task$gse_list[[1L]]), ",", fixed = TRUE)[[1L]]))
gse_list <- gse_list[nzchar(gse_list)]
if (length(gse_list) < 2L) stop("At least two cohorts are required")

single_dir <- file.path(as.character(task$single_output_dir[[1L]]),
                        as.character(task$single_output_label[[1L]]))
output_dir <- file.path(as.character(task$output_dir[[1L]]),
                        as.character(task$output_label[[1L]]))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

smoking_related <- isTRUE(opts$smoking_related) ||
  grepl("smok|tobacco|nicotine", phenotype_label, ignore.case = TRUE)
if ("smoking_related" %in% names(task) && !is.na(task$smoking_related[[1L]])) {
  smoking_related <- tolower(as.character(task$smoking_related[[1L]])[1L]) %in%
    c("1", "true", "yes", "y")
}

smoking_file <- opts$smoking_file
if (is.null(smoking_file) && "smoking_file" %in% names(task)) {
  smoking_file <- as.character(task$smoking_file[[1L]])
}
smoking_probes <- character()
if (!smoking_related) {
  if (is.null(smoking_file) || !nzchar(smoking_file) || !file.exists(smoking_file)) {
    stop("Provide --smoking_file for non-smoking phenotypes")
  }
  smoking <- fread(smoking_file, data.table = TRUE, check.names = FALSE)
  probe_col <- intersect(c("Probe ID", "Probe.ID", "cpg", "CpG", "Name"), names(smoking))[1L]
  p_col <- intersect(c("P-value", "P.value", "P_value", "p_value", "P"), names(smoking))[1L]
  if (is.na(probe_col) || is.na(p_col)) {
    stop("Smoking file must contain probe IDs and P-values")
  }
  smoking_p <- suppressWarnings(as.numeric(smoking[[p_col]]))
  smoking_probes <- unique(trimws(as.character(
    smoking[[probe_col]][is.finite(smoking_p) & smoking_p < 1e-7]
  )))
  smoking_probes <- smoking_probes[!is.na(smoking_probes) & nzchar(smoking_probes)]
}

output_files <- file.path(output_dir, paste0("meta_effectsize_", methods, "_", tissue, ".csv"))
if (isTRUE(opts$resume) && all(file.exists(output_files)) &&
    all(is.finite(file.info(output_files)$size) & file.info(output_files)$size > 0)) {
  cat("[RESUME] Existing meta results found; skipping.\n")
  quit(save = "no", status = 0)
}

read_single <- function(gse, data_type) {
  path <- file.path(single_dir, paste0(gse, "_", data_type, "_single_cpg_pVal_df.csv"))
  if (!file.exists(path)) {
    cat("[SKIP] Missing single-cohort result: ", basename(path), "\n", sep = "")
    return(NULL)
  }

  dt <- fread(path, data.table = TRUE, check.names = FALSE)
  count_cols <- c("n_samples", "n_cases", "n_controls")
  if (!all(count_cols %in% names(dt))) {
    stop("Missing sample-count columns in ", path,
         "; rerun the single-cohort analysis")
  }
  counts <- unique(dt[, ..count_cols])
  if (nrow(counts) != 1L) stop("Inconsistent sample counts in ", path)

  n_samples <- as.numeric(counts$n_samples[[1L]])
  n_cases <- as.numeric(counts$n_cases[[1L]])
  n_controls <- as.numeric(counts$n_controls[[1L]])
  binary <- is.finite(n_cases) && is.finite(n_controls)
  eligible <- if (binary) n_cases >= 50L && n_controls >= 50L else n_samples >= 50L
  if (!eligible) {
    count_text <- if (binary) {
      paste0("cases=", n_cases, ", controls=", n_controls)
    } else {
      paste0("complete samples=", n_samples)
    }
    cat("[SKIP] ", gse, "/", data_type, " does not meet the sample-size rule (",
        count_text, ").\n", sep = "")
    return(NULL)
  }

  cpg_col <- intersect(c("cpg", "cpg_id", "CpG", "Name"), names(dt))[1L]
  effect_col <- intersect(c("StdEffectSize", "Estimate_sd", "Estimate", "Eff"), names(dt))[1L]
  se_col <- intersect(c("StdSE", "SE_sd", "SE", "SEf"), names(dt))[1L]
  p_col <- intersect(c("pValue", "p_value", "pval", "P.Value", "Pvalue"), names(dt))[1L]
  if (anyNA(c(cpg_col, effect_col, se_col))) stop("Missing meta-analysis columns in ", path)

  effect <- suppressWarnings(as.numeric(dt[[effect_col]]))
  se <- suppressWarnings(as.numeric(dt[[se_col]]))
  p_value <- if (!is.na(p_col)) suppressWarnings(as.numeric(dt[[p_col]])) else
    rep(NA_real_, nrow(dt))
  derive <- !is.finite(p_value) | p_value < 0 | p_value > 1
  valid_effect <- is.finite(effect) & is.finite(se) & se > 0
  derive_ok <- derive & valid_effect
  p_value[derive_ok] <- 2 * stats::pnorm(-abs(effect[derive_ok] / se[derive_ok]))

  out <- data.table(cpg = trimws(as.character(dt[[cpg_col]])),
                    effect = effect, se = se, pValue = p_value)
  out <- out[!is.na(cpg) & nzchar(cpg) & is.finite(effect) & is.finite(se) & se > 0 &
               is.finite(pValue) & pValue >= 0 & pValue <= 1]
  if (anyDuplicated(out$cpg)) stop("Duplicate CpGs in ", path)
  setnames(out, c("effect", "se", "pValue"),
           paste0(gse, c("_effect", "_se", "_pValue")))
  out
}

meta_one <- function(row) {
  effect_cols <- grep("_effect$", names(row), value = TRUE)
  se_cols <- grep("_se$", names(row), value = TRUE)
  effects <- as.numeric(unlist(row[, ..effect_cols], use.names = FALSE))
  ses <- as.numeric(unlist(row[, ..se_cols], use.names = FALSE))
  ok <- is.finite(effects) & is.finite(ses) & ses > 0
  if (sum(ok) < 2L) return(NULL)

  fit <- tryCatch(
    meta::metagen(TE = effects[ok], seTE = ses[ok],
                  studlab = sub("_effect$", "", effect_cols[ok])),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  use_fixed <- is.na(fit$pval.Q) || fit$pval.Q >= 0.05
  data.table(
    cpg = row$cpg,
    n_cohorts = sum(ok),
    Estimate_fixed = fit$TE.common,
    SE_fixed = fit$seTE.common,
    pval_fixed = fit$pval.common,
    Estimate_random = fit$TE.random,
    SE_random = fit$seTE.random,
    pval_random = fit$pval.random,
    pval_Q = fit$pval.Q,
    I2 = fit$I2,
    meta_Estimate = if (use_fixed) fit$TE.common else fit$TE.random,
    meta_SE = if (use_fixed) fit$seTE.common else fit$seTE.random,
    pval_final = if (use_fixed) fit$pval.common else fit$pval.random
  )
}

run_method <- function(data_type) {
  pieces <- setNames(lapply(gse_list, read_single, data_type = data_type), gse_list)
  used_gse <- names(pieces)[!vapply(pieces, is.null, logical(1))]
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) < 2L) {
    stop("Fewer than two eligible cohorts for ", phenotype_label, "/", data_type)
  }

  merged <- Reduce(function(x, y) merge(x, y, by = "cpg", all = TRUE), pieces)
  blocks <- split(seq_len(nrow(merged)), ceiling(seq_len(nrow(merged)) / 5000L))
  meta_parts <- if (.Platform$OS.type != "windows" && meta_cores > 1L && length(blocks) > 1L) {
    parallel::mclapply(blocks, function(idx) lapply(idx, function(i) meta_one(merged[i])),
                       mc.cores = meta_cores)
  } else {
    list(lapply(seq_len(nrow(merged)), function(i) meta_one(merged[i])))
  }
  rows <- unlist(meta_parts, recursive = FALSE)
  out <- rbindlist(Filter(Negate(is.null), rows), fill = TRUE)
  if (!nrow(out)) stop("No valid meta-analysis rows for ", phenotype_label, "/", data_type)
  out[, c("method", "phenotype_label", "tissue", "cohorts") :=
        list(data_type, phenotype_label, tissue, paste(used_gse, collapse = ","))]
  out
}

results <- lapply(methods, run_method)
names(results) <- methods
out <- rbindlist(results, use.names = TRUE, fill = TRUE)
out[, fdr_method_specific := p.adjust(pval_final, method = "BH")]
out[, fdr := fdr_method_specific]
out[, fdr_family := if (length(methods) == 2L) "haEWAS_joint" else "EWAS"]

if (length(smoking_probes)) {
  n_removed <- sum(out$cpg %chin% smoking_probes)
  out <- out[!cpg %chin% smoking_probes]
  cat("[INFO] Removed ", n_removed, " smoking-associated CpGs from final results.\n", sep = "")
}

for (data_type in methods) {
  result <- out[method == data_type]
  setorder(result, pval_final)
  outfile <- file.path(output_dir, paste0("meta_effectsize_", data_type, "_", tissue, ".csv"))
  fwrite(result, outfile)
  cat("[DONE] ", data_type, " meta: ", basename(outfile), " (",
      nrow(result), " CpGs; ", unique(result$fdr_family), " FDR)\n", sep = "")
}

if ("m" %in% methods) {
  significant <- out[method == "m" & fdr < 0.05]
  setorder(significant, fdr, pval_final)
  outfile <- file.path(output_dir, paste0("meta_significant_EWAS_", tissue, ".csv"))
  fwrite(significant, outfile)
  cat("[DONE] EWAS significant CpGs: ", nrow(significant), "\n", sep = "")
}

if (all(c("chalm", "camda") %in% methods)) {
  ha_significant_cpgs <- unique(
    out[method %chin% c("chalm", "camda") & fdr < 0.05, cpg]
  )
  ha_significant_rows <- out[
    method %chin% c("chalm", "camda") & cpg %chin% ha_significant_cpgs
  ]
  ha_significant <- data.table(
    cpg = character(), SignificantMethods = character(),
    CHALM_MetaEstimate = numeric(), CHALM_P = numeric(), CHALM_FDR = numeric(),
    CAMDA_MetaEstimate = numeric(), CAMDA_P = numeric(), CAMDA_FDR = numeric()
  )
  if (nrow(ha_significant_rows)) {
    ha_significant <- ha_significant_rows[, .(
      SignificantMethods = paste(method[fdr < 0.05], collapse = ";"),
      CHALM_MetaEstimate = meta_Estimate[method == "chalm"][1L],
      CHALM_P = pval_final[method == "chalm"][1L],
      CHALM_FDR = fdr[method == "chalm"][1L],
      CAMDA_MetaEstimate = meta_Estimate[method == "camda"][1L],
      CAMDA_P = pval_final[method == "camda"][1L],
      CAMDA_FDR = fdr[method == "camda"][1L]
    ), by = "cpg"]
  }
  setorder(ha_significant, cpg)
  outfile <- file.path(output_dir, paste0("meta_significant_haEWAS_", tissue, ".csv"))
  fwrite(ha_significant, outfile)
  cat("[DONE] haEWAS significant CpGs: ", nrow(ha_significant),
      " unique loci after pooled BH\n", sep = "")
}
