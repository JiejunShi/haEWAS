# Run EWAS and haEWAS for selected cohorts.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
  library(EpiDISH)
  library(EpiSCORE)
  library(cets)
  library(dnaMethyAge)
  library(minfi)
  library(impute)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

methylation_root <- Sys.getenv("METHYLATION_ROOT")
if (!nzchar(methylation_root)) stop("Set METHYLATION_ROOT")

option_list <- list(
  make_option(c("--column"),                  type="character", help="Column name for filtering"),
  make_option(c("--selected_value"),          type="character", default=NULL, help="Value to select for binary phenotype"),
  make_option(c("--selected_control_value"),  type="character", default=NULL, help="Control group value for binary phenotype (optional)"),
  make_option(c("--subset_column"),           type="character", default=NULL, help="Optional phenotype column used to subset samples before continuous EWAS"),
  make_option(c("--subset_value"),            type="character", default=NULL, help="Optional phenotype value used with --subset_column"),
  make_option(c("--pheno_type"),              type="character", default="binary", help="Phenotype type: binary or continuous"),
  make_option(c("--tissue"),                  type="character", help="Tissue type (e.g., brain, whole_blood)"),
  make_option(c("--output_dir"),              type="character", help="Base output directory"),
  make_option(c("--output_label"),            type="character", default=NULL, help="Optional final output subdirectory name"),
  make_option(c("--gse_list"),                type="character", help="Comma-separated GSE IDs"),
  make_option(c("--phenotype_path"),          type="character",
              default=NULL,
              help="Path to phenotype CSV"),
  make_option(c("--beta_dir"),                type="character",
              default=file.path(methylation_root, "beta"),
              help="Directory for beta matrices"),
  make_option(c("--mvalue_dir"),              type="character", default=NULL,
              help="Directory for ComBat-adjusted M-value matrices"),
  make_option(c("--chalm_camda_dir"),         type="character",
              default=file.path(methylation_root, "0.CAMDA_new"),
              help="Directory for CHALM/CAMDA matrices"),
  make_option(c("--chalm_dir"),               type="character", default=NULL,
              help="Separate CHALM matrix directory; overrides --chalm_camda_dir for CHALM"),
  make_option(c("--camda_dir"),               type="character", default=NULL,
              help="Separate CAMDA matrix directory; overrides --chalm_camda_dir for CAMDA"),
  make_option(c("--cores"),                   type="integer",   default=8,   help="Number of cores for parallel processing"),
  make_option(c("--num"),                     type="integer",   default=-1,  help="Force PC number (-1=auto, 0=Force 0 PCs, or 5, 10...)")
)
opts <- parse_args(OptionParser(option_list=option_list))

is_blank <- function(x) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  x <- as.character(x[1])
  is.na(x) || trimws(x) == ""
}

required_args <- c("column","pheno_type","tissue","output_dir","gse_list","phenotype_path")
missing_args <- required_args[sapply(required_args, function(x) is.null(opts[[x]]))]
if (length(missing_args) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing_args, collapse=", ")))
}

column_name            <- opts$column
selected_value         <- opts$selected_value
selected_control_value <- opts$selected_control_value
subset_column          <- opts$subset_column
subset_value           <- opts$subset_value
pheno_type             <- tolower(opts$pheno_type)
tissue                 <- opts$tissue
output_base            <- opts$output_dir
output_label           <- opts$output_label
gse_list               <- trimws(strsplit(opts$gse_list, ",")[[1]])
gse_list               <- gse_list[nzchar(gse_list)]
phenotype_path         <- opts$phenotype_path
beta_dir               <- opts$beta_dir
mvalue_dir             <- if (is_blank(opts$mvalue_dir)) beta_dir else opts$mvalue_dir
chalm_camda_dir        <- opts$chalm_camda_dir
chalm_dir              <- if (is_blank(opts$chalm_dir)) chalm_camda_dir else opts$chalm_dir
camda_dir              <- if (is_blank(opts$camda_dir)) chalm_camda_dir else opts$camda_dir
cores                  <- opts$cores
num_fixed              <- opts$num

if (is_blank(output_label)) {
  output_label <- if (pheno_type == "continuous") {
    column_name
  } else if (!is_blank(selected_value)) {
    selected_value
  } else {
    column_name
  }
}

if (!is_blank(output_label)) {
  output_dir <- file.path(output_base, output_label)
} else {
  output_dir <- output_base
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

AGE_MISS_THR <- 0.10
MAX_RESID_PC <- suppressWarnings(as.integer(Sys.getenv("MAX_RESID_PC", "40")))
if (!is.finite(MAX_RESID_PC) || MAX_RESID_PC < 0L) MAX_RESID_PC <- 40L
PC_FULL_GRID <- identical(Sys.getenv("PC_FULL_GRID", unset="0"), "1")
PC_FULL_UNIT_DIR <- Sys.getenv("PC_FULL_UNIT_DIR", unset="")
PC_FULL_LEVELS <- if (PC_FULL_GRID) {
  x <- suppressWarnings(as.integer(trimws(strsplit(
    Sys.getenv("PC_FULL_GRID_LEVELS", unset="0,5,10,15,20,25,30,35,40"),
    ",", fixed=TRUE)[[1L]])))
  if (!length(x) || anyNA(x) || any(x < 0L) || any(x %% 5L != 0L) ||
      anyDuplicated(x) || !identical(x, sort(x)) ||
      max(x) > MAX_RESID_PC || !nzchar(PC_FULL_UNIT_DIR)) {
    stop("PC_FULL_GRID requires ordered unique PC levels (multiples of 5), ",
         "MAX_RESID_PC >= max level, and PC_FULL_UNIT_DIR")
  }
  x
} else integer()

to_m_value <- function(beta_mat) {
  epsilon <- 1e-6
  beta_mat[beta_mat < epsilon] <- epsilon
  beta_mat[beta_mat > (1 - epsilon)] <- 1 - epsilon
  log2(beta_mat / (1 - beta_mat))
}

impute_beta_knn <- function(mat, k = 5, rng.seed = 1, block_cols = 50000, maxp = 2000) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  mat[mat < 0] <- 0
  mat[mat > 1] <- 1
  if (!anyNA(mat)) return(mat)

  set.seed(rng.seed)
  n <- nrow(mat)
  p <- ncol(mat)
  out <- mat

  cat(sprintf("  [KNN] Imputing %d features x %d samples in blocks...\n", n, p))

  for (st in seq(1, p, by = block_cols)) {
    ed <- min(st + block_cols - 1, p)
    block <- t(mat[, st:ed, drop = FALSE])

    if (anyNA(block)) {
      row_na <- rowSums(!is.na(block)) == 0
      if (any(row_na)) block[row_na, ] <- 0.5

      col_na <- colSums(!is.na(block)) == 0
      if (any(col_na)) block[, col_na] <- 0.5

      block_imp <- tryCatch({
        impute::impute.knn(block, k=k, rowmax=1, colmax=1, maxp=maxp, rng.seed=rng.seed)$data
      }, error = function(e) {
        cm <- apply(block, 2, median, na.rm=TRUE)
        cm[is.na(cm)] <- 0.5
        inds <- which(is.na(block), arr.ind=TRUE)
        block[inds] <- cm[inds[,2]]
        block
      })

      block_imp[block_imp < 0] <- 0
      block_imp[block_imp > 1] <- 1
      out[, st:ed] <- t(block_imp)
    }
    rm(block)
    if (exists("block_imp")) rm(block_imp)
    gc(verbose = FALSE)
  }

  out
}

compute_lambda500 <- function(pv) {
  pv <- na.omit(pv)
  if (length(pv) == 0) return(NA_real_)
  chi2 <- qchisq(1 - pv, df = 1)
  median(chi2, na.rm=TRUE) / qchisq(0.5, 1)
}

normalize_sex <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[is.na(x) | x == ""] <- "Unknown"
  low <- tolower(x)
  x[low %in% c("f","female")] <- "Female"
  x[low %in% c("m","male")]   <- "Male"
  low2 <- tolower(x)
  x[low2 %in% c("unknown","unk","na","n/a")] <- "Unknown"
  x
}

fit_ols_matrix <- function(y_by_sample, ph, covs) {
  if (!is.matrix(y_by_sample)) y_by_sample <- as.matrix(y_by_sample)
  storage.mode(y_by_sample) <- "double"
  if (!nrow(y_by_sample) || !ncol(y_by_sample)) return(NULL)

  cc <- complete.cases(ph[, covs, drop = FALSE])
  if (!all(cc)) {
    ph <- ph[cc, , drop = FALSE]
    y_by_sample <- y_by_sample[cc, , drop = FALSE]
  }
  if (nrow(ph) < 5L || any(!is.finite(y_by_sample))) return(NULL)
  if (length(unique(ph$outcome)) < 2L) return(NULL)

  design_formula <- as.formula(paste("~", paste(covs, collapse = " + ")))
  design <- model.matrix(design_formula, data = ph)
  outcome_idx <- match("outcome", colnames(design))
  if (is.na(outcome_idx)) return(NULL)

  fit <- tryCatch(lm.fit(x = design, y = y_by_sample), error = function(e) NULL)
  if (is.null(fit) || fit$rank <= 0L || fit$df.residual <= 0) return(NULL)

  coef_mat <- fit$coefficients
  if (is.null(dim(coef_mat))) coef_mat <- matrix(coef_mat, ncol = 1L)
  estimate <- as.numeric(coef_mat[outcome_idx, ])
  if (any(!is.finite(estimate))) {
    bad <- !is.finite(estimate)
  } else {
    bad <- rep(FALSE, length(estimate))
  }

  se <- rep(NA_real_, ncol(y_by_sample))
  p_value <- rep(NA_real_, ncol(y_by_sample))
  piv <- fit$qr$pivot[seq_len(fit$rank)]
  if (outcome_idx %in% piv) {
    rmat <- fit$qr$qr[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
    cov_ind <- tryCatch(chol2inv(rmat), error = function(e) NULL)
    if (!is.null(cov_ind)) {
      outcome_pos <- match(outcome_idx, piv)
      var_outcome <- cov_ind[outcome_pos, outcome_pos]
      rss <- colSums(fit$residuals^2)
      mse <- rss / fit$df.residual
      se <- sqrt(pmax(mse * var_outcome, 0))
      ok <- is.finite(estimate) & is.finite(se) & se > 0
      p_value[ok] <- 2 * stats::pt(abs(estimate[ok] / se[ok]),
                                    df = fit$df.residual, lower.tail = FALSE)
    }
  }

  means <- colMeans(y_by_sample)
  y_sd <- sqrt(pmax((colSums(y_by_sample^2) - nrow(y_by_sample) * means^2) /
                      max(1, nrow(y_by_sample) - 1L), 0))
  list(
    estimate = estimate,
    se = se,
    p_value = p_value,
    y_sd = y_sd,
    invalid = bad
  )
}

fit_ewas_blocks <- function(mat, ph, covs, pheno_type, cpg_ids,
                            cores = 1L, block_size = 5000L) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  cc <- complete.cases(ph[, covs, drop = FALSE])
  if (!all(cc)) {
    ph <- ph[cc, , drop = FALSE]
    mat <- mat[, cc, drop = FALSE]
  }
  if (ncol(mat) < 5L || any(!is.finite(mat))) return(NULL)

  blocks <- split(seq_len(nrow(mat)), ceiling(seq_len(nrow(mat)) / block_size))
  fit_block <- function(idx) {
    stats <- fit_ols_matrix(t(mat[idx, , drop = FALSE]), ph, covs)
    if (is.null(stats)) return(NULL)
    valid <- is.finite(stats$estimate) & is.finite(stats$se) &
      is.finite(stats$p_value) & is.finite(stats$y_sd) & stats$y_sd > 0
    if (!any(valid)) return(NULL)
    sd_x <- if (pheno_type == "binary") {
      mu <- mean(ph$outcome)
      sqrt(mu * (1 - mu))
    } else {
      sd(ph$outcome, na.rm = TRUE)
    }
    std_eff <- if (is.finite(sd_x) && sd_x > 0) stats$estimate * sd_x / stats$y_sd else NA_real_
    std_se <- if (is.finite(sd_x) && sd_x > 0) stats$se * sd_x / stats$y_sd else NA_real_
    data.table(
      cpg = cpg_ids[idx][valid],
      Estimate = stats$estimate[valid],
      SE = stats$se[valid],
      pValue = stats$p_value[valid],
      StdEffectSize = std_eff[valid],
      StdSE = std_se[valid]
    )
  }

  if (.Platform$OS.type != "windows" && cores > 1L && length(blocks) > 1L) {
    pieces <- parallel::mclapply(blocks, fit_block, mc.cores = cores)
  } else {
    pieces <- lapply(blocks, fit_block)
  }
  valid_pieces <- Filter(function(x) is.data.table(x) && nrow(x) > 0L, pieces)
  if (!length(valid_pieces)) return(NULL)
  out <- rbindlist(valid_pieces, fill = TRUE)
  if (!nrow(out)) return(NULL)
  out
}

select_covariates_gse <- function(ph, cov_cells, age_miss_thr = AGE_MISS_THR) {
  cov_cells <- intersect(cov_cells, names(ph))

  keep_age <- FALSE
  age_miss <- NA_real_
  if ("Age" %in% names(ph)) {
    age_miss <- mean(is.na(ph$Age))
    keep_age <- is.finite(age_miss) &&
      age_miss <= age_miss_thr &&
      length(unique(na.omit(ph$Age))) > 1
  }

  keep_sex <- FALSE
  sex_levels <- NULL
  if ("Sex" %in% names(ph)) {
    sex_levels <- sort(unique(na.omit(ph$Sex)))
    if (identical(sort(sex_levels), c("Female","Unknown")) ||
        identical(sort(sex_levels), c("Male","Unknown"))) {
      keep_sex <- FALSE
    } else {
      keep_sex <- length(sex_levels) > 1
    }
  }

  candidates <- character(0)
  if (keep_age) candidates <- c(candidates, "Age")
  if (keep_sex) candidates <- c(candidates, "Sex")
  candidates <- c(candidates, cov_cells)
  candidates <- unique(candidates)

  final_covariates <- Filter(function(v) {
    if (!v %in% names(ph)) return(FALSE)
    vv <- ph[[v]]
    length(unique(na.omit(vv))) > 1
  }, candidates)

  list(
    final_covariates = final_covariates,
    keep_age = keep_age,
    age_missing = age_miss,
    keep_sex = keep_sex,
    sex_levels = sex_levels
  )
}

phenotype <- fread(phenotype_path, data.table=FALSE)
phenotype <- phenotype[!is.na(phenotype$GSE_ID) & phenotype$GSE_ID != "" &
                         !is.na(phenotype[[column_name]]) & phenotype[[column_name]] != "", ]
phenotype[[column_name]] <- trimws(phenotype[[column_name]])
phenotype <- phenotype[phenotype$GSE_ID %in% gse_list, ]

if (pheno_type == "binary") {
  if (is_blank(selected_value)) stop("For binary phenotype, --selected_value must be specified.")

  values <- if (!is_blank(selected_control_value)) {
    c(selected_value, selected_control_value)
  } else {
    c(selected_value, setdiff(unique(phenotype[[column_name]]), selected_value))
  }

  phenotype <- phenotype[phenotype[[column_name]] %in% values, ]
  phenotype$outcome <- ifelse(phenotype[[column_name]] == selected_value, 1, 0)
} else if (pheno_type == "continuous") {
  if (is_blank(subset_column) && !is_blank(selected_value)) {
    subset_column <- "Disease_State2"
    subset_value <- selected_value
    cat(sprintf("[INFO] Continuous mode: using legacy --selected_value as %s == %s\n",
                subset_column, subset_value))
  }

  if (!is_blank(subset_column) || !is_blank(subset_value)) {
    if (is_blank(subset_column) || is_blank(subset_value)) {
      stop("For continuous subset filtering, both --subset_column and --subset_value must be specified.")
    }
    if (!(subset_column %in% names(phenotype))) {
      stop(sprintf("subset_column not found in phenotype table: %s", subset_column))
    }
    phenotype[[subset_column]] <- trimws(as.character(phenotype[[subset_column]]))
    phenotype <- phenotype[phenotype[[subset_column]] == subset_value, ]
  }

  suppressWarnings({ phenotype$outcome <- as.numeric(as.character(phenotype[[column_name]])) })
  phenotype <- phenotype[!is.na(phenotype$outcome), ]
} else {
  stop(sprintf("Unsupported pheno_type: %s", pheno_type))
}

gsm_map <- setNames(
  lapply(gse_list, function(gse) sort(phenotype$GSM_ID[phenotype$GSE_ID == gse])),
  gse_list
)

data(centDHSbloodDMC.m)
data(cent12CT.m)
data(cent12CT450k.m)
data(centUniLIFE.m)
data(centEpiFibIC.m)
data(centEpiFibFatIC.m)
data(centBloodSub.m)

ref_dir <- Sys.getenv("EPISCORE_REF_DIR")
if (!nzchar(ref_dir)) stop("Set EPISCORE_REF_DIR")
ref_files <- list.files(path = ref_dir, pattern = "Ref\\.rda$", full.names = TRUE)
for (file in ref_files) load(file)

tissue2ref <- list(
  whole_blood       = list(ref1=centDHSbloodDMC.m, method="RPC", hier=FALSE),
  whole_blood_EPIC  = list(ref1=cent12CT.m,        method="RPC", hier=FALSE),
  whole_blood_450k  = list(ref1=cent12CT450k.m,    method="RPC", hier=FALSE),
  any_age_blood     = list(ref1=centUniLIFE.m,     method="RPC", hier=FALSE),
  solid_generic     = list(ref1=centEpiFibIC.m,    method="RPC", hier=FALSE),
  solid_with_immune = list(ref1=centEpiFibIC.m, ref2=centBloodSub.m, method="RPC", hier=TRUE, h.CT.idx=which(colnames(centEpiFibIC.m)=="IC")),
  breast            = list(ref1=mrefBreast.m,      method="wRPC", hier=FALSE),
  brain             = list(ref1=mrefBrain.m,       method="wRPC", hier=FALSE),
  prostate          = list(ref1=mrefProstate.m,    method="wRPC", hier=FALSE),
  skin              = list(ref1=mrefSkin.m,        method="wRPC", hier=FALSE),
  lung              = list(ref1=mrefLung.m,        method="wRPC", hier=FALSE),
  bladder           = list(ref1=mrefBladder.m,     method="wRPC", hier=FALSE),
  kidney            = list(ref1=Kidney_Mref.m,     method="wRPC", hier=FALSE),
  liver             = list(ref1=mrefLiver.m,       method="wRPC", hier=FALSE),
  esophagus         = list(ref1=Eso_Mref.m,        method="wRPC", hier=FALSE),
  pancreas          = list(ref1=mrefPancreas.m,    method="wRPC", hier=FALSE),
  colon             = list(ref1=Colon_Mref.m,      method="wRPC", hier=FALSE)
)

if (!(tissue %in% names(tissue2ref))) stop("Unknown tissue.")

compute_cell_props <- function(beta_mat, tissue) {
  nas <- is.na(beta_mat)
  if (any(nas)) {
    row_meds <- apply(beta_mat, 1, median, na.rm=TRUE)
    for (i in seq_len(nrow(beta_mat))) beta_mat[i, nas[i, ]] <- row_meds[i]
  }

  cfg <- tissue2ref[[tissue]]

  if (identical(cfg$method, "CETS")) {
    common <- intersect(rownames(beta_mat), rownames(cfg$profile))
    neu    <- estProportion(beta_mat[common, ], profile=cfg$profile)
    estF   <- t(rbind(Neurons=neu, Glia=1-neu))
  } else if (identical(cfg$method, "wRPC")) {
    beta_prom <- constAvBetaTSS(beta.m = beta_mat, type = "450k")
    common   <- intersect(rownames(beta_prom), rownames(cfg$ref1))
    beta_sub <- beta_prom[common, , drop = FALSE]
    ref_sub  <- cfg$ref1[common,  , drop = FALSE]
    estW     <- wRPC(beta_sub, ref_sub)
    estF     <- estW$estF
  } else if (!isTRUE(cfg$hier)) {
    estF <- epidish(beta.m=beta_mat, ref.m=cfg$ref1, method=cfg$method)$estF
  } else {
    estF <- hepidish(
      beta.m=beta_mat,
      ref1.m=cfg$ref1,
      ref2.m=cfg$ref2,
      h.CT.idx=cfg$h.CT.idx,
      method=cfg$method
    )
  }

  df <- as.data.frame(estF)
  df$GSM_ID <- rownames(df)
  keep <- colSums(df[setdiff(names(df), "GSM_ID")]) > 0
  df[, c("GSM_ID", names(keep)[keep]), drop=FALSE]
}

cell_prop_cache_file <- function(gse, tissue) {
  cache_dir <- Sys.getenv("CELL_PROP_CACHE_DIR", unset = "")
  if (!nzchar(cache_dir)) return("")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_.-]+", "_", gse), "__",
                              gsub("[^A-Za-z0-9_.-]+", "_", tissue), ".cellprops.tsv"))
}

load_or_compute_cell_props <- function(beta_mat, tissue, gse) {
  cache_file <- cell_prop_cache_file(gse, tissue)
  if (!nzchar(cache_file)) return(compute_cell_props(beta_mat, tissue))

  lock_dir <- paste0(cache_file, ".lock")
  while (!dir.create(lock_dir, showWarnings = FALSE)) Sys.sleep(1)
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

  if (file.exists(cache_file) && is.finite(file.info(cache_file)$size) && file.info(cache_file)$size > 0) {
    cached <- tryCatch(fread(cache_file, data.table = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.null(cached) && "GSM_ID" %in% names(cached)) {
      cat(sprintf("[CACHE] Reusing cell proportions for %s/%s from %s\n", gse, tissue, cache_file))
      return(cached)
    }
  }

  cat(sprintf("[CACHE] Computing cell proportions for %s/%s\n", gse, tissue))
  computed <- compute_cell_props(beta_mat, tissue)
  tmp_file <- paste0(cache_file, ".tmp.", Sys.getpid())
  fwrite(computed, tmp_file, sep = "\t")
  if (!file.rename(tmp_file, cache_file)) {
    unlink(tmp_file, force = TRUE)
    if (!file.exists(cache_file)) stop("Failed to write cell-proportion cache: ", cache_file)
  }
  computed
}

load_methylation_data_single <- function(data_dir, gse, gsm_map, data_type, phenotype,
                                         missing_value_threshold=0.7, value_range=c(0,1),
                                         missing_sentinel=NULL) {
  tpl <- switch(data_type,
                beta  = "%s_assay.tsv",
                mvalue = "%s_mvalue.tsv",
                chalm = "%s.chalm.filtered.tsv",
                camda = "%s.camda.filtered.tsv")

  fp <- file.path(data_dir, sprintf(tpl, gse))
  cat(sprintf("[%s] [STEP] Loading %s data for %s from %s\n", Sys.time(), data_type, gse, fp))

  if (!file.exists(fp)) {
    cat("  [WARN] file not found, skip\n")
    return(NULL)
  }

  header <- tryCatch(
    fread(fp, header=TRUE, sep="\t", nrows=0L, data.table=TRUE, showProgress=FALSE),
    error=function(e) NULL
  )
  if (is.null(header) || !ncol(header)) return(NULL)

  orig <- names(header)[-1]

  newn <- vapply(orig, function(x) {
    x <- sub("_.*$", "", x)
    if (grepl("^TCGA", x)) {
      x <- gsub("\\.", "-", x)
    } else if (grepl("^GSE", x)) {
      x <- sub("^[^.]+\\.", "", x)
    }
    x
  }, character(1))
  newn <- gsub("\\.", "-", newn)

  normalized_to_raw <- split(orig, newn)
  candidate_gsms <- intersect(gsm_map[[gse]], newn)
  raw_gsms <- vapply(candidate_gsms, function(x) normalized_to_raw[[x]][1L], character(1))
  if (!length(raw_gsms)) {
    cat("  [WARN] no overlapping GSMs, skip\n")
    return(NULL)
  }

  df <- tryCatch(
    fread(fp, header=TRUE, sep="\t", select=c(names(header)[1L], unname(raw_gsms)),
          data.table=TRUE, showProgress=FALSE),
    error=function(e) NULL
  )
  if (is.null(df)) return(NULL)
  setnames(df, c(names(df)), c("cpg_id", candidate_gsms))

  vals_mat <- as.matrix(df[, ..candidate_gsms])
  suppressWarnings(storage.mode(vals_mat) <- "double")
  vals_mat[is.infinite(vals_mat)] <- NA
  if (!is.null(missing_sentinel)) vals_mat[vals_mat == missing_sentinel] <- NA
  if (!is.null(value_range)) {
    vals_mat[vals_mat < value_range[1]] <- value_range[1]
    vals_mat[vals_mat > value_range[2]] <- value_range[2]
  }
  row_frac <- rowSums(is.na(vals_mat)) / ncol(vals_mat)
  keep <- row_frac < missing_value_threshold
  mat <- vals_mat[keep, , drop=FALSE]
  rownames(mat) <- as.character(df$cpg_id[keep])
  colnames(mat) <- candidate_gsms

  ph <- phenotype[phenotype$GSE_ID==gse & phenotype$GSM_ID %in% candidate_gsms, ]
  ph <- ph[match(candidate_gsms, ph$GSM_ID), ]

  cat(sprintf("  [INFO] Matched samples for %s/%s:\n", gse, data_type))
  if ("outcome" %in% names(ph)) {
    tbl <- table(ph$outcome)
    cat(sprintf("   -> Total: %d | Outcome Counts: %s\n",
                nrow(ph),
                paste(names(tbl), tbl, sep="=", collapse=", ")))
  }

  rm(header, df, vals_mat)
  gc(verbose = FALSE)

  list(methylation_matrix=mat, phenotype=ph)
}

perform_pca_and_lambda_check <- function(mat, ph_dt, matrix_type, pheno_type, cov_cells, gse, output_dir,
                                         fit_cores = 1L) {
  common_samples <- intersect(colnames(mat), ph_dt$GSM_ID)

  if (length(common_samples) < 5) {
    cat(sprintf("  [SKIP] %s (%s): Too few overlapping samples (%d).\n",
                gse, matrix_type, length(common_samples)))
    return(NULL)
  }

  mat0 <- mat[, common_samples, drop=FALSE]
  ph_t <- ph_dt[match(common_samples, ph_dt$GSM_ID), , drop=FALSE]
  rownames(ph_t) <- ph_t$GSM_ID

  valid_outcomes <- na.omit(ph_t$outcome)
  if (length(unique(valid_outcomes)) < 2) {
    cat(sprintf("  [SKIP] %s (%s): Outcome has only %d level(s). Skipping.\n",
                gse, matrix_type, length(unique(valid_outcomes))))
    return(NULL)
  }

  if (nrow(ph_t) < 5) {
    cat(sprintf("  [SKIP] %s (%s): Too few samples (%d). Skipping.\n",
                gse, matrix_type, nrow(ph_t)))
    return(NULL)
  }

  missing_frac <- rowMeans(is.na(mat0) | is.infinite(mat0))
  mat_imp <- mat0[missing_frac <= 0.10, , drop=FALSE]

  nas <- is.na(mat_imp)
  if (any(nas)) {
    medians <- apply(mat_imp, 1, median, na.rm=TRUE)
    medians[is.na(medians)] <- 0
    for (i in seq_len(nrow(mat_imp))) {
      if (any(nas[i, ])) mat_imp[i, nas[i, ]] <- medians[i]
    }
  }

  stopifnot(all(is.finite(mat_imp)))

  suppressWarnings({ ph_t$Age <- as.numeric(as.character(ph_t$Age)) })
  if ("Sex" %in% names(ph_t)) ph_t$Sex <- normalize_sex(ph_t$Sex)

  sel <- select_covariates_gse(ph_t, cov_cells, age_miss_thr = AGE_MISS_THR)

  cat(sprintf("  [DEBUG] %s (%s) Age missing = %s -> %s\n",
              gse, matrix_type,
              ifelse(is.na(sel$age_missing), "NA", sprintf("%.1f%%", 100*sel$age_missing)),
              ifelse(sel$keep_age, "KEEP", "DROP")))

  if (!is.null(sel$sex_levels)) {
    cat(sprintf("  [DEBUG] %s (%s) Sex levels = %s -> %s\n",
                gse, matrix_type,
                paste(sel$sex_levels, collapse=","),
                ifelse(sel$keep_sex, "KEEP", "DROP")))
  } else {
    cat(sprintf("  [DEBUG] %s (%s) Sex not available -> DROP\n", gse, matrix_type))
  }

  final_covariates <- sel$final_covariates

  cat(sprintf("  [DEBUG] %s (%s) Final covariates used for residualization: %s\n",
              gse, matrix_type,
              ifelse(length(final_covariates)==0, "(none)", paste(final_covariates, collapse=", "))))

  keep <- if (length(final_covariates)) {
    complete.cases(ph_t[, final_covariates, drop=FALSE])
  } else {
    rep(TRUE, nrow(ph_t))
  }

  n_before <- nrow(ph_t)
  ph_t     <- ph_t[keep, , drop=FALSE]
  mat_imp  <- mat_imp[, keep, drop=FALSE]
  n_after  <- nrow(ph_t)

  if (n_after < 5) {
    cat(sprintf("  [SKIP] %s (%s): Too few samples after covariate complete.cases (%d -> %d).\n",
                gse, matrix_type, n_before, n_after))
    return(NULL)
  }

  if (length(unique(ph_t$outcome)) < 2) {
    cat(sprintf("  [SKIP] %s (%s): Outcome < 2 levels after covariate cleaning.\n",
                gse, matrix_type))
    return(NULL)
  }

  if (length(final_covariates)) {
    design_known <- model.matrix(~ ., data=ph_t[, final_covariates, drop=FALSE])
    fit_k     <- lm.fit(x = design_known, y = t(mat_imp))
    resid_mat <- t(fit_k$residuals)
  } else {
    resid_mat <- mat_imp
  }

  k_max <- min(MAX_RESID_PC, ncol(resid_mat))
  PC_all <- NULL
  n_pc_available <- 0

  if (k_max > 0) {
    if (length(final_covariates)) {
      sample_crossprod <- crossprod(resid_mat)
    } else {
      row_centered <- sweep(resid_mat, 1L, rowMeans(resid_mat), FUN = "-")
      sample_crossprod <- crossprod(row_centered)
      rm(row_centered)
    }
    eig <- eigen(sample_crossprod, symmetric = TRUE, only.values = FALSE)
    pc_threshold <- if (PC_FULL_GRID) max(eig$values, na.rm=TRUE) * 1e-10 else 0
    positive <- which(is.finite(eig$values) & eig$values > pc_threshold)
    n_pc_available <- min(k_max, length(positive))
    if (n_pc_available > 0L) {
      keep_pc <- positive[seq_len(n_pc_available)]
      PC_all <- eig$vectors[, keep_pc, drop = FALSE] *
        rep(sqrt(eig$values[keep_pc]), each = nrow(eig$vectors))
    }
    rm(sample_crossprod, eig)
    gc(verbose = FALSE)
  }

  ms <- seq(0, n_pc_available, by=5)

  lam_map <- numeric(length(ms))
  ph_merge_list <- vector("list", length(ms))

  for (j in seq_along(ms)) {
    m <- ms[j]

    if (m > 0) {
      pcs <- PC_all[, 1:m, drop=FALSE]
      colnames(pcs) <- paste0("PC_resid", 1:m)
      ph_merge <- cbind(ph_t, pcs)
      covs <- c("outcome", final_covariates, colnames(pcs))
    } else {
      ph_merge <- ph_t
      covs <- c("outcome", final_covariates)
    }

    covs <- intersect(covs, names(ph_merge))

    if (PC_FULL_GRID) {
      ph_merge_list[[j]] <- ph_merge
      next
    }

    pvs <- tryCatch({
      blocks <- split(seq_len(nrow(mat_imp)), ceiling(seq_len(nrow(mat_imp)) / 5000L))
      get_block_p <- function(idx) {
        fit <- fit_ols_matrix(t(mat_imp[idx, , drop = FALSE]), ph_merge, covs)
        if (is.null(fit)) rep(NA_real_, length(idx)) else fit$p_value
      }
      pieces <- if (.Platform$OS.type != "windows" && fit_cores > 1L && length(blocks) > 1L) {
        parallel::mclapply(blocks, get_block_p, mc.cores = fit_cores)
      } else {
        lapply(blocks, get_block_p)
      }
      unlist(pieces, use.names = FALSE)
    }, error = function(e) {
      cat(sprintf("  [WARN] Fast lambda fit failed for %s/%s/m%d; using legacy fit: %s\n",
                  gse, matrix_type, m, conditionMessage(e)))
      apply(mat_imp, 1, function(y) {
        dfm <- data.frame(value = as.numeric(y), ph_merge[, covs, drop=FALSE], check.names=FALSE)
        dfm <- dfm[complete.cases(dfm), , drop=FALSE]
        if (nrow(dfm) < 5 || length(unique(dfm$outcome)) < 2) return(NA_real_)
        rhs <- paste(setdiff(covs, "outcome"), collapse="+")
        fmla <- if (nzchar(rhs)) as.formula(paste("value ~ outcome +", rhs)) else value ~ outcome
        tryCatch(coef(summary(lm(fmla, data=dfm)))["outcome","Pr(>|t|)"], error=function(e) NA_real_)
      })
    })

    lam500 <- compute_lambda500(pvs)
    lam_map[j] <- lam500
    ph_merge_list[[j]] <- ph_merge

    pval_df <- data.frame(cpg_id=rownames(mat_imp), p_value=pvs)
    out_file <- file.path(output_dir, sprintf("%s_%s_m%d_pVal.csv", gse, matrix_type, m))
    fwrite(pval_df, out_file)

    cat(sprintf("  Saved P-values for m=%d to %s (Lambda=%.3f)\n",
                m, out_file, lam500))
  }

  list(lambda = lam_map, ph_merge_list = ph_merge_list, ms = ms)
}

process_and_analyze_gse <- function(gse, gsm_map, phenotype, beta_dir, chalm_dir, camda_dir, tissue, pheno_type,
                                    selected_value, selected_control_value=NULL, cores=8,
                                    mvalue_dir=beta_dir) {
  cat(sprintf("\n===== Processing %s: cell proportions + M/CHALM/CAMDA EWAS =====\n", gse), flush=TRUE)

  gse_pheno <- phenotype[phenotype$GSE_ID == gse, ]
  cat(sprintf("[INFO] Raw Phenotype Counts for %s:\n", gse))
  print(table(gse_pheno$outcome))

  res_b <- load_methylation_data_single(
    data_dir  = beta_dir,
    gse       = gse,
    gsm_map   = gsm_map,
    data_type = "beta",
    phenotype = phenotype
  )

  if (is.null(res_b)) {
    cat(sprintf("[WARN] %s: beta matrix missing/failed. Cannot estimate cell proportions. Skip.\n", gse))
    return(NULL)
  }

  res_b$methylation_matrix <- impute_beta_knn(res_b$methylation_matrix, k=5, rng.seed=1)

  cell_df <- load_or_compute_cell_props(res_b$methylation_matrix, tissue, gse)

  ph_base <- merge(res_b$phenotype, cell_df, by="GSM_ID")

  if (!"Age" %in% names(ph_base)) ph_base$Age <- NA_real_
  suppressWarnings({ ph_base$Age <- as.numeric(as.character(ph_base$Age)) })
  if ("Sex" %in% names(ph_base)) ph_base$Sex <- normalize_sex(ph_base$Sex)

  cat(sprintf("[INFO] %s: beta is used for cell-proportion estimation; M/CHALM/CAMDA EWAS will be run.\n", gse))

  if (identical(Sys.getenv("CELL_PROP_ONLY", unset = "0"), "1")) {
    cat(sprintf("[DONE] Cell-proportion cache prepared for %s/%s\n", gse, tissue))
    return(invisible(TRUE))
  }

  run_method_ewas <- function(method, method_matrix=NULL, already_m_value=FALSE) {
    cat(sprintf("\n----- %s/%s EWAS -----\n", gse, method), flush=TRUE)

    if (is.null(method_matrix)) {
      data_type <- if (method == "m") "mvalue" else method
      data_dir <- switch(method,
                         m = mvalue_dir,
                         chalm = chalm_dir,
                         camda = camda_dir)
      res_m <- load_methylation_data_single(
        data_dir  = data_dir,
        gse       = gse,
        gsm_map   = gsm_map,
        data_type = data_type,
        phenotype = ph_base,
        value_range = if (data_type == "mvalue") NULL else c(0, 1)
      )

      if (is.null(res_m)) {
        cat(sprintf("[WARN] %s: %s matrix missing/failed. Skip %s EWAS.\n",
                    gse, toupper(method), method))
        return(FALSE)
      }

      method_matrix <- res_m$methylation_matrix
      rm(res_m)
    }

    if (!already_m_value) {
      method_matrix <- impute_beta_knn(method_matrix, k=5, rng.seed=1)
      method_matrix <- to_m_value(method_matrix)
    } else {
      storage.mode(method_matrix) <- "double"
    }

    method_res <- perform_pca_and_lambda_check(
      mat         = method_matrix,
      ph_dt       = ph_base,
      matrix_type = method,
      pheno_type  = pheno_type,
      cov_cells   = names(cell_df)[-1],
      gse         = gse,
      output_dir  = output_dir,
      fit_cores   = cores
    )

    if (is.null(method_res)) {
      cat(sprintf("[WARN] %s: %s PCA/lambda check failed. Skip.\n",
                  gse, toupper(method)))
      return(FALSE)
    }

    if (PC_FULL_GRID) {
      unavailable <- setdiff(PC_FULL_LEVELS, method_res$ms)
      if (length(unavailable)) {
        stop("Requested PC level(s) unavailable for ", gse, "/", method,
             ": ", paste(unavailable, collapse=","),
             "; available=", paste(method_res$ms, collapse=","))
      }
      for (pc in PC_FULL_LEVELS) {
        pc_dir <- file.path(PC_FULL_UNIT_DIR, paste0("pc_", pc))
        single_dir <- file.path(pc_dir, "single")
        audit_dir <- file.path(pc_dir, "audit")
        dir.create(single_dir, recursive=TRUE, showWarnings=FALSE)
        dir.create(audit_dir, recursive=TRUE, showWarnings=FALSE)
        outfile <- file.path(single_dir, sprintf(
          "%s_%s_single_cpg_pVal_df.csv", gse, method))
        marker <- file.path(audit_dir, sprintf("%s_%s_single.complete", gse, method))
        if (file.exists(marker) && file.info(marker)$size > 0L &&
            file.exists(outfile) && file.info(outfile)$size > 0L) {
          first_row <- fread(outfile, nrows=1L, check.names=FALSE)
          if (nrow(first_row) == 1L &&
              all(c("cpg", "StdEffectSize", "StdSE", "pValue", "pc_num") %in%
                  names(first_row)) && identical(as.integer(first_row$pc_num[[1L]]),
                                                 as.integer(pc))) {
            cat(sprintf("[REUSE] %s/%s PC=%d full single\n", gse, method, pc))
            next
          }
        }

        ph_dt <- method_res$ph_merge_list[[match(pc, method_res$ms)]]
        common <- intersect(colnames(method_matrix), ph_dt$GSM_ID)
        if (length(common) < 5L) stop("Too few samples for ", gse, "/", method,
                                      " PC=", pc)
        mat2 <- method_matrix[, common, drop=FALSE]
        ph2 <- ph_dt[match(common, ph_dt$GSM_ID), , drop=FALSE]
        suppressWarnings({ ph2$Age <- as.numeric(as.character(ph2$Age)) })
        if ("Sex" %in% names(ph2)) ph2$Sex <- normalize_sex(ph2$Sex)
        sel2 <- select_covariates_gse(ph2, names(cell_df)[-1],
                                      age_miss_thr=AGE_MISS_THR)
        pc_cols <- grep("^PC_resid", names(ph2), value=TRUE)
        if (length(pc_cols) != pc) stop("PC column count mismatch for ", gse,
                                        "/", method, ": requested ", pc,
                                        ", observed ", length(pc_cols))
        covs <- intersect(c("outcome", sel2$final_covariates, pc_cols), names(ph2))
        out <- fit_ewas_blocks(mat=mat2, ph=ph2, covs=covs,
                               pheno_type=pheno_type, cpg_ids=rownames(mat2),
                               cores=cores, block_size=5000L)
        if (is.null(out) || !nrow(out)) stop("No valid full single result for ",
                                            gse, "/", method, " PC=", pc)
        sample_keep <- complete.cases(ph2[, covs, drop=FALSE])
        n_samples <- sum(sample_keep)
        n_cases <- n_controls <- NA_integer_
        if (pheno_type == "binary") {
          n_cases <- sum(ph2$outcome[sample_keep] == 1)
          n_controls <- sum(ph2$outcome[sample_keep] == 0)
        }
        out[, c("n_samples", "n_cases", "n_controls") :=
              list(n_samples, n_cases, n_controls)]
        out[, FDR := p.adjust(pValue, method="BH")]
        out[, pc_num := as.integer(pc)]
        tmp <- tempfile(pattern=paste0(".", basename(outfile), "."),
                        tmpdir=single_dir, fileext=".tmp")
        fwrite(out, tmp)
        if (!file.rename(tmp, outfile)) {
          unlink(tmp)
          stop("Could not publish PC single file: ", outfile)
        }
        audit <- data.table(gse=gse, method=method, pc_requested=pc,
                            pc_used=length(pc_cols), n_samples=nrow(ph2),
                            n_valid_cpg=nrow(out),
                            lambda=compute_lambda500(out$pValue),
                            source_matrix=method,
                            single_file=outfile)
        fwrite(audit, file.path(audit_dir, sprintf(
          "%s_%s_single_audit.csv", gse, method)))
        writeLines(c(paste0("pc_num=", pc), paste0("n_cpg=", nrow(out)),
                     paste0("file=", outfile)), marker)
        cat(sprintf("[DONE] %s/%s PC=%d full single: %d CpGs\n",
                    gse, method, pc, nrow(out)))
        rm(out, mat2)
        gc(verbose=FALSE)
      }
      rm(method_res)
      gc(verbose=FALSE)
      return(TRUE)
    }

    lambda_df <- data.frame(PC_num = method_res$ms, lambda = method_res$lambda)
    names(lambda_df)[2] <- method

    cat(sprintf("\n======= %s Lambda =======\n", toupper(method)))
    print(lambda_df)
    cat("\n")

    fwrite(lambda_df, file.path(output_dir, sprintf("%s_%s_lambda_table.csv", gse, method)))

    if (num_fixed >= 0) {
      if (num_fixed %in% lambda_df$PC_num) {
        best_m <- num_fixed
      } else {
        best_m <- lambda_df$PC_num[which.min(abs(lambda_df$PC_num - num_fixed))]
        cat(sprintf("[WARN] --num=%d not in calculated steps for %s; using nearest PC_num=%d\n",
                    num_fixed, method, best_m))
      }

      cat(sprintf("Using user-specified %s PC number: %d\n", toupper(method), best_m))
    } else {
      lambda_values <- lambda_df[[method]]
      ok <- lambda_values < 1.3
      ok[is.na(ok)] <- FALSE

      if (sum(ok) >= 1) {
        best_m <- min(lambda_df$PC_num[ok])
      } else {
        valid_idx <- which(!is.na(lambda_values))
        best_m <- if (length(valid_idx) > 0) {
          lambda_df$PC_num[valid_idx[which.min(lambda_values[valid_idx])]]
        } else {
          0
        }
      }

      cat(sprintf("Auto-selected %s PC number: %d\n", toupper(method), best_m))
    }

    ph_dt <- method_res$ph_merge_list[[which(method_res$ms == best_m)]]
    rm(method_res)
    gc(verbose = FALSE)
    mat2  <- method_matrix

    if (is.null(ph_dt) || is.null(mat2)) {
      cat(sprintf("[WARN] %s/%s: phenotype or matrix missing for final EWAS. Skip.\n", gse, method))
      return(FALSE)
    }

    common <- intersect(colnames(mat2), ph_dt$GSM_ID)

    if (length(common) < 5) {
      cat(sprintf("[WARN] %s/%s: too few overlapping samples for final EWAS (%d). Skip.\n",
                  gse, method, length(common)))
      return(FALSE)
    }

    mat2 <- mat2[, common, drop=FALSE]
    ph2  <- ph_dt[match(common, ph_dt$GSM_ID), , drop=FALSE]

    suppressWarnings({ ph2$Age <- as.numeric(as.character(ph2$Age)) })
    if ("Sex" %in% names(ph2)) ph2$Sex <- normalize_sex(ph2$Sex)

    pc_cols <- grep("^PC_resid", names(ph2), value=TRUE)

    sel2 <- select_covariates_gse(ph2, names(cell_df)[-1], age_miss_thr = AGE_MISS_THR)

    cat(sprintf("[DEBUG] %s/%s final EWAS: Age missing=%s -> %s | Sex levels=%s -> %s\n",
                gse,
                method,
                ifelse(is.na(sel2$age_missing), "NA", sprintf("%.1f%%", 100*sel2$age_missing)),
                ifelse(sel2$keep_age, "KEEP", "DROP"),
                ifelse(is.null(sel2$sex_levels), "NA", paste(sel2$sex_levels, collapse=",")),
                ifelse(sel2$keep_sex, "KEEP", "DROP")))

    covs <- c("outcome", sel2$final_covariates, pc_cols)
    covs <- intersect(covs, names(ph2))

    out <- fit_ewas_blocks(
      mat = mat2,
      ph = ph2,
      covs = covs,
      pheno_type = pheno_type,
      cpg_ids = rownames(mat2),
      cores = cores,
      block_size = 5000L
    )

    if (!is.null(out) && nrow(out) > 0L) {
      sample_keep <- complete.cases(ph2[, covs, drop=FALSE])
      n_samples <- sum(sample_keep)
      n_cases <- n_controls <- NA_integer_
      if (pheno_type == "binary") {
        n_cases <- sum(ph2$outcome[sample_keep] == 1)
        n_controls <- sum(ph2$outcome[sample_keep] == 0)
      }
      out[, c("n_samples", "n_cases", "n_controls") :=
            list(n_samples, n_cases, n_controls)]
      out[, FDR := p.adjust(pValue, method="fdr")]

      fwrite(out, file.path(output_dir, sprintf("%s_%s_single_cpg_pVal_df.csv", gse, method)))

      cat(sprintf("[INFO] Wrote %s result for %s (n=%d)\n", toupper(method), gse, nrow(out)))
      return(TRUE)
    }

    cat(sprintf("[INFO] No valid CpGs produced for %s/%s after QC.\n", gse, method))
    FALSE
  }

  requested_methods <- trimws(strsplit(Sys.getenv("REVISION_METHODS", unset=""), "[:,]", perl=TRUE)[[1]])
  requested_methods <- requested_methods[nzchar(requested_methods)]
  if (!length(requested_methods)) requested_methods <- c("m", "chalm", "camda")
  allowed_methods <- c("m", "chalm", "camda")
  requested_methods <- unique(requested_methods[requested_methods %in% allowed_methods])
  if (!length(requested_methods)) stop("Choose at least one method from m, chalm and camda")

  method_status <- list()
  for (method in intersect(requested_methods, c("chalm", "camda"))) {
    method_status[[method]] <- run_method_ewas(method)
  }
  if ("m" %in% requested_methods) {
    method_status[["m"]] <- run_method_ewas("m", already_m_value=TRUE)
  }
  method_status <- unlist(method_status)

  if (!all(method_status)) {
    failed_methods <- names(method_status)[!method_status]
    cat(sprintf("[WARN] %s: method(s) failed or produced no output: %s\n",
                gse, paste(failed_methods, collapse=", ")))
  }

  return(invisible(method_status))

  res_d <- load_methylation_data_single(
    data_dir  = camda_dir,
    gse       = gse,
    gsm_map   = gsm_map,
    data_type = "camda",
    phenotype = ph_base
  )

  if (is.null(res_d)) {
    cat(sprintf("[WARN] %s: CAMDA matrix missing/failed. Skip.\n", gse))
    return(NULL)
  }

  res_d$methylation_matrix <- impute_beta_knn(res_d$methylation_matrix, k=5, rng.seed=1)

  res_d$methylation_matrix <- to_m_value(res_d$methylation_matrix)

  camda_res <- perform_pca_and_lambda_check(
    mat         = res_d$methylation_matrix,
    ph_dt       = ph_base,
    matrix_type = "camda",
    pheno_type  = pheno_type,
    cov_cells   = names(cell_df)[-1],
    gse         = gse,
    output_dir  = output_dir
  )

  if (is.null(camda_res)) {
    cat(sprintf("[WARN] %s: CAMDA PCA/lambda check failed. Skip.\n", gse))
    return(NULL)
  }

  lambda_df <- data.frame(
    PC_num = camda_res$ms,
    camda  = camda_res$lambda
  )

  cat("\n======= CAMDA Lambda =======\n")
  print(lambda_df)
  cat("\n")

  fwrite(lambda_df, file.path(output_dir, sprintf("%s_camda_lambda_table.csv", gse)))

  if (num_fixed >= 0) {
    if (num_fixed %in% lambda_df$PC_num) {
      best_m <- num_fixed
    } else {
      best_m <- lambda_df$PC_num[which.min(abs(lambda_df$PC_num - num_fixed))]
      cat(sprintf("[WARN] --num=%d not in calculated steps; using nearest PC_num=%d\n", num_fixed, best_m))
    }

    cat(sprintf("Using user-specified CAMDA PC number: %d\n", best_m))
  } else {
    ok <- lambda_df$camda < 1.3
    ok[is.na(ok)] <- FALSE

    if (sum(ok) >= 1) {
      best_m <- min(lambda_df$PC_num[ok])
    } else {
      valid_idx <- which(!is.na(lambda_df$camda))
      best_m <- if (length(valid_idx) > 0) {
        lambda_df$PC_num[valid_idx[which.min(lambda_df$camda[valid_idx])]]
      } else {
        0
      }
    }

    cat(sprintf("Auto-selected CAMDA PC number: %d\n", best_m))
  }

  ph_dt <- camda_res$ph_merge_list[[which(camda_res$ms == best_m)]]
  mat2  <- res_d$methylation_matrix

  if (is.null(ph_dt) || is.null(mat2)) {
    cat(sprintf("[WARN] %s/camda: phenotype or matrix missing for final EWAS. Skip.\n", gse))
    return(NULL)
  }

  common <- intersect(colnames(mat2), ph_dt$GSM_ID)

  if (length(common) < 5) {
    cat(sprintf("[WARN] %s/camda: too few overlapping samples for final EWAS (%d). Skip.\n",
                gse, length(common)))
    return(NULL)
  }

  mat2 <- mat2[, common, drop=FALSE]
  ph2  <- ph_dt[match(common, ph_dt$GSM_ID), , drop=FALSE]

  suppressWarnings({ ph2$Age <- as.numeric(as.character(ph2$Age)) })
  if ("Sex" %in% names(ph2)) ph2$Sex <- normalize_sex(ph2$Sex)

  pc_cols <- grep("^PC_resid", names(ph2), value=TRUE)

  sel2 <- select_covariates_gse(ph2, names(cell_df)[-1], age_miss_thr = AGE_MISS_THR)

  cat(sprintf("[DEBUG] %s/camda final EWAS: Age missing=%s -> %s | Sex levels=%s -> %s\n",
              gse,
              ifelse(is.na(sel2$age_missing), "NA", sprintf("%.1f%%", 100*sel2$age_missing)),
              ifelse(sel2$keep_age, "KEEP", "DROP"),
              ifelse(is.null(sel2$sex_levels), "NA", paste(sel2$sex_levels, collapse=",")),
              ifelse(sel2$keep_sex, "KEEP", "DROP")))

  covs <- c("outcome", sel2$final_covariates, pc_cols)
  covs <- intersect(covs, names(ph2))

  results <- mclapply(seq_len(nrow(mat2)), function(i) {
    y <- as.numeric(mat2[i, ])
    if (all(is.na(y))) return(NULL)

    dfm <- data.frame(cpg = y, ph2[, covs, drop=FALSE], check.names=FALSE)

    dfm <- dfm[complete.cases(dfm), , drop=FALSE]

    if (nrow(dfm) < 5) return(NULL)

    y_sd <- sd(dfm$cpg, na.rm = TRUE)
    if (is.na(y_sd) || y_sd == 0) return(NULL)

    covs_i <- Filter(function(v) {
      v %in% names(dfm) && length(unique(na.omit(dfm[[v]]))) > 1
    }, covs)

    if (!("outcome" %in% covs_i)) return(NULL)

    rhs <- paste(setdiff(covs_i, "outcome"), collapse="+")
    fmla_i <- if (nzchar(rhs)) {
      as.formula(paste("cpg ~ outcome +", rhs))
    } else {
      cpg ~ outcome
    }

    fit <- tryCatch(lm(fmla_i, data=dfm), error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    sm <- summary(fit)
    if (!("outcome" %in% rownames(sm$coefficients))) return(NULL)

    cm <- sm$coefficients["outcome",]
    est <- unname(cm["Estimate"])
    se  <- unname(cm["Std. Error"])
    p   <- unname(cm["Pr(>|t|)"])

    sd_x <- if (pheno_type=="binary") {
      mu <- mean(dfm$outcome)
      sqrt(mu*(1-mu))
    } else {
      sd(dfm$outcome, na.rm=TRUE)
    }

    std_eff <- if (!is.na(sd_x) && sd_x>0 && y_sd>0) est*(sd_x/y_sd) else NA_real_
    std_se  <- if (!is.na(sd_x) && sd_x>0 && y_sd>0) se*(sd_x/y_sd) else NA_real_

    data.frame(
      cpg=rownames(mat2)[i],
      Estimate=est,
      SE=se,
      pValue=p,
      StdEffectSize=std_eff,
      StdSE=std_se,
      stringsAsFactors=FALSE
    )
  }, mc.cores=cores)

  res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, results)

  if (length(res_list) > 0) {
    out <- rbindlist(res_list, fill=TRUE)
    out[, FDR := p.adjust(pValue, method="fdr")]

    fwrite(out, file.path(output_dir, sprintf("%s_camda_single_cpg_pVal_df.csv", gse)))

    cat(sprintf("[INFO] Wrote CAMDA result for %s (n=%d)\n", gse, nrow(out)))
  } else {
    cat(sprintf("[INFO] No valid CpGs produced for %s/camda after QC.\n", gse))
  }
}

for (gse in gse_list) {
  process_and_analyze_gse(
    gse = gse,
    gsm_map = gsm_map,
    phenotype = phenotype,
    beta_dir = beta_dir,
    mvalue_dir = mvalue_dir,
    chalm_dir = chalm_dir,
    camda_dir = camda_dir,
    tissue = tissue,
    pheno_type = pheno_type,
    selected_value = selected_value,
    selected_control_value = selected_control_value,
    cores = cores
  )
}
