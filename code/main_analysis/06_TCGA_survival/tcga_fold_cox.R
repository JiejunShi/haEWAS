# Run fold-specific TCGA Cox scans.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(survival)
  library(EpiDISH)
  library(EpiSCORE)
  library(cets)
  library(parallel)
})
option_list <- list(
  make_option(c("--projects"), type="character", default=NULL,
              help="One TCGA project ID"),
  make_option(c("--outdir"), type="character", default=file.path(tempdir(), "tcga_survival"),
              help="Base output directory [default %default]"),
  make_option(c("--cores"), type="integer", default=2,
              help="Cores for parallel Cox [default %default]"),
  make_option(c("--tissue"), type="character", default=NULL,
              help="Override tissue for cell proportion (optional). If not set, auto map by project."),
  make_option(c("--survival_phenotype"), type="character", default=NULL,
              help="Merged survival phenotype used to define the sample rows"),
  make_option(c("--camda_dir"), type="character",
              default=NULL, help="Updated CAMDA matrix directory"),
  make_option(c("--phenotype_dir"), type="character", default=NULL,
              help="GDC phenotype directory"),
  make_option(c("--beta_dir"), type="character", default=NULL,
              help="Beta matrix directory"),
  make_option(c("--chalm_dir"), type="character", default=NULL,
              help="CHALM matrix directory"),
  make_option(c("--ref_dir"), type="character", default=NULL,
              help="EpiSCORE reference directory")
)
opts <- parse_args(OptionParser(option_list=option_list))
required_dirs <- list(opts$camda_dir, opts$phenotype_dir, opts$beta_dir,
                      opts$chalm_dir, opts$ref_dir)
if (any(vapply(required_dirs, function(x) is.null(x) || !dir.exists(x), logical(1)))) {
  stop("Provide existing CAMDA, phenotype, beta, CHALM and reference directories")
}
if (is.null(opts$survival_phenotype) || !file.exists(opts$survival_phenotype)) {
  stop("Merged survival phenotype file is required")
}
base_out_dir  <- opts$outdir
dir.create(base_out_dir, showWarnings=FALSE, recursive=TRUE)
tcga_pheno_dir   <- opts$phenotype_dir
beta_dir         <- opts$beta_dir
chalm_camda_dir  <- opts$chalm_dir
drop_sex_projects  <- c("BRCA","PRAD")
norm_barcode <- function(x, keep_first_16 = FALSE) {
  y <- toupper(gsub("\\.", "-", trimws(as.character(x))))
  if (keep_first_16) {
    y <- sub("^([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]{3}).*$", "\\1", y)
  }
  y
}

tcga_stage_to_factor <- function(stage_chr) {
  s <- tolower(trimws(as.character(stage_chr)))
  s <- gsub("[_ -]+", " ", s)
  s <- sub("^stage ", "", s); s <- sub("^stage:", "", s)
  s[s %in% c("not reported","unknown","x","tx","nx","na","n/a","")] <- NA
  map_roman <- function(z){
    if (is.na(z)) return(NA_integer_)
    z <- gsub("[abcduv]+$","", z); z <- trimws(z)
    if (z %in% c("i","1"))   return(1L)
    if (z %in% c("ii","2"))  return(2L)
    if (z %in% c("iii","3")) return(3L)
    if (z %in% c("iv","4"))  return(4L)
    NA_integer_
  }
  main <- sub("^(i{1,3}|iv|[1-4]).*$", "\\1", s)
  v <- vapply(main, map_roman, 0L)
  factor(ifelse(is.na(v), "Unknown", paste0("S", v)),
         levels = c("S1","S2","S3","S4","Unknown"))
}

.map_sample_type <- function(x, warn_unknown = TRUE) {
  s <- tolower(trimws(as.character(x)))

  fine <- ifelse(s == "primary tumor",    "Primary Tumor",
          ifelse(s == "recurrent tumor",  "Recurrent Tumor",
          ifelse(s == "solid tissue normal","Solid Tissue Normal", NA_character_)))
  fine[is.na(fine) & grepl("normal", s)]                    <- "Solid Tissue Normal"
  fine[is.na(fine) & grepl("tumou?r|metastatic|cancer", s)] <- "Other Tumor"
  fine[is.na(fine)] <- "Other"
  fine <- factor(fine, levels = c("Primary Tumor","Recurrent Tumor","Solid Tissue Normal","Other Tumor","Other"))

  group <- ifelse(fine == "Primary Tumor", "Tumor",
           ifelse(fine == "Solid Tissue Normal", "Normal", "Other"))
  group <- factor(group, levels = c("Tumor","Normal","Other"))

  if (isTRUE(warn_unknown)) {
    known_exact <- c("primary tumor","recurrent tumor","solid tissue normal")
    covered_by_regex <- grepl("normal|tumou?r|metastatic|cancer", s)
    unknown_raw <- sort(unique(s[!(s %in% known_exact) & !covered_by_regex & nzchar(s)]))
    if (length(unknown_raw)) message("[WARN] Unmapped sample_type.samples values: ", paste(unknown_raw, collapse = ", "))
  }

  list(fine = fine, group = group)
}

read_gdc_phenotype <- function(project_code, pheno_dir = tcga_pheno_dir, keep_first_16 = FALSE) {
  fp <- file.path(pheno_dir, sprintf("%s.GDC_phenotype.tsv", project_code))
  if (!file.exists(fp)) stop("[ERR] GDC phenotype not found: ", fp)
  ph <- data.table::fread(fp, sep="\t", data.table=FALSE, check.names=FALSE)

  getcol <- function(nm) if (nm %in% names(ph)) ph[[nm]] else rep(NA, nrow(ph))
  sample  <- norm_barcode(getcol("submitter_id.samples"), keep_first_16)
  age     <- suppressWarnings(as.numeric(getcol("age_at_index.demographic")))
  sex     <- tolower(trimws(getcol("gender.demographic")))
  vital   <- tolower(trimws(getcol("vital_status.demographic")))
  dlfu    <- suppressWarnings(as.numeric(getcol("days_to_last_follow_up.diagnoses")))
  ddeath  <- if ("days_to_death.demographic" %in% names(ph))
               suppressWarnings(as.numeric(getcol("days_to_death.demographic"))) else rep(NA, nrow(ph))
  stage_t <- getcol("tumor_stage.diagnoses")
  samp_type <- getcol("sample_type.samples")

  time_days  <- ifelse(is.finite(ddeath), ddeath, dlfu)
  time_years <- time_days / 365.25
  status     <- ifelse(vital == "dead", 1L, ifelse(vital == "alive", 0L, NA_integer_))

  Sex_fac <- factor(sub("^male$","Male", sub("^female$","Female", sex)),
                    levels = c("Male","Female"))
  Sex_fac <- addNA(Sex_fac); levels(Sex_fac)[is.na(levels(Sex_fac))] <- "Unknown"

  Stage_fac <- tcga_stage_to_factor(stage_t)
  mt <- .map_sample_type(samp_type)

  out <- data.frame(
    Sample              = sample,
    Age                 = age,
    Sex_fac             = factor(Sex_fac, levels=c("Male","Female","Unknown")),
    Stage_fac           = Stage_fac,
    time                = time_years,
    status              = status,
    sample_type.samples = samp_type,
    SampleType_fac      = mt$fine,
    SampleType_group    = mt$group,
    stringsAsFactors = FALSE
  )
  stage_value <- as.character(out$Stage_fac)
  sex_value <- as.character(out$Sex_fac)
  keep <- is.finite(out$time) & out$time > 0 & !is.na(out$status) & is.finite(out$Age) &
    sex_value %in% c("Male", "Female") & stage_value %in% paste0("S", 1:4)
  out <- out[keep, ]
  out
}
load_tcga_methyl_one <- function(project_id,
                                 data_type = c("beta","chalm","camda"),
                                 ph_tumor,
                                 beta_dir_in = NULL,
                                 chalm_camda_dir_in = NULL,
                                 na_row_thresh = 0.7,
                                 na_col_thresh = 0.7) {

  if (is.null(beta_dir_in))        beta_dir_in        <- get0("beta_dir",        inherits = TRUE)
  if (is.null(chalm_camda_dir_in)) chalm_camda_dir_in <- get0("chalm_camda_dir", inherits = TRUE)

  data_type <- match.arg(data_type)
  tpl <- switch(data_type,
                beta  = "%s_assay.tsv",
                chalm = "%s.chalm.filtered.tsv",
                camda = "%s.camda.filtered.tsv")
  data_dir <- if (data_type == "beta") beta_dir_in else chalm_camda_dir_in
  fp <- file.path(data_dir, sprintf(tpl, project_id))

  cat("[STEP] Loading", data_type, "for", project_id, "from", fp, "\n")
  if (!file.exists(fp)) { cat("  [SKIP] file not found\n"); return(NULL) }

  df <- data.table::fread(fp, data.table = FALSE, check.names = FALSE)
  colnames(df) <- gsub("\\.", "-", colnames(df))
  names(df)[1] <- "cpg_id"
  if (!"cpg_id" %in% names(df)) { cat("  [ERR] first column must be CpG IDs\n"); return(NULL) }

  vals <- as.data.frame(lapply(df[, -1, drop = FALSE], as.numeric), check.names = FALSE)
  m <- as.matrix(vals)
  m[!is.finite(m)] <- NA
  m[m < 0] <- 0
  m[m > 1] <- 1

  rownames(m) <- as.character(df$cpg_id)
  dup <- duplicated(rownames(m))
  if (any(dup)) {
    m <- m[!dup, , drop = FALSE]
    cat("  [INFO] removed", sum(dup), "duplicated CpGs\n")
  }

  if (nrow(m) > 0) {
    row_frac <- rowSums(is.na(m)) / ncol(m)
    keep_r   <- row_frac < na_row_thresh
    if (sum(keep_r) == 0) { cat("  [WARN] no rows passed NA filter\n"); return(NULL) }
    if (sum(keep_r) < nrow(m)) cat("  [INFO] kept", sum(keep_r), "/", nrow(m), "CpGs after NA-row filter\n")
    m <- m[keep_r, , drop = FALSE]
  }
  if (nrow(m) == 0) { cat("  [WARN] no CpGs left after filtering\n"); return(NULL) }

  col_frac <- colSums(is.na(m)) / nrow(m)
  keep_c   <- col_frac < na_col_thresh
  if (sum(keep_c) == 0) { cat("  [WARN] no samples passed NA-col filter\n"); return(NULL) }
  if (sum(keep_c) < ncol(m)) cat("  [INFO] kept", sum(keep_c), "/", ncol(m), "samples after NA-col filter\n")
  m <- m[, keep_c, drop = FALSE]
  tumor_ids <- unique(as.character(ph_tumor$Sample))
  keep_ids  <- intersect(colnames(m), tumor_ids)

  if (length(keep_ids) < 5) {
    cat("  [WARN]", project_id, data_type, "too few samples after intersection:", length(keep_ids), "\n")
    return(NULL)
  }

  m_keep <- m[, keep_ids, drop = FALSE]
  ph_keep <- ph_tumor[match(keep_ids, ph_tumor$Sample), , drop = FALSE]
  stopifnot(identical(colnames(m_keep), ph_keep$Sample))

  cat(sprintf("  [OK] %s matrix: %d CpGs x %d samples (after filters & intersections)\n",
              data_type, nrow(m_keep), ncol(m_keep)))

  list(mat = m_keep, ph = ph_keep)
}
data(centDHSbloodDMC.m); data(cent12CT.m); data(cent12CT450k.m)
data(centUniLIFE.m); data(centEpiFibIC.m); data(centBloodSub.m)

ref_dir <- opts$ref_dir
ref_files <- list.files(path = ref_dir, pattern = "Ref\\.rda$", full.names = TRUE)
for (file in ref_files) load(file)

tissue2ref <- list(
  solid_generic     = list(ref1=centEpiFibIC.m,    method="RPC",  hier=FALSE),
  brain             = list(ref1=mrefBrain.m,       method="wRPC", hier=FALSE),
  breast            = list(ref1=mrefBreast.m,      method="wRPC", hier=FALSE),
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

project2tissue <- function(project_code) {
  short <- sub("^TCGA-", "", project_code)
  if (short %in% c("GBM","LGG")) return("brain")
  if (short %in% c("BRCA"))     return("breast")
  if (short %in% c("PRAD"))     return("prostate")
  if (short %in% c("SKCM"))     return("skin")
  if (short %in% c("LUAD","LUSC")) return("lung")
  if (short %in% c("BLCA"))     return("bladder")
  if (short %in% c("KIRC","KIRP","KICH")) return("kidney")
  if (short %in% c("LIHC"))     return("liver")
  if (short %in% c("PAAD"))     return("pancreas")
  if (short %in% c("COAD","READ")) return("colon")
  if (short %in% c("ESCA"))     return("esophagus")
  "solid_generic"
}

make_folds <- function(samples, ph, output_dir) {
  if (anyDuplicated(ph$GSM_ID)) stop("Duplicated survival sample IDs")
  ph <- ph[match(samples, ph$GSM_ID), , drop=FALSE]
  if (anyNA(ph$GSM_ID)) stop("Survival phenotype does not cover the beta sample rows")

  rows <- list()
  for (repeat_id in 1:3) {
    set.seed(41 + repeat_id)
    order <- sample(seq_along(samples))
    sizes <- rep(length(samples) %/% 5, 5)
    sizes[seq_len(length(samples) %% 5)] <- sizes[seq_len(length(samples) %% 5)] + 1L
    first <- 1L
    for (fold in 1:5) {
      test <- order[first:(first + sizes[fold] - 1L)]
      scope <- paste0("r", repeat_id, "_f", fold)
      rows[[length(rows) + 1L]] <- data.frame(
        scope=scope, repeat_id=repeat_id, outer_fold=fold, Sample=samples,
        partition=ifelse(seq_along(samples) %in% test, "test", "train"),
        time=ph$time, status=ph$status, stringsAsFactors=FALSE
      )
      first <- first + sizes[fold]
    }
  }
  assignments <- rbindlist(rows)
  path <- file.path(output_dir, "folds", "split_assignments.csv")
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  fwrite(assignments, path)
  assignments
}

compute_cell_props <- function(beta_mat, tissue_key) {
  cfg <- tissue2ref[[tissue_key]]
  if (is.null(cfg)) cfg <- tissue2ref[["solid_generic"]]

  nas <- is.na(beta_mat)
  if (any(nas)) {
    row_meds <- apply(beta_mat, 1, median, na.rm=TRUE)
    for (i in seq_len(nrow(beta_mat))) beta_mat[i, nas[i, ]] <- row_meds[i]
  }

  if (identical(cfg$method, "wRPC")) {
    beta_prom <- constAvBetaTSS(beta.m=beta_mat, type="450k")
    common    <- intersect(rownames(beta_prom), rownames(cfg$ref1))
    estW      <- wRPC(beta_prom[common,,drop=FALSE], cfg$ref1[common,,drop=FALSE])
    estF      <- estW$estF
  } else {
    estF <- epidish(beta.m=beta_mat, ref.m=cfg$ref1, method=cfg$method)$estF
  }
  df <- as.data.frame(estF); df$Sample <- rownames(df); df
}

celltype_map <- list(
  IC  = c("IC","Macro","Lym","Gran","Mono","Plasma","Tcell","Mye"),
  EC  = c("Endo"),
  Fib = c("Fib","Stellate","CM","SMC","Fat","Stromal")
)

derive_purity_immune <- function(cell_df, map) {
  avail <- setdiff(names(cell_df), "Sample")
  pc <- function(keys, parent=NULL){
    cols <- intersect(keys, avail)
    if (!is.null(parent)) cols <- setdiff(cols, parent)
    cols
  }
  ic  <- pc(map$IC,  "IC");  ec <- pc(map$EC);  fib <- pc(map$Fib, "Fib")
  Immune_total <- if (length(ic))  rowSums(cell_df[, ic,  drop=FALSE], na.rm=TRUE) else 0
  Endothelial  <- if (length(ec))  rowSums(cell_df[, ec,  drop=FALSE], na.rm=TRUE) else 0
  Fibro        <- if (length(fib)) rowSums(cell_df[, fib, drop=FALSE], na.rm=TRUE) else 0
  Purity <- 1 - (Immune_total + Endothelial + Fibro)
  Purity[Purity<0] <- 0; Purity[Purity>1] <- 1
  out <- data.frame(Sample=cell_df$Sample, Immune_total=Immune_total, Purity=Purity, stringsAsFactors=FALSE)
  out
}
if (is.null(opts$projects) || !nzchar(trimws(opts$projects))) stop("Provide one project ID")
pick_projects <- unique(trimws(strsplit(opts$projects, ",", fixed=TRUE)[[1]]))
if (length(pick_projects) != 1L) stop("Run one project per call")
for (proj in pick_projects) {

  cat("\n========== PROCESSING Project:", proj, "==========\n")
  ph_all <- read_gdc_phenotype(proj, pheno_dir = tcga_pheno_dir, keep_first_16 = FALSE)
  tumor_ids <- ph_all$Sample[ph_all$SampleType_group == "Tumor"]
  ph_tumor  <- ph_all[ph_all$Sample %in% tumor_ids, , drop=FALSE]
  cat("[INFO]", proj, "tumor phenotype n =", nrow(ph_tumor), "\n")
  if (nrow(ph_tumor) < 10) { cat("  [SKIP] too few tumor samples\n"); next }
  res_beta <- load_tcga_methyl_one(proj, "beta", ph_tumor)
  if (is.null(res_beta)) { cat("  [WARN] no BETA usable; skip project\n"); next }

  tissue_key <- if (!is.null(opts$tissue)) opts$tissue else project2tissue(proj)
  cat("[INFO]", proj, "tissue =", tissue_key, "\n")

  cell_df <- compute_cell_props(res_beta$mat, tissue_key)
  colnames(cell_df)[colnames(cell_df)=="Sample"] <- "Sample"
  pur_df  <- derive_purity_immune(cell_df, celltype_map)
  ph0 <- merge(res_beta$ph, cell_df, by.x="Sample", by.y="Sample")
  ph0 <- merge(ph0, pur_df, by="Sample")
  rownames(ph0) <- ph0$Sample
  cat("[INFO] merged pheno with cell props; n =", nrow(ph0), "\n")
  short <- sub("^TCGA-","",proj)
  include_sex <- !(short %in% drop_sex_projects)
  if (!include_sex) cat("[INFO] in sex project blacklist → drop Sex\n")
  use_stage <- length(unique(as.character(ph0$Stage_fac))) > 1L
  if (!use_stage) cat("[INFO] Stage has no variation → omit constant term from Cox\n")
  samples <- colnames(res_beta$mat)
  surv <- fread(opts$survival_phenotype, data.table=FALSE, check.names=FALSE)
  surv$time <- suppressWarnings(as.numeric(surv$`Overall_Survival(Year)`))
  vital <- tolower(trimws(as.character(surv$Vital_Status)))
  surv$status <- ifelse(vital == "dead", 1L, ifelse(vital == "alive", 0L, NA_integer_))
  surv$GSM_ID <- norm_barcode(surv$GSM_ID)
  surv <- surv[!is.na(surv$GSE_ID) & surv$GSE_ID == proj &
               is.finite(surv$time) & surv$time > 0 & !is.na(surv$status), , drop=FALSE]
  if (anyDuplicated(surv$GSM_ID)) stop("Duplicated survival sample IDs")
  samples <- samples[samples %in% surv$GSM_ID]
  if (length(samples) < 10L) stop("Too few beta samples with survival data")
  fold_assignments <- make_folds(samples, surv, base_out_dir)
  fold_scopes <- unique(as.character(fold_assignments$scope))
  for (dt in c("beta","chalm","camda")) {
    res_dt <- if (dt=="beta") res_beta else if (dt=="camda")
      load_tcga_methyl_one(proj, dt, ph_tumor,
                           chalm_camda_dir_in=opts$camda_dir) else
      load_tcga_methyl_one(proj, dt, ph_tumor)
    if (is.null(res_dt)) { cat("  [WARN] no data for", dt, "skip\n"); next }
    mat_full <- res_dt$mat
    if (!all(samples %in% colnames(mat_full))) stop("Assay sample rows differ for ", dt)
    for (scope in fold_scopes) {
    train_ids <- fold_assignments$Sample[fold_assignments$scope==scope &
                                         fold_assignments$partition=="train"]
    if (!length(train_ids)) stop("Empty train set: ", scope)
    keep <- which(colnames(mat_full) %in% train_ids)
    mat <- mat_full[, keep, drop=FALSE]
    ph  <- ph0[match(colnames(mat), ph0$Sample), , drop = FALSE]
    if (anyNA(ph$Sample) || any(!ph$Sample %in% train_ids)) {
      stop("Test sample or missing phenotype in training Cox: ", scope)
    }
    stopifnot(identical(colnames(mat), ph$Sample)) 
    cat("[INFO]", dt, scope, "matrix dims:", dim(mat), "\n")

    res_list <- mclapply(seq_len(nrow(mat)), function(i){
      y <- suppressWarnings(as.numeric(mat[i,])); sd_x <- sd(y, na.rm=TRUE)
      if (!is.finite(sd_x) || all(is.na(y))) return(NULL)

      dfm <- data.frame(
        methyl=y,
        time=ph$time, status=ph$status,
        Age=ph$Age, Stage_fac=ph$Stage_fac,
        Immune_total=ph$Immune_total, Purity=ph$Purity
      )
      if (include_sex && "Sex_fac" %in% names(ph)) dfm$Sex_fac <- ph$Sex_fac

      rhs <- c("Immune_total","Purity","Age", if (isTRUE(use_stage)) "Stage_fac",
               if (include_sex && "Sex_fac" %in% names(dfm)) "Sex_fac")
      rhs <- rhs[rhs %in% names(dfm)]
      fmla <- as.formula(paste("Surv(time,status) ~", paste(c("methyl", rhs), collapse=" + ")))

      need_cols <- unique(c("methyl","time","status", rhs))
      ok <- complete.cases(dfm[, need_cols])
      df_ok <- dfm[ok, , drop=FALSE]
      if (nrow(df_ok) < 5) return(NULL)

      fit <- tryCatch(coxph(fmla, data=df_ok, ties="efron"), error=function(e) NULL)
      if (is.null(fit)) return(NULL)

      cm <- tryCatch(summary(fit)$coefficients, error=function(e) NULL)
      if (is.null(cm) || (!"methyl" %in% rownames(cm))) return(NULL)

      beta_r <- suppressWarnings(cm["methyl","coef"])
      se_r   <- suppressWarnings(cm["methyl","se(coef)"])
      pval   <- suppressWarnings(cm["methyl","Pr(>|z|)"])
      if (!is.finite(beta_r)) return(NULL)

      hr     <- exp(beta_r); ci_l <- exp(beta_r-1.96*se_r); ci_h <- exp(beta_r+1.96*se_r)
      beta_s <- beta_r*sd_x; se_s <- se_r*sd_x; stdhr <- exp(beta_s)
      cil_s  <- exp(beta_s-1.96*se_s); cih_s <- exp(beta_s+1.96*se_s)

      list(
        cpg=rownames(mat)[i],
        Estimate=beta_r, StdErr=se_r, HR=hr, CI_lower=ci_l, CI_upper=ci_h,
        StdEffectSize=beta_s, SE_StdHR=se_s, StdHR=stdhr,
        CI_lower_std=cil_s, CI_upper_std=cih_s,
        SDx=sd_x, N=sum(ok), pValue=if (is.finite(pval)) pval else NA_real_
      )
    }, mc.cores=opts$cores)

    res_list <- Filter(Negate(is.null), res_list)
    dt_res <- if (length(res_list)>0) rbindlist(res_list) else
      data.table(cpg=character(), Estimate=numeric(), StdErr=numeric(), HR=numeric(),
                 CI_lower=numeric(), CI_upper=numeric(),
                 StdEffectSize=numeric(), SE_StdHR=numeric(), StdHR=numeric(),
                 CI_lower_std=numeric(), CI_upper_std=numeric(),
                 SDx=numeric(), N=integer(), pValue=numeric())
    if (nrow(dt_res)>0) {
      suppressWarnings(dt_res[, FDR := p.adjust(pValue, method="fdr")])
    } else {
      dt_res[, FDR := numeric()]
    }
    setcolorder(dt_res, c("cpg","Estimate","StdErr","HR","CI_lower","CI_upper",
                          "StdEffectSize","SE_StdHR","StdHR","CI_lower_std","CI_upper_std",
                          "SDx","N","pValue","FDR"))
    dt_res <- dt_res[order(pValue)]

    out_csv <- file.path(base_out_dir, "cox", proj, scope,
                         sprintf("%s_%s_cox_purity_immune.csv", proj, dt))
    dir.create(dirname(out_csv), recursive=TRUE, showWarnings=FALSE)
    fwrite(dt_res, file=out_csv)
    cat("[DONE]", dt, "Cox CSV saved to:", out_csv, "\n")

    }
  }

  cat("[CV DONE]", proj, "15 folds x 3 methylation measures\n")

}

cat("\n[ALL DONE] Outputs under:", base_out_dir, "\n")
