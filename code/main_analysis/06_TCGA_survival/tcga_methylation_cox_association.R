# Run TCGA methylation Cox scans.

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
  make_option("--project-map", type = "character"),
  make_option("--phenotype-dir", type = "character"),
  make_option("--beta-dir", type = "character"),
  make_option("--chalm-dir", type = "character"),
  make_option("--camda-dir", type = "character"),
  make_option("--reference-dir", type = "character"),
  make_option("--output-dir", type = "character"),
  make_option("--projects", type = "character", default = ""),
  make_option("--drop-sex-projects", type = "character", default = "TCGA-BRCA,TCGA-PRAD"),
  make_option("--cores", type = "integer", default = 1L),
  make_option("--na-row-threshold", type = "double", default = 0.7),
  make_option("--na-column-threshold", type = "double", default = 0.7)
)
opts <- parse_args(OptionParser(option_list = option_list))
required <- c("project_map", "phenotype_dir", "beta_dir", "chalm_dir", "camda_dir", "reference_dir", "output_dir")
if (any(vapply(required, function(name) is.null(opts[[name]]) || !nzchar(opts[[name]]), logical(1))) || opts$cores < 1L) {
  stop("Project map, input roots, output directory, and positive core count are required.")
}

project_map <- fread(opts$project_map, data.table = TRUE, check.names = FALSE)
if (!all(c("project", "tissue") %in% names(project_map))) stop("Project map requires project and tissue columns.")
project_map[, project := trimws(as.character(project))]
project_map[, tissue := tolower(trimws(as.character(tissue)))]
if (anyDuplicated(project_map$project) || anyNA(project_map$project) || anyNA(project_map$tissue) ||
    any(!nzchar(project_map$project)) || any(!nzchar(project_map$tissue))) {
  stop("Project map has empty or duplicate entries.")
}
projects <- trimws(strsplit(opts$projects, ",", fixed = TRUE)[[1L]])
projects <- projects[nzchar(projects)]
if (!length(projects)) projects <- project_map$project
if (length(setdiff(projects, project_map$project))) stop("A requested project is absent from the project map.")
drop_sex_projects <- unique(c(
  "TCGA-BRCA", "TCGA-PRAD",
  trimws(strsplit(opts$drop_sex_projects, ",", fixed = TRUE)[[1L]])
))
drop_sex_projects <- drop_sex_projects[nzchar(drop_sex_projects)]

ref_env <- new.env(parent = globalenv())
for (object in c("centDHSbloodDMC.m", "cent12CT.m", "cent12CT450k.m", "centUniLIFE.m",
                 "centEpiFibIC.m", "centEpiFibFatIC.m", "centBloodSub.m")) {
  suppressWarnings(data(list = object, package = "EpiDISH", envir = ref_env))
}
reference_files <- list.files(opts$reference_dir, pattern = "Ref\\.rda$", full.names = TRUE)
if (!length(reference_files)) stop("No tissue reference RData files were found.")
for (path in reference_files) load(path, envir = ref_env)

tissue_references <- list(
  solid_generic = list(object = "centEpiFibIC.m", method = "RPC", hierarchical = FALSE),
  brain = list(object = "mrefBrain.m", method = "wRPC", hierarchical = FALSE),
  breast = list(object = "mrefBreast.m", method = "wRPC", hierarchical = FALSE),
  prostate = list(object = "mrefProstate.m", method = "wRPC", hierarchical = FALSE),
  skin = list(object = "mrefSkin.m", method = "wRPC", hierarchical = FALSE),
  lung = list(object = "mrefLung.m", method = "wRPC", hierarchical = FALSE),
  bladder = list(object = "mrefBladder.m", method = "wRPC", hierarchical = FALSE),
  kidney = list(object = "Kidney_Mref.m", method = "wRPC", hierarchical = FALSE),
  liver = list(object = "mrefLiver.m", method = "wRPC", hierarchical = FALSE),
  esophagus = list(object = "Eso_Mref.m", method = "wRPC", hierarchical = FALSE),
  pancreas = list(object = "mrefPancreas.m", method = "wRPC", hierarchical = FALSE),
  colon = list(object = "Colon_Mref.m", method = "wRPC", hierarchical = FALSE)
)

normalize_id <- function(x, keep_first_16 = FALSE) {
  value <- toupper(gsub("\\.", "-", trimws(as.character(x))))
  if (keep_first_16) value <- sub("^([A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]+-[A-Z0-9]{3}).*$", "\\1", value)
  value
}

tcga_stage <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value <- gsub("[_ -]+", " ", value)
  value <- sub("^stage ", "", value)
  value <- sub("^stage:", "", value)
  value[value %in% c("not reported", "unknown", "x", "tx", "nx", "na", "n/a", "")] <- NA_character_
  main <- sub("^(i{1,3}|iv|[1-4]).*$", "\\1", value)
  stage <- ifelse(main %in% c("i", "1"), "S1",
           ifelse(main %in% c("ii", "2"), "S2",
           ifelse(main %in% c("iii", "3"), "S3",
           ifelse(main %in% c("iv", "4"), "S4", "Unknown"))))
  factor(stage, levels = c("S1", "S2", "S3", "S4", "Unknown"))
}

sample_type_group <- function(x) {
  value <- tolower(trimws(as.character(x)))
  group <- ifelse(value == "primary tumor", "Tumor",
           ifelse(value == "solid tissue normal", "Normal",
           ifelse(grepl("normal", value), "Normal",
           ifelse(grepl("tumou?r|metastatic|cancer", value), "Other Tumor", "Other"))))
  factor(group, levels = c("Tumor", "Normal", "Other Tumor", "Other"))
}

read_phenotype <- function(project) {
  path <- file.path(opts$phenotype_dir, paste0(project, ".GDC_phenotype.tsv"))
  if (!file.exists(path)) stop("Missing GDC phenotype table: ", path)
  frame <- fread(path, sep = "\t", data.table = FALSE, check.names = FALSE)
  get_column <- function(name) if (name %in% names(frame)) frame[[name]] else rep(NA, nrow(frame))
  sample <- normalize_id(get_column("submitter_id.samples"))
  days_death <- suppressWarnings(as.numeric(get_column("days_to_death.demographic")))
  days_followup <- suppressWarnings(as.numeric(get_column("days_to_last_follow_up.diagnoses")))
  vital <- tolower(trimws(as.character(get_column("vital_status.demographic"))))
  time <- ifelse(is.finite(days_death), days_death, days_followup) / 365.25
  status <- ifelse(vital == "dead", 1L, ifelse(vital == "alive", 0L, NA_integer_))
  sex <- tolower(trimws(as.character(get_column("gender.demographic"))))
  sex <- factor(ifelse(sex %in% c("male", "female"), tools::toTitleCase(sex), "Unknown"),
                levels = c("Male", "Female", "Unknown"))
  out <- data.frame(
    Sample = sample,
    Age = suppressWarnings(as.numeric(get_column("age_at_index.demographic"))),
    Sex_fac = sex,
    Stage_fac = tcga_stage(get_column("tumor_stage.diagnoses")),
    time = time,
    status = status,
    SampleType_group = sample_type_group(get_column("sample_type.samples")),
    stringsAsFactors = FALSE
  )
  stage_value <- as.character(out$Stage_fac)
  sex_value <- as.character(out$Sex_fac)
  keep <- is.finite(out$time) & out$time > 0 & !is.na(out$status) & is.finite(out$Age) &
    sex_value %in% c("Male", "Female") & stage_value %in% paste0("S", 1:4)
  out <- out[keep, , drop = FALSE]
  out[!duplicated(out$Sample), , drop = FALSE]
}

read_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  frame <- fread(path, sep = "\t", data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  if (ncol(frame) < 2L) stop("Methylation matrix has fewer than two columns: ", path)
  feature <- trimws(as.character(frame[[1L]]))
  matrix <- as.matrix(data.frame(lapply(frame[, -1L, drop = FALSE], function(x) suppressWarnings(as.numeric(x))),
                                 check.names = FALSE))
  rownames(matrix) <- feature
  colnames(matrix) <- normalize_id(colnames(frame)[-1L])
  matrix <- matrix[!duplicated(rownames(matrix)), !duplicated(colnames(matrix)), drop = FALSE]
  matrix[!is.finite(matrix)] <- NA_real_
  matrix[matrix < 0] <- 0
  matrix[matrix > 1] <- 1
  row_fraction <- rowSums(is.na(matrix)) / ncol(matrix)
  matrix <- matrix[row_fraction < opts$na_row_threshold, , drop = FALSE]
  if (!nrow(matrix)) return(NULL)
  column_fraction <- colSums(is.na(matrix)) / nrow(matrix)
  matrix <- matrix[, column_fraction < opts$na_column_threshold, drop = FALSE]
  matrix
}

load_matrix <- function(project, method) {
  if (method == "beta") {
    path <- file.path(opts$beta_dir, paste0(project, "_assay.tsv"))
  } else if (method == "chalm") {
    path <- file.path(opts$chalm_dir, paste0(project, ".chalm.filtered.tsv"))
  } else {
    path <- file.path(opts$camda_dir, paste0(project, ".camda.filtered.tsv"))
  }
  read_matrix(path)
}

reference_for <- function(tissue) {
  spec <- tissue_references[[tissue]]
  if (is.null(spec)) stop("Unsupported tissue reference: ", tissue)
  if (!exists(spec$object, envir = ref_env, inherits = FALSE)) stop("Missing reference object: ", spec$object)
  get(spec$object, envir = ref_env, inherits = FALSE)
}

estimate_cells <- function(beta, tissue) {
  spec <- tissue_references[[tissue]]
  reference <- reference_for(tissue)
  missing <- !is.finite(beta)
  if (any(missing)) {
    medians <- apply(beta, 1L, median, na.rm = TRUE)
    medians[!is.finite(medians)] <- 0.5
    for (i in which(rowSums(missing) > 0L)) beta[i, missing[i, ]] <- medians[i]
  }
  if (spec$method == "wRPC") {
    promoted <- constAvBetaTSS(beta.m = beta, type = "450k")
    common <- intersect(rownames(promoted), rownames(reference))
    if (!length(common)) stop("No shared probes with the requested tissue reference.")
    estimates <- wRPC(promoted[common, , drop = FALSE], reference[common, , drop = FALSE])$estF
  } else {
    estimates <- epidish(beta.m = beta, ref.m = reference, method = spec$method)$estF
  }
  as.data.frame(estimates)
}

derive_tme <- function(cell) {
  immune_names <- c("Macro", "Lym", "Gran", "Mono", "Plasma", "Tcell", "Mye")
  endothelial_names <- "Endo"
  stromal_names <- c("Stellate", "CM", "SMC", "Fat", "Stromal")
  columns <- setdiff(names(cell), "Sample")
  immune <- intersect(immune_names, columns)
  endothelial <- intersect(endothelial_names, columns)
  stromal <- intersect(stromal_names, columns)
  immune_total <- if (length(immune)) rowSums(cell[, immune, drop = FALSE], na.rm = TRUE) else rep(0, nrow(cell))
  endothelial_total <- if (length(endothelial)) rowSums(cell[, endothelial, drop = FALSE], na.rm = TRUE) else rep(0, nrow(cell))
  stromal_total <- if (length(stromal)) rowSums(cell[, stromal, drop = FALSE], na.rm = TRUE) else rep(0, nrow(cell))
  data.frame(Sample = cell$Sample, Immune_total = immune_total,
             Purity = pmin(pmax(1 - immune_total - endothelial_total - stromal_total, 0), 1))
}

fit_one <- function(feature, methyl, pheno, rhs) {
  methyl_value <- suppressWarnings(as.numeric(methyl))
  if (!is.finite(sd(methyl_value, na.rm = TRUE)) || sd(methyl_value, na.rm = TRUE) < 1e-8) return(NULL)
  frame <- data.frame(methyl = methyl_value, time = pheno$time, status = pheno$status,
                      Age = pheno$Age, Stage_fac = pheno$Stage_fac,
                      Immune_total = pheno$Immune_total, Purity = pheno$Purity)
  if ("Sex_fac" %in% names(pheno) && "Sex_fac" %in% rhs) frame$Sex_fac <- pheno$Sex_fac
  used_rhs <- intersect(rhs, names(frame))
  keep <- complete.cases(frame[, unique(c("methyl", "time", "status", used_rhs)), drop = FALSE])
  frame <- droplevels(frame[keep, , drop = FALSE])
  if (nrow(frame) < 5L) return(NULL)
  formula <- as.formula(paste("Surv(time, status) ~", paste(c("methyl", used_rhs), collapse = " + ")))
  fit <- tryCatch(coxph(formula, data = frame, ties = "efron"), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  summary_fit <- summary(fit)$coefficients
  if (!"methyl" %in% rownames(summary_fit)) return(NULL)
  estimate <- unname(summary_fit["methyl", "coef"])
  se <- unname(summary_fit["methyl", "se(coef)"])
  pvalue <- unname(summary_fit["methyl", "Pr(>|z|)"])
  if (!is.finite(estimate) || !is.finite(se) || se <= 0) return(NULL)
  sd_x <- sd(methyl_value, na.rm = TRUE)
  data.table(
    cpg = feature, Estimate = estimate, StdErr = se,
    HR = exp(estimate), CI_lower = exp(estimate - 1.96 * se), CI_upper = exp(estimate + 1.96 * se),
    StdEffectSize = estimate * sd_x, SE_StdHR = se * sd_x,
    StdHR = exp(estimate * sd_x), CI_lower_std = exp(estimate * sd_x - 1.96 * se * sd_x),
    CI_upper_std = exp(estimate * sd_x + 1.96 * se * sd_x),
    SDx = sd_x, N = sum(keep), pValue = pvalue
  )
}

scan_method <- function(project, method, pheno, include_sex) {
  matrix <- load_matrix(project, method)
  if (is.null(matrix)) return(NULL)
  sample_ids <- intersect(colnames(matrix), pheno$Sample)
  if (length(sample_ids) < 5L) return(NULL)
  matrix <- matrix[, sample_ids, drop = FALSE]
  ph <- pheno[match(sample_ids, pheno$Sample), , drop = FALSE]
  stage_levels <- unique(as.character(pheno$Stage_fac))
  use_stage <- length(stage_levels) > 1L
  rhs <- c("Immune_total", "Purity", "Age", if (use_stage) "Stage_fac", if (include_sex) "Sex_fac")
  rhs <- rhs[rhs %in% names(ph)]
  results <- mclapply(seq_len(nrow(matrix)), function(i) fit_one(rownames(matrix)[i], matrix[i, ], ph, rhs),
                      mc.cores = opts$cores)
  result <- rbindlist(Filter(Negate(is.null), results), fill = TRUE)
  if (!nrow(result)) return(result)
  setorder(result, pValue)
  result
}

for (project in projects) {
  results <- list()
  project_tissue <- project_map$tissue[match(project, project_map$project)]
  phenotype <- read_phenotype(project)
  tumor <- phenotype[SampleType_group == "Tumor"]
  if (nrow(tumor) < 10L) next
  beta <- load_matrix(project, "beta")
  if (is.null(beta)) next
  beta_samples <- intersect(colnames(beta), tumor$Sample)
  if (length(beta_samples) < 5L) next
  beta <- beta[, beta_samples, drop = FALSE]
  tumor <- tumor[match(beta_samples, Sample)]
  cell <- estimate_cells(beta, project_tissue)
  cell$Sample <- rownames(cell)
  tme <- derive_tme(cell)
  pheno <- merge(tumor, tme, by = "Sample", sort = FALSE)
  include_sex <- !project %in% drop_sex_projects
  for (method in c("beta", "chalm", "camda")) {
    result <- scan_method(project, method, pheno, include_sex)
    if (is.null(result)) next
    if (is.null(results[[project]])) results[[project]] <- list()
    results[[project]][[method]] <- result
  }

  if (is.null(results[[project]]) ||
      !all(c("beta", "chalm", "camda") %in% names(results[[project]])) ||
      any(vapply(results[[project]][c("beta", "chalm", "camda")], nrow, integer(1)) == 0L)) {
    stop("All three methylation representations are required for ", project)
  }

  if (!is.null(results[[project]]$beta) && nrow(results[[project]]$beta)) {
    results[[project]]$beta[, `:=`(
      FDR = p.adjust(pValue, method = "BH"),
      FDR_family = "EWAS"
    )]
  }

  ha_tables <- results[[project]][c("chalm", "camda")]
  if (!all(vapply(ha_tables, function(x) !is.null(x) && nrow(x) > 0L, logical(1)))) {
    stop("Both CHALM and CAMDA Cox results are required for joint haEWAS FDR: ", project)
  }
  ha_pooled <- rbindlist(ha_tables, idcol = "method")
  ha_pooled[, `:=`(
    FDR = p.adjust(pValue, method = "BH"),
    FDR_family = "haEWAS_joint"
  )]
  for (method in c("chalm", "camda")) {
    method_result <- ha_pooled[ha_pooled$method == method]
    method_result[, method := NULL]
    results[[project]][[method]] <- method_result
  }

  project_dir <- file.path(opts$output_dir, project)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  for (method in c("beta", "chalm", "camda")) {
    result <- results[[project]][[method]]
    if (is.null(result) || !nrow(result)) next
    output <- file.path(project_dir, paste0(project, "_", method, "_cox_purity_immune.csv"))
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    fwrite(result, output)
  }

  beta_significant <- results[[project]]$beta[FDR < 0.05]
  fwrite(beta_significant, file.path(project_dir, paste0(project, "_ewas_significant_loci.csv")))

  ha_significant <- rbindlist(lapply(c("chalm", "camda"), function(method) {
    result <- results[[project]][[method]][FDR < 0.05, .(cpg, pValue, FDR)]
    result[, method := method]
    result
  }), use.names = TRUE, fill = TRUE)
  ha_union <- ha_significant[, .(
    SignificantMethods = paste(toupper(sort(unique(method))), collapse = ";"),
    CHALM_P = if (any(method == "chalm")) pValue[match("chalm", method)] else NA_real_,
    CHALM_FDR = if (any(method == "chalm")) FDR[match("chalm", method)] else NA_real_,
    CAMDA_P = if (any(method == "camda")) pValue[match("camda", method)] else NA_real_,
    CAMDA_FDR = if (any(method == "camda")) FDR[match("camda", method)] else NA_real_
  ), by = "cpg"]
  fwrite(ha_union, file.path(project_dir, paste0(project, "_haewas_significant_loci.csv")))
}
