#!/usr/bin/env Rscript
# Prepare binary brain-meta candidates, London blood M-values and EpiDISH covariates for nested MRS modelling.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--binary_meta_dir", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--pair_file", type = "character"),
  make_option("--blood_cell_file", type = "character"),
  make_option("--beta_file", type = "character"),
  make_option("--chalm_file", type = "character"),
  make_option("--camda_file", type = "character"),
  make_option("--candidate_p", type = "double", default = 0.05),
  make_option("--late_cutoff", type = "double", default = 5),
  make_option("--epsilon", type = "double", default = 1e-6),
  make_option("--methods", type = "character", default = "beta,chalm,camda")
)

opts <- parse_args(OptionParser(option_list = option_list))
required <- c("binary_meta_dir", "output_dir", "pair_file", "blood_cell_file",
              "beta_file", "chalm_file", "camda_file")
missing <- required[vapply(required, function(x) {
  is.null(opts[[x]]) || !nzchar(opts[[x]])
}, logical(1))]
if (length(missing)) stop("Missing required options: ", paste(missing, collapse = ", "))
if (!is.finite(opts$candidate_p) || opts$candidate_p <= 0 || opts$candidate_p >= 1 ||
    !is.finite(opts$late_cutoff) || opts$epsilon <= 0 || opts$epsilon >= 0.5) {
  stop("Invalid candidate P-value, Braak cutoff or epsilon.")
}

DIRS <- list(
  root = opts$output_dir,
  source = file.path(opts$output_dir, "source_data"),
  audit = file.path(opts$output_dir, "audit")
)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

METHODS <- unique(tolower(trimws(strsplit(opts$methods, ",", fixed = TRUE)[[1L]])))
METHODS <- METHODS[nzchar(METHODS)]
if (!length(METHODS) || any(!METHODS %in% c("beta", "chalm", "camda"))) {
  stop("--methods must contain only beta, chalm and/or camda.")
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
normalize_id <- function(x) gsub("\\.", "-", trimws(as.character(x)))
first_col <- function(x, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit)) return(hit[[1L]])
  if (required) stop("Missing column; tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

find_meta_file <- function(method) {
  files <- list.files(
    opts$binary_meta_dir,
    pattern = paste0("_", method, "_raw_filtered\\.csv$"),
    full.names = TRUE, ignore.case = TRUE
  )
  if (length(files) != 1L) {
    stop("Expected one binary brain-meta file for ", method, "; found ", length(files), ".")
  }
  files[[1L]]
}

read_meta <- function(method) {
  path <- find_meta_file(method)
  x <- fread(path, data.table = TRUE, check.names = FALSE)
  cpg_col <- first_col(x, c("cpg", "CpG", "Name"))
  effect_col <- first_col(x, c("meta_Estimate", "meta_EffectSize", "Estimate_fixed"))
  p_col <- first_col(x, c("pval_final", "pval_fixed", "pval_random"))
  fdr_col <- first_col(x, c("fdr_all", "fdr", "fdr_U"), required = FALSE)
  fdr <- if (is.na(fdr_col)) rep(NA_real_, nrow(x)) else to_num(x[[fdr_col]])
  out <- data.table(
    endpoint = "binary",
    cpg = trimws(as.character(x[[cpg_col]])),
    value_method = method,
    brain_meta_weight = to_num(x[[effect_col]]),
    brain_meta_p = to_num(x[[p_col]]),
    brain_meta_fdr = fdr,
    brain_meta_file = path
  )
  out <- out[
    !is.na(cpg) & nzchar(cpg) & is.finite(brain_meta_weight) &
      is.finite(brain_meta_p) & brain_meta_p >= 0 & brain_meta_p < opts$candidate_p
  ]
  out[, feature_id := paste(cpg, value_method, sep = "__")]
  unique(out, by = c("endpoint", "feature_id"))
}

features <- rbindlist(lapply(METHODS, read_meta), use.names = TRUE, fill = TRUE)
if (!nrow(features)) stop("No binary brain-meta locus passed P < ", opts$candidate_p, ".")
features[, p_rank := frank(brain_meta_p, ties.method = "first")]
fwrite(features, file.path(DIRS$audit, "binary_brain_meta_p005_candidates.csv"))

read_blood_samples <- function() {
  x <- fread(opts$pair_file, data.table = TRUE, check.names = FALSE)
  sample_col <- first_col(x, c("blood_sample", "blood_sample_id", "Sample", "sample"))
  stage_col <- first_col(x, c("blood_stage", "braak_stage", "stage"))
  age_col <- first_col(x, c("blood_age", "blood_age.blood", "age"))
  sex_col <- first_col(x, c("blood_sex", "blood_sex.blood", "sex"))
  slide_col <- first_col(x, c("blood_slide", "blood_slide.blood", "slide", "Sentrix_ID", "array_id"))
  value_or_na <- function(column) if (is.na(column)) rep(NA_character_, nrow(x)) else x[[column]]
  out <- data.table(
    blood_sample = normalize_id(x[[sample_col]]),
    stage = to_num(x[[stage_col]]),
    blood_age = to_num(value_or_na(age_col)),
    blood_sex = as.character(value_or_na(sex_col)),
    blood_slide = as.character(value_or_na(slide_col))
  )
  out <- out[!is.na(blood_sample) & nzchar(blood_sample) & is.finite(stage)]
  if (anyDuplicated(out$blood_sample)) stop("Duplicate blood sample IDs in the phenotype file.")
  out[, late_braak := as.integer(stage >= opts$late_cutoff)]
  out
}

merge_blood_cells <- function(samples) {
  x <- fread(opts$blood_cell_file, data.table = TRUE, check.names = FALSE)
  id_col <- first_col(x, c("GSM_ID", "blood_sample", "sample", "Sample"))
  x[, join_sample := normalize_id(get(id_col))]
  value_cols <- setdiff(names(x), c(id_col, "join_sample"))
  numeric_cols <- value_cols[vapply(x[, ..value_cols], function(z) {
    any(is.finite(to_num(z)))
  }, logical(1))]
  for (column in numeric_cols) set(x, j = column, value = to_num(x[[column]]))
  names_out <- ifelse(
    startsWith(numeric_cols, "blood_cell_"), numeric_cols,
    paste0("blood_cell_", make.names(numeric_cols))
  )
  setnames(x, numeric_cols, names_out)
  x <- unique(x[, c("join_sample", names_out), with = FALSE], by = "join_sample")
  merge(samples, x, by.x = "blood_sample", by.y = "join_sample", all.x = TRUE, sort = FALSE)
}

collapse_blood_cells <- function(x) {
  cell_columns <- grep("^blood_cell_", names(x), value = TRUE)
  cell_keys <- setNames(
    gsub("[^a-z0-9]", "", tolower(sub("^blood_cell_", "", cell_columns))),
    cell_columns
  )
  definitions <- list(
    blood_major_B = c("b", "bcell", "bcells"),
    blood_major_NK = c("nk", "nkcell", "nkcells"),
    blood_major_CD4T = c("cd4", "cd4t", "cd4cell", "cd4cells", "cd4tcell", "cd4tcells"),
    blood_major_CD8T = c("cd8", "cd8t", "cd8cell", "cd8cells", "cd8tcell", "cd8tcells"),
    blood_major_Mono = c("mono", "monocyte", "monocytes"),
    blood_major_Neutro = c("neu", "neutro", "neutrophil", "neutrophils"),
    blood_major_Eos = c("eos", "eosinophil", "eosinophils")
  )
  mapping <- lapply(names(definitions), function(target) {
    sources <- names(cell_keys)[cell_keys %in% definitions[[target]]]
    if (length(sources) != 1L) {
      stop("EpiDISH cell file must have exactly one estimate for ", target, ".")
    }
    values <- as.matrix(x[, ..sources])
    storage.mode(values) <- "numeric"
    finite_count <- rowSums(is.finite(values))
    total <- rowSums(values, na.rm = TRUE)
    total[finite_count == 0L] <- NA_real_
    set(x, j = target, value = total)
    data.table(cell_type = target, source_columns = paste(sources, collapse = ";"))
  })
  list(data = x, mapping = rbindlist(mapping))
}

samples <- merge_blood_cells(read_blood_samples())
cells <- collapse_blood_cells(samples)
samples <- unique(cells$data, by = "blood_sample")
setorder(samples, blood_sample)
required_complete <- c("blood_age", "blood_sex", "blood_slide", "late_braak",
                       "blood_major_B", "blood_major_NK", "blood_major_CD4T",
                       "blood_major_CD8T", "blood_major_Mono", "blood_major_Neutro",
                       "blood_major_Eos")
if (any(!complete.cases(samples[, ..required_complete]))) {
  stop("London blood samples have missing age, sex, slide, outcome or EpiDISH estimates.")
}
if (nrow(samples) != 71L || sum(samples$late_braak == 1L) != 49L ||
    sum(samples$late_braak == 0L) != 22L) {
  stop("Expected 71 blood samples: 49 Braak V-VI and 22 Braak 0-IV.")
}
fwrite(cells$mapping, file.path(DIRS$audit, "epidish_blood_cell_mapping.csv"))
fwrite(
  data.table(
    n_samples = nrow(samples),
    n_advanced = sum(samples$late_braak == 1L),
    n_nonadvanced = sum(samples$late_braak == 0L)
  ),
  file.path(DIRS$audit, "binary_blood_sample_counts.csv")
)

matrix_paths <- list(beta = opts$beta_file, chalm = opts$chalm_file, camda = opts$camda_file)
read_matrix <- function(path, method, sample_ids, cpgs) {
  if (!file.exists(path)) stop("Matrix not found: ", path)
  x <- fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  ids <- trimws(as.character(x[[1L]]))
  selected_rows <- which(ids %in% unique(cpgs))
  if (!length(selected_rows)) {
    return(list(
      matrix = matrix(numeric(), nrow = 0L, ncol = length(sample_ids),
                      dimnames = list(character(), sample_ids)),
      scale = "no_candidates", fraction_in_0_1 = NA_real_
    ))
  }
  mat <- as.matrix(x[selected_rows, -1L, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- ids[selected_rows]
  colnames(mat) <- normalize_id(colnames(mat))
  if (anyDuplicated(rownames(mat)) || anyDuplicated(colnames(mat))) {
    stop("Duplicate CpG or sample IDs in ", method, " matrix.")
  }
  missing_samples <- setdiff(sample_ids, colnames(mat))
  if (length(missing_samples)) stop(method, " matrix lacks one or more blood samples.")
  finite <- mat[is.finite(mat)]
  fraction <- mean(finite >= 0 & finite <= 1)
  if (!is.finite(fraction)) stop(method, " matrix contains no finite values.")
  if (fraction >= 0.95) {
    mat[is.finite(mat)] <- pmin(pmax(mat[is.finite(mat)], opts$epsilon), 1 - opts$epsilon)
    mat <- log2(mat / (1 - mat))
    scale <- "proportion_to_M"
  } else {
    scale <- "retained_M"
  }
  list(matrix = mat[, sample_ids, drop = FALSE], scale = scale, fraction_in_0_1 = fraction)
}

candidate_counts <- features[, .(n_features = .N, n_unique_cpg = uniqueN(cpg)), by = "value_method"]
fwrite(candidate_counts, file.path(DIRS$audit, "binary_brain_meta_candidate_counts.csv"))
common_samples <- samples$blood_sample
method_matrices <- list()
scale_audit <- list()
for (method in METHODS) {
  cpgs <- features[value_method == method, cpg]
  result <- read_matrix(matrix_paths[[method]], method, common_samples, cpgs)
  method_matrices[[method]] <- result$matrix
  scale_audit[[method]] <- data.table(
    method = method, scale_decision = result$scale,
    fraction_in_0_1 = result$fraction_in_0_1,
    n_features = nrow(result$matrix), n_samples = ncol(result$matrix)
  )
}

available_features <- unlist(lapply(names(method_matrices), function(method) {
  paste(rownames(method_matrices[[method]]), method, sep = "__")
}), use.names = FALSE)
features <- features[feature_id %in% available_features]
if (!nrow(features)) stop("No P<0.05 brain-meta candidate was present in the blood matrices.")
fwrite(features, file.path(DIRS$source, "binary_brain_meta_p005_mrs_candidates.csv"))
fwrite(rbindlist(scale_audit), file.path(DIRS$audit, "blood_matrix_scale_audit.csv"))

saveRDS(
  list(
    feature_metadata = features,
    sample_data = samples,
    method_matrices = method_matrices,
    late_braak_cutoff = opts$late_cutoff
  ),
  file.path(DIRS$source, "brain_meta_blood_mrs_input.rds"),
  compress = FALSE
)
message("[DONE] Binary brain-meta candidates and London blood M-values prepared: ",
        file.path(DIRS$source, "brain_meta_blood_mrs_input.rds"))
