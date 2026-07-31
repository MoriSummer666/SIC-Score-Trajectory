#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

settings <- list(
  cdic_input = file.path("data", "CDIC_D1D3_original_modeling_patients_with_trajectory.xlsx"),
  mimic_input = file.path("data", "MIMIC_D1D3_original_modeling_patients_with_trajectory.xlsx"),
  input_sheet = "Sheet1",
  out_dir = file.path("results", "03_random_forest_imputation"),
  seed = 20260524L,
  ntree = 100L,
  maxiter = 10L,
  variablewise = TRUE,
  backend = "ranger",
  parallelize = "no",
  num_threads = max(1L, parallel::detectCores(logical = TRUE) - 1L, na.rm = TRUE),
  missing_threshold = 0.40,
  winsor_probs = c(0.01, 0.99),
  missing_tokens = c("", "NA", "N/A", ".", "NULL", "NaN"),
  imputation_vars = c(
    "death_28d",
    "survival_time_28d_days",
    "age",
    "sex_male",
    "bmi_kg_m2",
    "temp_max_d1",
    "rr_max_d1",
    "hr_max_d1",
    "arterial_ph_d1",
    "lactate_d1",
    "bun_d1",
    "alt_d1",
    "ast_d1",
    "alp_d1",
    "glucose_d1",
    "wbc_d1",
    "hemoglobin_d1",
    "observed_cci",
    "heparin_group_preD3",
    "trajectory"
  ),
  predictor_only_vars = c(
    "death_28d",
    "survival_time_28d_days",
    "age",
    "sex_male",
    "observed_cci",
    "heparin_group_preD3",
    "trajectory"
  ),
  winsor_vars = c(
    "bmi_kg_m2",
    "temp_max_d1",
    "rr_max_d1",
    "hr_max_d1",
    "arterial_ph_d1",
    "lactate_d1",
    "bun_d1",
    "alt_d1",
    "ast_d1",
    "alp_d1",
    "glucose_d1",
    "wbc_d1",
    "hemoglobin_d1"
  ),
  binary_vars = c("death_28d", "sex_male"),
  nominal_levels = list(
    heparin_group_preD3 = c("No_by_D3", "LE24h", "H24_72h"),
    trajectory = c("T1", "T2", "T3", "T4", "T5")
  )
)

settings$target_vars <- setdiff(settings$imputation_vars, settings$predictor_only_vars)

args <- commandArgs(trailingOnly = TRUE)
aliases <- c(
  "cdic-input" = "cdic_input",
  "mimic-input" = "mimic_input",
  "input-sheet" = "input_sheet",
  "out-dir" = "out_dir",
  "threads" = "num_threads"
)
for (arg in args) {
  if (!startsWith(arg, "--")) next
  kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  key <- kv[1]
  if (key %in% names(aliases)) key <- aliases[[key]]
  key <- gsub("-", "_", key)
  value <- if (length(kv) > 1L) paste(kv[-1], collapse = "=") else "TRUE"
  if (!key %in% names(settings)) next
  old <- settings[[key]]
  if (is.logical(old) && length(old) == 1L) {
    settings[[key]] <- tolower(value) %in% c("true", "t", "1", "yes", "y")
  } else if (is.integer(old) && length(old) == 1L) {
    settings[[key]] <- as.integer(value)
  } else if (is.numeric(old) && length(old) == 1L) {
    settings[[key]] <- as.numeric(value)
  } else if (is.character(old) && length(old) == 1L) {
    settings[[key]] <- value
  }
}

required_packages <- c("readxl", "missForest", "ranger", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

set.seed(settings$seed)
Sys.setenv(
  OMP_NUM_THREADS = as.character(settings$num_threads),
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
dir.create(settings$out_dir, recursive = TRUE, showWarnings = FALSE)

declared_missing <- function(x) {
  z <- trimws(as.character(x))
  is.na(x) | is.na(z) | z == "" | toupper(z) %in% toupper(settings$missing_tokens)
}

numeric_strict <- function(x, variable) {
  z <- trimws(as.character(x))
  missing <- declared_missing(x)
  value <- suppressWarnings(as.numeric(z))
  invalid <- !missing & is.na(value)
  if (any(invalid)) {
    stop(
      "Unexpected nonnumeric value in ",
      variable,
      ": ",
      paste(head(unique(z[invalid]), 10), collapse = ", ")
    )
  }
  value[missing] <- NA_real_
  value
}

binary_strict <- function(x, variable) {
  value <- numeric_strict(x, variable)
  if (any(!is.na(value) & !value %in% c(0, 1))) stop(variable, " must contain only 0, 1, or missing values.")
  as.integer(value)
}

nominal_strict <- function(x, variable, levels) {
  value <- trimws(as.character(x))
  value[declared_missing(x)] <- NA_character_
  invalid <- !is.na(value) & !value %in% levels
  if (any(invalid)) {
    stop(
      "Unexpected level in ",
      variable,
      ": ",
      paste(head(unique(value[invalid]), 10), collapse = ", ")
    )
  }
  value
}

standardize_variables <- function(data) {
  categorical <- c(settings$binary_vars, names(settings$nominal_levels))
  for (variable in setdiff(settings$imputation_vars, categorical)) {
    data[[variable]] <- numeric_strict(data[[variable]], variable)
  }
  for (variable in settings$binary_vars) {
    data[[variable]] <- binary_strict(data[[variable]], variable)
  }
  for (variable in names(settings$nominal_levels)) {
    data[[variable]] <- nominal_strict(
      data[[variable]],
      variable,
      settings$nominal_levels[[variable]]
    )
  }
  data
}

validate_data <- function(data, id_vars) {
  required <- unique(c(id_vars, settings$imputation_vars))
  if (length(setdiff(required, names(data)))) {
    stop("Missing columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  }
  if (any(is.na(data$row_id_mf)) || anyDuplicated(data$row_id_mf)) {
    stop("row_id_mf must be complete and unique.")
  }
  if (any(is.na(data$patient_id_cleaned)) || anyDuplicated(data$patient_id_cleaned)) {
    stop("patient_id_cleaned must be complete and unique.")
  }
  missing_proportion <- vapply(
    data[, settings$imputation_vars, drop = FALSE],
    function(x) mean(is.na(x)),
    numeric(1)
  )
  if (any(missing_proportion > settings$missing_threshold + 1e-12)) {
    failed <- names(missing_proportion)[missing_proportion > settings$missing_threshold + 1e-12]
    stop(
      "Variables exceeding the 40% missingness threshold: ",
      paste0(failed, "=", sprintf("%.1f%%", 100 * missing_proportion[failed]), collapse = "; ")
    )
  }
  incomplete_predictors <- settings$predictor_only_vars[
    vapply(data[, settings$predictor_only_vars, drop = FALSE], anyNA, logical(1))
  ]
  if (length(incomplete_predictors)) {
    stop("Predictor-only variables contain missing values: ", paste(incomplete_predictors, collapse = ", "))
  }
  unique_counts <- vapply(
    data[, settings$imputation_vars, drop = FALSE],
    function(x) length(unique(stats::na.omit(x))),
    integer(1)
  )
  if (any(unique_counts <= 1L)) {
    stop("All-missing or constant variables: ", paste(names(unique_counts)[unique_counts <= 1L], collapse = ", "))
  }
  invisible(TRUE)
}

winsorize <- function(data) {
  for (variable in settings$winsor_vars) {
    value <- as.numeric(data[[variable]])
    observed <- value[!is.na(value)]
    if (length(observed) < 20L) stop("Too few observed values for winsorization: ", variable)
    limits <- quantile(
      observed,
      probs = settings$winsor_probs,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    if (any(!is.finite(limits)) || limits[1] > limits[2]) stop("Invalid winsorization limits: ", variable)
    value[!is.na(value) & value < limits[1]] <- limits[1]
    value[!is.na(value) & value > limits[2]] <- limits[2]
    data[[variable]] <- value
  }
  data
}

prepare_imputation_matrix <- function(data) {
  xmis <- data[, settings$imputation_vars, drop = FALSE]
  for (variable in settings$binary_vars) {
    xmis[[variable]] <- factor(xmis[[variable]], levels = c(0, 1))
  }
  for (variable in names(settings$nominal_levels)) {
    xmis[[variable]] <- factor(
      xmis[[variable]],
      levels = settings$nominal_levels[[variable]]
    )
  }
  numeric_variables <- setdiff(
    settings$imputation_vars,
    c(settings$binary_vars, names(settings$nominal_levels))
  )
  for (variable in numeric_variables) xmis[[variable]] <- as.numeric(xmis[[variable]])
  xmis
}

run_missforest <- function(xmis) {
  function_arguments <- names(formals(missForest::missForest))
  arguments <- list(
    xmis = xmis,
    maxiter = settings$maxiter,
    ntree = settings$ntree,
    variablewise = settings$variablewise,
    verbose = TRUE
  )
  if ("parallelize" %in% function_arguments) arguments$parallelize <- settings$parallelize
  if ("backend" %in% function_arguments) arguments$backend <- settings$backend
  if ("num.threads" %in% function_arguments) arguments$num.threads <- settings$num_threads
  if ("decreasing" %in% function_arguments) arguments$decreasing <- FALSE
  do.call(missForest::missForest, arguments)
}

restore_storage <- function(ximp) {
  output <- as.data.frame(ximp, check.names = FALSE, stringsAsFactors = FALSE)
  for (variable in settings$binary_vars) {
    output[[variable]] <- as.integer(as.character(output[[variable]]))
  }
  for (variable in names(settings$nominal_levels)) {
    output[[variable]] <- as.character(output[[variable]])
  }
  numeric_variables <- setdiff(
    settings$imputation_vars,
    c(settings$binary_vars, names(settings$nominal_levels))
  )
  for (variable in numeric_variables) output[[variable]] <- as.numeric(output[[variable]])
  output[, settings$imputation_vars, drop = FALSE]
}

equal_vector <- function(x, y, tolerance = 1e-8) {
  if (length(x) != length(y)) return(FALSE)
  if ((is.numeric(x) || is.integer(x)) && (is.numeric(y) || is.integer(y))) {
    same_missing <- is.na(x) & is.na(y)
    both <- !is.na(x) & !is.na(y)
    result <- same_missing
    result[both] <- abs(as.numeric(x[both]) - as.numeric(y[both])) <= tolerance
    return(all(result))
  }
  all((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & as.character(x) == as.character(y)))
}

impute_cohort <- function(cohort, input_file, id_vars, output_file) {
  if (!file.exists(input_file)) stop(cohort, " input file not found: ", input_file)
  raw <- readxl::read_excel(
    input_file,
    sheet = settings$input_sheet,
    na = settings$missing_tokens,
    guess_max = 100000
  ) |>
    as.data.frame(check.names = FALSE, stringsAsFactors = FALSE)
  original_names <- names(raw)
  original_rows <- nrow(raw)
  standardized <- standardize_variables(raw)
  validate_data(standardized, id_vars)
  winsorized <- winsorize(standardized)
  xmis <- prepare_imputation_matrix(winsorized)
  set.seed(settings$seed)
  fit <- run_missforest(xmis)
  imputed <- restore_storage(fit$ximp)
  if (nrow(imputed) != original_rows || !identical(names(imputed), settings$imputation_vars)) {
    stop(cohort, " missForest output structure changed.")
  }
  completed <- winsorized
  for (variable in settings$target_vars) {
    missing <- is.na(winsorized[[variable]])
    completed[[variable]][missing] <- imputed[[variable]][missing]
  }
  if (
    nrow(completed) != original_rows ||
    !identical(names(completed), original_names) ||
    !identical(raw$row_id_mf, completed$row_id_mf)
  ) {
    stop(cohort, " row or column identity changed.")
  }
  if (anyNA(completed[, settings$imputation_vars, drop = FALSE])) {
    stop(cohort, " completed data still contain missing imputation-model values.")
  }
  unchanged_predictors <- vapply(
    settings$predictor_only_vars,
    function(variable) equal_vector(standardized[[variable]], completed[[variable]]),
    logical(1)
  )
  if (!all(unchanged_predictors)) {
    stop(cohort, " predictor-only variables changed: ", paste(names(unchanged_predictors)[!unchanged_predictors], collapse = ", "))
  }
  openxlsx::write.xlsx(
    completed,
    file.path(settings$out_dir, output_file),
    overwrite = TRUE,
    keepNA = FALSE
  )
  invisible(completed)
}

impute_cohort(
  "CDIC",
  settings$cdic_input,
  c("row_id_mf", "patient_id_cleaned"),
  "CDIC_D1D3_covariate_missForest_completed_dataset_survival_used.xlsx"
)

impute_cohort(
  "MIMIC",
  settings$mimic_input,
  c("row_id_mf", "patient_id_cleaned"),
  "MIMIC_D1D3_covariate_missForest_completed_dataset_survival_used.xlsx"
)
