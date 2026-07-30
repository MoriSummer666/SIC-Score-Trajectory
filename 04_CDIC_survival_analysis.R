#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

settings <- list(
  cohort = "CDIC",
  input_file = file.path(
    "results",
    "03_random_forest_imputation",
    "CDIC_D1D3_covariate_missForest_completed_dataset_survival_used.xlsx"
  ),
  input_sheet = 1,
  posterior_file = file.path("results", "01_CDIC_trajectory_modeling", "CDIC_posterior_assignments.csv"),
  classification_file = file.path("results", "01_CDIC_trajectory_modeling", "CDIC_classification_summary.csv"),
  profile_file = file.path("results", "01_CDIC_trajectory_modeling", "CDIC_observed_trajectory_profiles.csv"),
  out_dir = file.path("results", "04_CDIC_survival_analysis"),
  landmark_day = 3,
  index14_day = 14,
  index28_day = 28,
  early_cut_day = 7,
  ties_method = "breslow",
  posterior_thresholds = c(0.70, 0.80),
  model1_covars = c("age", "sex_male", "bmi_kg_m2", "observed_cci"),
  model2_covars = c(
    "age", "sex_male", "bmi_kg_m2", "observed_cci",
    "temp_max_d1", "rr_max_d1", "hr_max_d1"
  ),
  model3_covars = c(
    "age", "sex_male", "bmi_kg_m2", "observed_cci",
    "temp_max_d1", "rr_max_d1", "hr_max_d1",
    "lactate_d1", "bun_d1", "ast_d1", "glucose_d1",
    "wbc_d1", "hemoglobin_d1"
  )
)

override <- getOption("sic.survival.settings")
if (is.list(override)) settings <- utils::modifyList(settings, override)

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (!startsWith(arg, "--")) next
  kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  key <- gsub("-", "_", kv[1])
  value <- if (length(kv) > 1L) paste(kv[-1], collapse = "=") else "TRUE"
  if (key %in% c(
    "cohort", "input_file", "input_sheet", "posterior_file",
    "classification_file", "profile_file", "out_dir"
  )) {
    settings[[key]] <- value
  }
}

required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "tibble",
  "stringr", "survival", "broom", "survRM2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(survival)
  library(broom)
  library(survRM2)
})

dir.create(settings$out_dir, recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

as_int01 <- function(x) {
  z <- as_num(x)
  out <- rep(NA_integer_, length(z))
  out[!is.na(z) & z %in% c(0, 1)] <- as.integer(z[!is.na(z) & z %in% c(0, 1)])
  out
}

id_to_character <- function(x) {
  if (is.numeric(x) || is.integer(x)) return(format(x, scientific = FALSE, trim = TRUE, digits = 22))
  trimws(as.character(x))
}

normalize_trajectory <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- ifelse(grepl("^[1-5]$", z), paste0("T", z), z)
  z[!z %in% paste0("T", 1:5)] <- NA_character_
  z
}

format_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

format_hr <- function(hr, low, high, p) {
  ifelse(
    is.na(hr),
    "Reference",
    sprintf("%.2f (%.2f-%.2f), p%s", hr, low, high, ifelse(p < 0.001, "<0.001", paste0("=", sprintf("%.3f", p))))
  )
}

format_rmst <- function(estimate, low, high, p) {
  sprintf(
    "%.3f (%.3f-%.3f), p%s",
    estimate,
    low,
    high,
    ifelse(p < 0.001, "<0.001", paste0("=", sprintf("%.3f", p)))
  )
}

read_table <- function(path, sheet = 1) {
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("xlsx", "xls")) {
    return(as.data.frame(read_excel(path, sheet = sheet, guess_max = 100000)))
  }
  if (extension == "csv") {
    return(as.data.frame(read_csv(path, show_col_types = FALSE, progress = FALSE)))
  }
  if (extension == "rds") {
    x <- readRDS(path)
    if (!is.data.frame(x)) stop("RDS input must contain a data.frame.")
    return(as.data.frame(x))
  }
  stop("Unsupported input extension: ", extension)
}

if (!file.exists(settings$input_file)) stop("Input file not found: ", settings$input_file)
raw <- read_table(settings$input_file, settings$input_sheet)
names(raw) <- trimws(names(raw))

trajectory_column <- intersect(c("trajectory", "trajectory_group"), names(raw))[1]
if (is.na(trajectory_column)) stop("Missing trajectory or trajectory_group.")

required_columns <- c(
  "patient_id_cleaned",
  "death_14d", "survival_time_14d_days",
  "death_28d", "survival_time_28d_days",
  settings$model3_covars
)
if (length(setdiff(required_columns, names(raw)))) {
  stop("Missing columns: ", paste(setdiff(required_columns, names(raw)), collapse = ", "))
}

df <- raw
df$patient_id_cleaned <- id_to_character(df$patient_id_cleaned)
numeric_columns <- unique(c(
  "death_14d", "survival_time_14d_days",
  "death_28d", "survival_time_28d_days",
  settings$model3_covars
))
for (variable in numeric_columns) df[[variable]] <- as_num(df[[variable]])
df$trajectory_group <- factor(normalize_trajectory(df[[trajectory_column]]), levels = paste0("T", 1:5))
df$death_14d <- as_int01(df$death_14d)
df$death_28d <- as_int01(df$death_28d)

if (anyDuplicated(df$patient_id_cleaned)) stop("patient_id_cleaned must be unique.")
if (any(is.na(df$trajectory_group))) stop("Invalid or missing trajectory group.")
if (any(is.na(df$death_14d)) || any(is.na(df$death_28d))) stop("Mortality outcomes must be coded 0/1.")
if (any(df$death_14d > df$death_28d)) stop("death_14d cannot exceed death_28d.")
if (any(df$survival_time_14d_days <= settings$landmark_day) ||
    any(df$survival_time_28d_days <= settings$landmark_day)) {
  stop("All patients must survive beyond the D3 landmark.")
}
if (any(df$survival_time_14d_days > settings$index14_day) ||
    any(df$survival_time_28d_days > settings$index28_day)) {
  stop("Follow-up time exceeds the prespecified horizon.")
}
if (any(!complete.cases(df[, settings$model3_covars, drop = FALSE]))) {
  stop("The missForest-completed dataset still contains missing Model III covariates.")
}
if (any(vapply(settings$model3_covars, function(v) length(unique(df[[v]])) <= 1L, logical(1)))) {
  stop("At least one Model III covariate is non-varying.")
}

tau14 <- settings$index14_day - settings$landmark_day
tau28 <- settings$index28_day - settings$landmark_day
df <- df |>
  mutate(
    time_14day_landmark = pmin(survival_time_14d_days - settings$landmark_day, tau14),
    event_14day_landmark = death_14d,
    time_28day_landmark = pmin(survival_time_28d_days - settings$landmark_day, tau28),
    event_28day_landmark = death_28d
  )

if (!all(file.exists(c(settings$classification_file, settings$profile_file)))) {
  stop("Classification summary or observed trajectory profile file not found.")
}
classification_summary <- read_table(settings$classification_file) |>
  mutate(trajectory_group = factor(normalize_trajectory(trajectory_group), levels = paste0("T", 1:5)))
observed_profile <- read_table(settings$profile_file) |>
  mutate(
    trajectory_group = factor(normalize_trajectory(trajectory_group), levels = paste0("T", 1:5)),
    time0 = as_num(time0),
    mean_sic = as_num(mean_sic)
  )
profile_wide <- observed_profile |>
  select(trajectory_group, time0, mean_sic) |>
  pivot_wider(names_from = time0, values_from = mean_sic, names_prefix = "Mean_SIC_D") |>
  transmute(
    trajectory_group,
    Mean_SIC_D1 = Mean_SIC_D0,
    Mean_SIC_D2 = Mean_SIC_D1,
    Mean_SIC_D3 = Mean_SIC_D2
  )
outcome_summary <- df |>
  group_by(trajectory_group) |>
  summarise(
    N = n(),
    Proportion = n() / nrow(df),
    Mortality_14day = mean(event_14day_landmark),
    Mortality_28day = mean(event_28day_landmark),
    .groups = "drop"
  ) |>
  left_join(
    classification_summary |>
      transmute(trajectory_group, APP = as_num(APP)),
    by = "trajectory_group"
  ) |>
  left_join(profile_wide, by = "trajectory_group") |>
  mutate(
    Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
    .before = 1
  )
outcome_summary_publication <- outcome_summary |>
  transmute(
    Trajectory = paste0(
      trajectory_group,
      " ",
      c(
        T1 = "Low-stable",
        T2 = "Early-rising",
        T3 = "Moderate-decreasing",
        T4 = "High-decreasing",
        T5 = "Persistent-high"
      )[as.character(trajectory_group)]
    ),
    `Class size, n (%)` = sprintf("%s (%.1f%%)", format(N, big.mark = ",", trim = TRUE), 100 * Proportion),
    APP = sprintf("%.3f", APP),
    `Mean SIC score on D1` = sprintf("%.2f", Mean_SIC_D1),
    `Mean SIC score on D2` = sprintf("%.2f", Mean_SIC_D2),
    `Mean SIC score on D3` = sprintf("%.2f", Mean_SIC_D3),
    `14-day mortality (%)` = sprintf("%.1f%%", 100 * Mortality_14day),
    `28-day mortality (%)` = sprintf("%.1f%%", 100 * Mortality_28day)
  )
summary_filename <- if (settings$cohort == "MIMIC") {
  "Table_S7_MIMIC_trajectory_summary.csv"
} else {
  "Table_S6_CDIC_trajectory_summary.csv"
}
write_csv(outcome_summary_publication, file.path(settings$out_dir, summary_filename), na = "")

model_sets <- list(
  Crude = character(0),
  `Model I` = settings$model1_covars,
  `Model II` = settings$model2_covars,
  `Model III` = settings$model3_covars
)

outcomes <- list(
  `14-day mortality` = list(time = "time_14day_landmark", event = "event_14day_landmark", tau = tau14),
  `28-day mortality` = list(time = "time_28day_landmark", event = "event_28day_landmark", tau = tau28)
)

extract_ph <- function(fit, outcome, model) {
  zph <- cox.zph(fit, transform = "km", terms = TRUE, singledf = FALSE)
  tab <- as.data.frame(zph$table)
  tab$Tested_term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab)[grep("^chisq$|chi", names(tab), ignore.case = TRUE)[1]] <- "Chi_square"
  names(tab)[grep("^df$", names(tab), ignore.case = TRUE)[1]] <- "df"
  names(tab)[grep("^p$|pvalue|p.value", names(tab), ignore.case = TRUE)[1]] <- "P_value"
  tab |>
    filter(Tested_term %in% c("trajectory_group", "GLOBAL")) |>
    transmute(
      Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
      Outcome = outcome,
      Model = model,
      Tested_term,
      Chi_square = as_num(Chi_square),
      df = as_num(df),
      P_value = as_num(P_value)
    )
}

fit_cox <- function(data, outcome, time_var, event_var, model, covars) {
  needed <- c(time_var, event_var, "trajectory_group", covars)
  analysis_data <- data[complete.cases(data[, needed, drop = FALSE]), , drop = FALSE]
  analysis_data$trajectory_group <- factor(analysis_data$trajectory_group, levels = paste0("T", 1:5))
  formula <- as.formula(
    paste0("Surv(", time_var, ", ", event_var, ") ~ ", paste(c("trajectory_group", covars), collapse = " + "))
  )
  fit <- coxph(
    formula,
    data = analysis_data,
    ties = settings$ties_method,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  effects <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(str_detect(term, "^trajectory_group")) |>
    transmute(
      Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
      Outcome = outcome,
      Trajectory = str_remove(term, "^trajectory_group"),
      Model = model,
      N = fit$n,
      Events = fit$nevent,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
  reference <- tibble(
    Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
    Outcome = outcome,
    Trajectory = "T1",
    Model = model,
    N = fit$n,
    Events = fit$nevent,
    HR = 1,
    CI_low = NA_real_,
    CI_high = NA_real_,
    P_value = NA_real_
  )
  list(results = bind_rows(reference, effects), ph = extract_ph(fit, outcome, model))
}

cox_runs <- list()
for (outcome in names(outcomes)) {
  info <- outcomes[[outcome]]
  for (model in names(model_sets)) {
    cox_runs[[paste(outcome, model, sep = "__")]] <- fit_cox(
      df,
      outcome,
      info$time,
      info$event,
      model,
      model_sets[[model]]
    )
  }
}

cox_results <- bind_rows(lapply(cox_runs, `[[`, "results")) |>
  mutate(
    Outcome = factor(Outcome, levels = names(outcomes)),
    Trajectory = factor(Trajectory, levels = paste0("T", 1:5)),
    Model = factor(Model, levels = names(model_sets))
  ) |>
  arrange(Outcome, Trajectory, Model)

ph_results <- bind_rows(lapply(cox_runs, `[[`, "ph")) |>
  mutate(
    Outcome = factor(Outcome, levels = names(outcomes)),
    Model = factor(Model, levels = names(model_sets)),
    Tested_term = factor(Tested_term, levels = c("trajectory_group", "GLOBAL"))
  ) |>
  arrange(Outcome, Model, Tested_term)

cox_filename <- if (settings$cohort == "MIMIC") {
  "Table_S10_MIMIC_Cox_models.csv"
} else {
  "Table_3_CDIC_Cox_models.csv"
}
trajectory_names <- c(
  T1 = "Low-stable",
  T2 = "Early-rising",
  T3 = "Moderate-decreasing",
  T4 = "High-decreasing",
  T5 = "Persistent-high"
)
cox_publication <- cox_results |>
  mutate(
    Trajectory = paste0(as.character(Trajectory), " ", trajectory_names[as.character(Trajectory)]),
    Result = format_hr(ifelse(HR == 1 & is.na(CI_low), NA_real_, HR), CI_low, CI_high, P_value)
  ) |>
  select(Outcome, Trajectory, Model, Result) |>
  pivot_wider(names_from = Model, values_from = Result) |>
  arrange(Outcome, factor(str_extract(Trajectory, "^T[1-5]"), levels = paste0("T", 1:5)))
write_csv(cox_publication, file.path(settings$out_dir, cox_filename), na = "")
write_csv(ph_results, file.path(settings$out_dir, "Table_S11_Schoenfeld_PH_tests.csv"), na = "")

fit_piecewise_pair <- function(target_group) {
  needed <- c(
    "trajectory_group", "time_28day_landmark", "event_28day_landmark",
    settings$model3_covars
  )
  pair <- df |>
    filter(trajectory_group %in% c("T1", target_group)) |>
    select(all_of(needed)) |>
    filter(complete.cases(.)) |>
    mutate(exposed = as.integer(trajectory_group == target_group))
  rhs <- paste(c("exposed", settings$model3_covars), collapse = " + ")
  early_data <- pair |>
    mutate(
      time_piece = pmin(time_28day_landmark, settings$early_cut_day),
      event_piece = as.integer(event_28day_landmark == 1 & time_28day_landmark <= settings$early_cut_day)
    )
  early_fit <- coxph(
    as.formula(paste0("Surv(time_piece, event_piece) ~ ", rhs)),
    data = early_data,
    ties = settings$ties_method
  )
  early <- tidy(early_fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "exposed") |>
    transmute(
      Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
      Contrast = paste0(target_group, " vs T1"),
      Time_window = paste0("0-", settings$early_cut_day, " days"),
      N = nrow(early_data),
      Events = sum(early_data$event_piece),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
  late_data <- pair |>
    filter(time_28day_landmark > settings$early_cut_day) |>
    mutate(
      time_piece = time_28day_landmark - settings$early_cut_day,
      event_piece = event_28day_landmark
    )
  late_fit <- coxph(
    as.formula(paste0("Surv(time_piece, event_piece) ~ ", rhs)),
    data = late_data,
    ties = settings$ties_method
  )
  late <- tidy(late_fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term == "exposed") |>
    transmute(
      Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
      Contrast = paste0(target_group, " vs T1"),
      Time_window = paste0(settings$early_cut_day, "-", tau28, " days"),
      N = nrow(late_data),
      Events = sum(late_data$event_piece),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
  bind_rows(early, late)
}

piecewise_results <- bind_rows(lapply(paste0("T", 2:5), fit_piecewise_pair))
piecewise_filename <- if (settings$cohort == "MIMIC") {
  "Table_S13_MIMIC_piecewise_Cox.csv"
} else {
  "Table_S12_CDIC_piecewise_Cox.csv"
}
piecewise_publication <- piecewise_results |>
  transmute(
    Contrast,
    `Time window after the D3 landmark` = str_replace_all(Time_window, " days", " d"),
    `HR (95% CI), P-value` = format_hr(HR, CI_low, CI_high, P_value)
  )
write_csv(piecewise_publication, file.path(settings$out_dir, piecewise_filename), na = "")

extract_rmst_difference <- function(result, adjustment, outcome, target, tau, n, events) {
  tab <- as.data.frame(result)
  tab$Metric <- rownames(tab)
  rownames(tab) <- NULL
  row <- tab |>
    filter(
      grepl("RMST", Metric, ignore.case = TRUE),
      grepl("arm=1", Metric, fixed = TRUE),
      grepl("arm=0", Metric, fixed = TRUE),
      grepl("-", Metric, fixed = TRUE)
    )
  if (!nrow(row)) stop("Could not identify the RMST-difference row.")
  find_column <- function(patterns) {
    hits <- unique(unlist(lapply(patterns, function(pattern) grep(pattern, names(row), ignore.case = TRUE, value = TRUE))))
    if (!length(hits)) NA_character_ else hits[1]
  }
  estimate_column <- find_column(c("^Est\\.$", "^Est$", "Estimate"))
  lower_column <- find_column("lower")
  upper_column <- find_column("upper")
  p_column <- find_column(c("^p$", "p.value", "p_value"))
  tibble(
    Cohort = ifelse(settings$cohort == "MIMIC", "MIMIC-IV", settings$cohort),
    Outcome = outcome,
    Contrast = paste0(target, " vs T1"),
    Tau = tau,
    Adjustment = adjustment,
    N = n,
    Events = events,
    RMST_difference_days = as_num(row[[estimate_column]][1]),
    CI_low = as_num(row[[lower_column]][1]),
    CI_high = as_num(row[[upper_column]][1]),
    P_value = as_num(row[[p_column]][1])
  )
}

run_rmst_pair <- function(target, outcome, time_var, event_var, tau) {
  needed <- c(time_var, event_var, "trajectory_group", settings$model3_covars)
  pair <- df |>
    filter(trajectory_group %in% c("T1", target)) |>
    select(all_of(needed)) |>
    filter(complete.cases(.)) |>
    mutate(arm = as.integer(trajectory_group == target))
  n_events <- sum(pair[[event_var]])
  unadjusted_fit <- rmst2(
    time = pair[[time_var]],
    status = pair[[event_var]],
    arm = pair$arm,
    tau = tau
  )
  covariate_matrix <- model.matrix(
    reformulate(settings$model3_covars),
    data = pair
  )
  covariate_matrix <- covariate_matrix[, colnames(covariate_matrix) != "(Intercept)", drop = FALSE]
  adjusted_fit <- rmst2(
    time = pair[[time_var]],
    status = pair[[event_var]],
    arm = pair$arm,
    tau = tau,
    covariates = covariate_matrix
  )
  bind_rows(
    extract_rmst_difference(
      unadjusted_fit$unadjusted.result,
      "Unadjusted",
      outcome,
      target,
      tau,
      nrow(pair),
      n_events
    ),
    extract_rmst_difference(
      adjusted_fit$adjusted.result,
      "Model III",
      outcome,
      target,
      tau,
      nrow(pair),
      n_events
    )
  )
}

rmst_results <- bind_rows(
  lapply(paste0("T", 2:5), function(target) {
    run_rmst_pair(target, "14-day mortality", "time_14day_landmark", "event_14day_landmark", tau14)
  }),
  lapply(paste0("T", 2:5), function(target) {
    run_rmst_pair(target, "28-day mortality", "time_28day_landmark", "event_28day_landmark", tau28)
  })
)
rmst_filename <- if (settings$cohort == "MIMIC") {
  "Table_S15_MIMIC_RMST.csv"
} else {
  "Table_S14_CDIC_RMST.csv"
}
rmst_publication <- rmst_results |>
  mutate(Result = format_rmst(RMST_difference_days, CI_low, CI_high, P_value)) |>
  select(Outcome, Contrast, Adjustment, Result) |>
  pivot_wider(names_from = Adjustment, values_from = Result) |>
  rename(
    `Unadjusted RMST difference (days)` = Unadjusted,
    `Model III-adjusted RMST difference (days)` = `Model III`
  )
write_csv(rmst_publication, file.path(settings$out_dir, rmst_filename), na = "")

if (!file.exists(settings$posterior_file)) stop("Posterior file not found: ", settings$posterior_file)
posterior <- read_table(settings$posterior_file)
names(posterior) <- trimws(names(posterior))
if (!all(c("patient_id_cleaned", "max_postprob") %in% names(posterior))) {
  stop("Posterior file must contain patient_id_cleaned and max_postprob.")
}
posterior <- posterior |>
  transmute(
    patient_id_cleaned = id_to_character(patient_id_cleaned),
    max_postprob = as_num(max_postprob)
  ) |>
  arrange(patient_id_cleaned, desc(max_postprob)) |>
  distinct(patient_id_cleaned, .keep_all = TRUE)

df_posterior <- df |>
  left_join(posterior, by = "patient_id_cleaned")
if (any(is.na(df_posterior$max_postprob))) stop("Not all survival records matched a posterior probability.")

run_posterior_threshold <- function(threshold) {
  selected <- df_posterior |>
    filter(max_postprob >= threshold)
  needed <- c(
    "time_28day_landmark", "event_28day_landmark",
    "trajectory_group", settings$model3_covars
  )
  selected <- selected[complete.cases(selected[, needed, drop = FALSE]), , drop = FALSE]
  fit <- coxph(
    as.formula(
      paste0(
        "Surv(time_28day_landmark, event_28day_landmark) ~ ",
        paste(c("trajectory_group", settings$model3_covars), collapse = " + ")
      )
    ),
    data = selected,
    ties = settings$ties_method
  )
  effects <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(str_detect(term, "^trajectory_group")) |>
    transmute(
      Trajectory = str_remove(term, "^trajectory_group"),
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
  counts <- selected |>
    group_by(trajectory_group) |>
    summarise(
      Retained_N = n(),
      Events = sum(event_28day_landmark),
      .groups = "drop"
    ) |>
    mutate(Trajectory = as.character(trajectory_group)) |>
    select(-trajectory_group) |>
    left_join(
      df |>
        count(trajectory_group, name = "Original_N") |>
        mutate(Trajectory = as.character(trajectory_group)) |>
        select(-trajectory_group),
      by = "Trajectory"
    ) |>
    mutate(Retained_percent = 100 * Retained_N / Original_N)
  bind_rows(
    tibble(Trajectory = "T1", HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P_value = NA_real_),
    effects
  ) |>
    left_join(counts, by = "Trajectory") |>
    transmute(
      `Sensitivity analysis` = paste0("Maximum posterior probability≥", threshold),
      Outcome = "28-day mortality",
      `Trajectory group` = Trajectory,
      `Patients retained, n (%)` = sprintf("%s (%.1f%%)", format(Retained_N, big.mark = ",", trim = TRUE), Retained_percent),
      `Events, n` = Events,
      `HR (95% CI), P-value` = format_hr(HR, CI_low, CI_high, P_value)
    )
}

posterior_results <- bind_rows(lapply(settings$posterior_thresholds, run_posterior_threshold))
posterior_filename <- if (settings$cohort == "MIMIC") {
  "Table_S17_MIMIC_high_posterior_subgroups.csv"
} else {
  "Table_S16_CDIC_high_posterior_subgroups.csv"
}
write_csv(posterior_results, file.path(settings$out_dir, posterior_filename), na = "")
