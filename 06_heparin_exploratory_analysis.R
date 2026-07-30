#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

settings <- list(
  input_file = file.path("data", "CDIC_MIMIC_pooled_heparin_timing_survival_analysis_dataset.xlsx"),
  input_sheet = "pooled_analysis_dataset",
  out_dir = file.path("results", "06_heparin_exploratory_analysis"),
  seed = 20260524L,
  early_cut_day = 7,
  tau14 = 11,
  tau28 = 25,
  heparin_levels = c("No_by_D3", "D1_0_24h", "D2D3_24_72h"),
  heparin_km_levels = c("D1_0_24h", "D2D3_24_72h", "No_by_D3"),
  heparin_labels = c(
    No_by_D3 = "No heparin by D3",
    D1_0_24h = "0-24h initiation",
    D2D3_24_72h = "24-72h initiation"
  ),
  heparin_colors = c(D1_0_24h = "#E41A1C", D2D3_24_72h = "#FF7F00", No_by_D3 = "#1F77B4"),
  trajectory_levels = paste0("T", 1:5),
  trajectory_names = c(
    T1 = "Low-stable",
    T2 = "Early-rising",
    T3 = "Moderate-decreasing",
    T4 = "High-decreasing",
    T5 = "Persistent-high"
  ),
  trajectory_colors = c(T1 = "#1F77B4", T2 = "#FF7F0E", T3 = "#2CA02C", T4 = "#9467BD", T5 = "#D62728"),
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
  ),
  balance_labels = c(
    age = "Age",
    sex_male = "Male sex",
    bmi_kg_m2 = "BMI",
    observed_cci = "Age-adjusted CCI",
    temp_max_d1 = "D1 maximum temperature",
    rr_max_d1 = "D1 maximum respiratory rate",
    hr_max_d1 = "D1 maximum heart rate",
    lactate_d1 = "D1 lactate",
    bun_d1 = "D1 BUN",
    ast_d1 = "D1 AST",
    glucose_d1 = "D1 glucose",
    wbc_d1 = "D1 WBC",
    hemoglobin_d1 = "D1 hemoglobin",
    cohort = "MIMIC-IV cohort indicator"
  )
)

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (!startsWith(arg, "--")) next
  kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  key <- gsub("-", "_", kv[1])
  value <- if (length(kv) > 1L) paste(kv[-1], collapse = "=") else "TRUE"
  if (key %in% c("input_file", "input_sheet", "out_dir")) settings[[key]] <- value
}

required_packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "tibble", "stringr",
  "survival", "broom", "ggplot2", "survminer", "scales",
  "survRM2", "WeightIt", "cobalt"
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
  library(ggplot2)
  library(survminer)
  library(scales)
  library(survRM2)
  library(WeightIt)
  library(cobalt)
})

set.seed(settings$seed)
dir.create(settings$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(settings$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

as_int01 <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), NA_integer_, as.integer(x)))
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(z))
  out[z %in% c("1", "yes", "y", "true", "t", "male", "m", "event", "dead")] <- 1L
  out[z %in% c("0", "no", "n", "false", "f", "female", "alive", "censored")] <- 0L
  numeric_value <- suppressWarnings(as.numeric(z))
  fallback <- is.na(out) & !is.na(numeric_value) & numeric_value %in% c(0, 1)
  out[fallback] <- as.integer(numeric_value[fallback])
  out
}

normalize_cohort <- function(x) {
  z <- gsub("[^A-Z0-9]", "", toupper(trimws(as.character(x))))
  ifelse(z == "CDIC", "CDIC", ifelse(z %in% c("MIMIC", "MIMICIV"), "MIMIC", NA_character_))
}

normalize_trajectory <- function(x) {
  str_extract(toupper(trimws(as.character(x))), "T[1-5]")
}

normalize_heparin <- function(x) {
  z <- trimws(as.character(x))
  u <- toupper(gsub("[^A-Z0-9]", "", z))
  case_when(
    u %in% c("NOBYD3", "NOHEPARINBYD3", "NOBEFORED3") ~ "No_by_D3",
    u %in% c("LE24H", "D1024H", "D10TO24H", "024H", "024") ~ "D1_0_24h",
    u %in% c("H2472H", "D2D32472H", "D2D324TO72H", "2472H", "2472") ~ "D2D3_24_72h",
    grepl("NO.*D3|NO.*HEPARIN", z, ignore.case = TRUE) ~ "No_by_D3",
    grepl("0\\s*[-_/]?\\s*24|LE24|D1", z, ignore.case = TRUE) ~ "D1_0_24h",
    grepl("24\\s*[-_/]?\\s*72|H24.*72|D2D3|D2-D3", z, ignore.case = TRUE) ~ "D2D3_24_72h",
    TRUE ~ NA_character_
  )
}

make_rhs <- function(terms) paste(unique(terms), collapse = " + ")

format_hr <- function(hr, low, high, p) {
  sprintf(
    "%.2f (%.2f-%.2f), p%s",
    hr,
    low,
    high,
    ifelse(p < 0.001, "<0.001", paste0("=", sprintf("%.3f", p)))
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

if (!file.exists(settings$input_file)) stop("Input file not found: ", settings$input_file)
raw <- read_excel(settings$input_file, sheet = settings$input_sheet, guess_max = 100000) |>
  as.data.frame(check.names = FALSE)

if (!"trajectory_group" %in% names(raw) && "trajectory" %in% names(raw)) raw$trajectory_group <- raw$trajectory
if (!"heparin_group_preD3" %in% names(raw) && "heparin_timing_common" %in% names(raw)) {
  raw$heparin_group_preD3 <- raw$heparin_timing_common
}
if (!"pooled_patient_id_cleaned" %in% names(raw)) {
  if (!all(c("cohort", "patient_id_cleaned") %in% names(raw))) stop("Missing pooled_patient_id_cleaned.")
  raw$pooled_patient_id_cleaned <- paste(raw$cohort, raw$patient_id_cleaned, sep = "_")
}
if (!"time_14day_landmark" %in% names(raw)) {
  raw$time_14day_landmark <- pmin(pmax(as_num(raw$survival_time_14d_days) - 3, 0), settings$tau14)
}
if (!"event_14day_landmark" %in% names(raw)) raw$event_14day_landmark <- as_int01(raw$death_14d)
if (!"time_28day_landmark" %in% names(raw)) {
  raw$time_28day_landmark <- pmin(pmax(as_num(raw$survival_time_28d_days) - 3, 0), settings$tau28)
}
if (!"event_28day_landmark" %in% names(raw)) raw$event_28day_landmark <- as_int01(raw$death_28d)

required_columns <- unique(c(
  "cohort", "pooled_patient_id_cleaned", "trajectory_group", "heparin_group_preD3",
  "time_14day_landmark", "event_14day_landmark",
  "time_28day_landmark", "event_28day_landmark",
  settings$model3_covars
))
if (length(setdiff(required_columns, names(raw)))) {
  stop("Missing columns: ", paste(setdiff(required_columns, names(raw)), collapse = ", "))
}

for (variable in setdiff(settings$model3_covars, "sex_male")) raw[[variable]] <- as_num(raw[[variable]])

d <- raw |>
  mutate(
    cohort = factor(normalize_cohort(cohort), levels = c("CDIC", "MIMIC")),
    pooled_patient_id_cleaned = as.character(pooled_patient_id_cleaned),
    trajectory_group = factor(normalize_trajectory(trajectory_group), levels = settings$trajectory_levels),
    heparin_timing = factor(normalize_heparin(heparin_group_preD3), levels = settings$heparin_levels),
    heparin_timing_km = factor(normalize_heparin(heparin_group_preD3), levels = settings$heparin_km_levels),
    sex_male = factor(as_int01(sex_male), levels = c(0, 1), labels = c("Female", "Male")),
    time_14day_landmark = as_num(time_14day_landmark),
    event_14day_landmark = as_int01(event_14day_landmark),
    time_28day_landmark = as_num(time_28day_landmark),
    event_28day_landmark = as_int01(event_28day_landmark)
  ) |>
  filter(
    if_all(all_of(required_columns), ~ !is.na(.x)),
    time_14day_landmark >= 0,
    time_14day_landmark <= settings$tau14,
    time_28day_landmark >= 0,
    time_28day_landmark <= settings$tau28
  ) |>
  droplevels()

if (anyDuplicated(d$pooled_patient_id_cleaned)) stop("pooled_patient_id_cleaned must be unique.")
if (nlevels(d$cohort) != 2L) stop("Both cohorts are required.")
if (nlevels(d$trajectory_group) != 5L) stop("All five trajectories are required.")
if (nlevels(d$heparin_timing) != 3L) stop("All three heparin groups are required.")
if (any(vapply(settings$model3_covars, function(v) length(unique(d[[v]])) <= 1L, logical(1)))) {
  stop("At least one Model III covariate is non-varying.")
}

outcomes <- list(
  `14-day mortality` = list(time = "time_14day_landmark", event = "event_14day_landmark", tau = settings$tau14),
  `28-day mortality` = list(time = "time_28day_landmark", event = "event_28day_landmark", tau = settings$tau28)
)

model_sets <- list(
  Crude = character(0),
  `Model I` = settings$model1_covars,
  `Model II` = settings$model2_covars,
  `Model III` = settings$model3_covars,
  `Model III + trajectory` = c(settings$model3_covars, "trajectory_group")
)

fit_heparin_cox <- function(data, outcome, time_var, event_var, model, covars, weights = NULL) {
  terms <- c("heparin_timing", covars, "strata(cohort)")
  if (!is.null(weights)) terms <- c(terms, "cluster(pooled_patient_id_cleaned)")
  formula <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", make_rhs(terms)))
  fit <- if (is.null(weights)) {
    coxph(formula, data = data, ties = "breslow", x = TRUE, y = TRUE)
  } else {
    coxph(
      formula,
      data = data,
      ties = "breslow",
      weights = data[[weights]],
      robust = TRUE,
      x = TRUE,
      y = TRUE
    )
  }
  results <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term %in% c("heparin_timingD1_0_24h", "heparin_timingD2D3_24_72h")) |>
    transmute(
      Outcome = outcome,
      Model = model,
      Comparison = ifelse(
        term == "heparin_timingD1_0_24h",
        "0-24 h initiation vs No heparin",
        "24-72 h initiation vs No heparin"
      ),
      N = fit$n,
      Events = fit$nevent,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    )
  list(fit = fit, results = results)
}

main_runs <- list()
for (outcome in names(outcomes)) {
  info <- outcomes[[outcome]]
  for (model in names(model_sets)) {
    main_runs[[paste(outcome, model, sep = "__")]] <- fit_heparin_cox(
      d,
      outcome,
      info$time,
      info$event,
      model,
      model_sets[[model]]
    )
  }
}

standardize_for_ps <- function(data, covars) {
  output <- data
  modeled <- covars
  for (i in seq_along(covars)) {
    variable <- covars[i]
    value <- output[[variable]]
    numeric_value <- as_num(value)
    if ((is.numeric(value) || is.integer(value)) &&
        length(unique(na.omit(numeric_value))) > 5L &&
        sd(numeric_value) > 0) {
      new_name <- paste0("z__", variable)
      output[[new_name]] <- as.numeric(scale(numeric_value))
      modeled[i] <- new_name
    }
  }
  list(data = output, covars = modeled)
}

ps_data <- d[complete.cases(d[, unique(c(
  "heparin_timing", "cohort", "pooled_patient_id_cleaned",
  "time_14day_landmark", "event_14day_landmark",
  "time_28day_landmark", "event_28day_landmark",
  settings$model3_covars
)), drop = FALSE]), , drop = FALSE] |>
  droplevels()
standardized <- standardize_for_ps(ps_data, settings$model3_covars)
ps_data <- standardized$data
ps_formula <- as.formula(paste0(
  "heparin_timing ~ ",
  make_rhs(c(standardized$covars, "cohort"))
))
weight_model <- weightit(
  ps_formula,
  data = ps_data,
  method = "glm",
  estimand = "ATO"
)
ps_data$w_ato <- as.numeric(weight_model$weights)
ps_data <- ps_data[is.finite(ps_data$w_ato) & ps_data$w_ato > 0, , drop = FALSE]

ato_runs <- lapply(names(outcomes), function(outcome) {
  info <- outcomes[[outcome]]
  fit_heparin_cox(
    ps_data,
    outcome,
    info$time,
    info$event,
    "ATO weighted",
    character(0),
    weights = "w_ato"
  )
})
names(ato_runs) <- names(outcomes)

main_results <- bind_rows(
  lapply(main_runs, `[[`, "results"),
  lapply(ato_runs, `[[`, "results")
) |>
  mutate(
    Outcome = factor(Outcome, levels = names(outcomes)),
    Model = factor(Model, levels = c(names(model_sets), "ATO weighted")),
    Comparison = factor(
      Comparison,
      levels = c("0-24 h initiation vs No heparin", "24-72 h initiation vs No heparin")
    )
  ) |>
  arrange(Outcome, Model, Comparison)

main_publication <- main_results |>
  mutate(Result = format_hr(HR, CI_low, CI_high, P_value)) |>
  select(Outcome, Model, Comparison, Result) |>
  pivot_wider(names_from = Comparison, values_from = Result)
write_csv(main_publication, file.path(settings$out_dir, "Table_S19_heparin_pooled_Cox.csv"), na = "")

balance_object <- bal.tab(
  weight_model,
  un = TRUE,
  pairwise = TRUE,
  multi.summary = TRUE,
  abs = TRUE,
  continuous = "std",
  binary = "std",
  quick = FALSE
)

pick_column <- function(names_value, exact, pattern) {
  direct <- exact[exact %in% names_value]
  if (length(direct)) return(direct[1])
  hits <- grep(pattern, names_value, value = TRUE, ignore.case = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

pair_balance_fallback <- function(balance_object) {
  pair_names <- names(balance_object$Pair.Balance)
  rows <- bind_rows(lapply(seq_along(balance_object$Pair.Balance), function(i) {
    tab <- as.data.frame(balance_object$Pair.Balance[[i]]$Balance, check.names = FALSE) |>
      rownames_to_column("variable")
    unweighted <- pick_column(names(tab), c("Diff.Un"), "^Diff.*Un$")
    weighted <- pick_column(names(tab), c("Diff.Adj"), "^Diff.*Adj$")
    if (is.na(unweighted) || is.na(weighted)) return(tibble())
    tibble(
      variable = tab$variable,
      pair = if (length(pair_names)) pair_names[i] else paste0("pair_", i),
      unweighted = abs(as_num(tab[[unweighted]])),
      weighted = abs(as_num(tab[[weighted]]))
    )
  }))
  rows |>
    group_by(variable) |>
    summarise(
      Max_Diff_Unweighted = max(unweighted, na.rm = TRUE),
      Max_Diff_ATO = max(weighted, na.rm = TRUE),
      .groups = "drop"
    )
}

if (!is.null(balance_object$Balance.Across.Pairs)) {
  balance_raw <- as.data.frame(balance_object$Balance.Across.Pairs, check.names = FALSE) |>
    rownames_to_column("variable")
  unweighted_column <- pick_column(
    names(balance_raw),
    c("Max.Diff.Un", "Diff.Un", "M.Diff.Un", "Mean.Diff.Un"),
    "^Max.*Diff.*Un$|^Diff.*Un$"
  )
  weighted_column <- pick_column(
    names(balance_raw),
    c("Max.Diff.Adj", "Diff.Adj", "M.Diff.Adj", "Mean.Diff.Adj"),
    "^Max.*Diff.*Adj$|^Diff.*Adj$"
  )
  balance_data <- balance_raw |>
    transmute(
      variable,
      Max_Diff_Unweighted = abs(as_num(.data[[unweighted_column]])),
      Max_Diff_ATO = abs(as_num(.data[[weighted_column]]))
    )
} else {
  balance_data <- pair_balance_fallback(balance_object)
}

identify_original <- function(x) {
  z <- sub("^z__", "", gsub("`", "", x, fixed = TRUE))
  candidates <- c(settings$model3_covars, "cohort")
  candidates <- candidates[order(nchar(candidates), decreasing = TRUE)]
  out <- rep(NA_character_, length(z))
  for (i in seq_along(z)) {
    hits <- candidates[vapply(candidates, function(v) {
      identical(z[i], v) || startsWith(z[i], paste0(v, "_")) ||
        startsWith(z[i], paste0(v, ".")) || startsWith(z[i], paste0(v, ":")) ||
        (v %in% c("sex_male", "cohort") && startsWith(z[i], v))
    }, logical(1))]
    if (length(hits)) out[i] <- hits[1]
  }
  out
}

balance_plot_data <- balance_data |>
  mutate(
    original_variable = identify_original(variable),
    label = unname(settings$balance_labels[original_variable])
  ) |>
  filter(!is.na(label), is.finite(Max_Diff_Unweighted), is.finite(Max_Diff_ATO)) |>
  group_by(original_variable, label) |>
  summarise(
    Max_Diff_Unweighted = max(Max_Diff_Unweighted),
    Max_Diff_ATO = max(Max_Diff_ATO),
    .groups = "drop"
  ) |>
  arrange(desc(Max_Diff_Unweighted)) |>
  mutate(label = factor(label, levels = rev(label)))

balance_long <- bind_rows(
  balance_plot_data |> transmute(label, Sample = "Unweighted", SMD = Max_Diff_Unweighted),
  balance_plot_data |> transmute(label, Sample = "Overlap weighted", SMD = Max_Diff_ATO)
) |>
  mutate(Sample = factor(Sample, levels = c("Unweighted", "Overlap weighted")))

p_balance <- ggplot() +
  geom_segment(
    data = balance_plot_data,
    aes(x = Max_Diff_Unweighted, xend = Max_Diff_ATO, y = label, yend = label),
    color = "grey68",
    linewidth = 0.8
  ) +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "grey35") +
  geom_point(
    data = balance_long,
    aes(SMD, label, color = Sample, shape = Sample),
    size = 4
  ) +
  scale_color_manual(values = c(Unweighted = "#1F77B4", `Overlap weighted` = "#FF7F0E")) +
  scale_shape_manual(values = c(Unweighted = 16, `Overlap weighted` = 15)) +
  labs(
    x = "Maximum absolute standardized mean difference",
    y = NULL,
    color = NULL,
    shape = NULL,
    title = "Covariate balance before and after overlap weighting"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(settings$out_dir, "figures", "Figure_S2_ATO_covariate_balance.pdf"),
  p_balance,
  width = 11.2,
  height = max(7.2, 2.3 + 0.46 * nrow(balance_plot_data))
)

extract_lrt_p <- function(x) {
  tab <- as.data.frame(x)
  column <- grep("P\\(|Pr\\(", names(tab), value = TRUE)[1]
  tail(na.omit(as_num(tab[[column]])), 1)
}

find_interaction_term <- function(coefficient_names, exposure, trajectory) {
  candidates <- c(
    paste0(exposure, ":trajectory_group", trajectory),
    paste0("trajectory_group", trajectory, ":", exposure)
  )
  hit <- candidates[candidates %in% coefficient_names]
  if (!length(hit)) NA_character_ else hit[1]
}

run_interaction <- function(outcome, time_var, event_var) {
  needed <- unique(c(
    time_var, event_var, "heparin_timing", "trajectory_group",
    "cohort", settings$model3_covars
  ))
  analysis_data <- d[complete.cases(d[, needed, drop = FALSE]), , drop = FALSE] |>
    droplevels()
  main_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ ",
    make_rhs(c("heparin_timing", "trajectory_group", settings$model3_covars, "strata(cohort)"))
  ))
  interaction_formula <- as.formula(paste0(
    "Surv(", time_var, ", ", event_var, ") ~ ",
    make_rhs(c("heparin_timing * trajectory_group", settings$model3_covars, "strata(cohort)"))
  ))
  main_fit <- coxph(main_formula, data = analysis_data, ties = "breslow", x = TRUE, y = TRUE)
  interaction_fit <- coxph(interaction_formula, data = analysis_data, ties = "breslow", x = TRUE, y = TRUE)
  interaction_p <- extract_lrt_p(anova(main_fit, interaction_fit, test = "LRT"))
  coefficients <- coef(interaction_fit)
  covariance <- vcov(interaction_fit)
  exposure_definitions <- c(
    "0-24 h initiation vs No heparin" = "heparin_timingD1_0_24h",
    "24-72 h initiation vs No heparin" = "heparin_timingD2D3_24_72h"
  )
  bind_rows(lapply(settings$trajectory_levels, function(trajectory) {
    bind_rows(lapply(names(exposure_definitions), function(comparison) {
      exposure <- exposure_definitions[[comparison]]
      contrast <- setNames(rep(0, length(coefficients)), names(coefficients))
      contrast[exposure] <- 1
      if (trajectory != "T1") {
        interaction_term <- find_interaction_term(names(coefficients), exposure, trajectory)
        if (is.na(interaction_term)) stop("Interaction coefficient not found.")
        contrast[interaction_term] <- 1
      }
      log_hr <- sum(contrast * coefficients)
      standard_error <- sqrt(as.numeric(t(contrast) %*% covariance %*% contrast))
      tibble(
        Outcome = outcome,
        Trajectory = trajectory,
        Trajectory_label = unname(settings$trajectory_names[trajectory]),
        Comparison = comparison,
        N = nrow(analysis_data),
        Events = sum(analysis_data[[event_var]]),
        HR = exp(log_hr),
        CI_low = exp(log_hr - 1.96 * standard_error),
        CI_high = exp(log_hr + 1.96 * standard_error),
        P_value = 2 * pnorm(-abs(log_hr / standard_error)),
        Global_interaction_P = interaction_p
      )
    }))
  }))
}

interaction_results <- bind_rows(lapply(names(outcomes), function(outcome) {
  info <- outcomes[[outcome]]
  run_interaction(outcome, info$time, info$event)
}))
interaction_publication <- interaction_results |>
  group_by(Outcome) |>
  arrange(Outcome, factor(Trajectory, levels = settings$trajectory_levels), Comparison, .by_group = TRUE) |>
  mutate(
    `Global interaction P` = ifelse(row_number() == 1L, sprintf("%.3f", Global_interaction_P), ""),
    `Simple-effect HR (95% CI), P-value` = format_hr(HR, CI_low, CI_high, P_value)
  ) |>
  ungroup() |>
  transmute(
    Outcome,
    `Trajectory group` = paste0(Trajectory, " ", Trajectory_label),
    `Heparin timing comparison` = Comparison,
    `Simple-effect HR (95% CI), P-value`,
    `Global interaction P`
  )
write_csv(interaction_publication, file.path(settings$out_dir, "Table_S22_heparin_trajectory_interaction.csv"), na = "")

run_piecewise <- function() {
  needed <- unique(c(
    "time_28day_landmark", "event_28day_landmark",
    "heparin_timing", "cohort", "pooled_patient_id_cleaned",
    settings$model3_covars
  ))
  analysis_data <- d[complete.cases(d[, needed, drop = FALSE]), , drop = FALSE] |>
    droplevels()
  split_data <- survSplit(
    Surv(time_28day_landmark, event_28day_landmark) ~ .,
    data = analysis_data,
    cut = settings$early_cut_day,
    episode = "period",
    start = "tstart",
    end = "tstop",
    event = "event_piece"
  ) |>
    mutate(period = factor(period, levels = c(1, 2), labels = c("0-7 days", "7-25 days")))
  bind_rows(lapply(levels(split_data$period), function(period) {
    part <- split_data |> filter(.data$period == .env$period) |> droplevels()
    formula <- as.formula(paste0(
      "Surv(tstart, tstop, event_piece) ~ ",
      make_rhs(c("heparin_timing", settings$model3_covars, "strata(cohort)", "cluster(pooled_patient_id_cleaned)"))
    ))
    fit <- coxph(formula, data = part, ties = "breslow", robust = TRUE, x = TRUE, y = TRUE)
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term %in% c("heparin_timingD1_0_24h", "heparin_timingD2D3_24_72h")) |>
      transmute(
        Time_window = period,
        Comparison = ifelse(
          term == "heparin_timingD1_0_24h",
          "0-24 h initiation vs No heparin",
          "24-72 h initiation vs No heparin"
        ),
        N = fit$n,
        Events = fit$nevent,
        HR = estimate,
        CI_low = conf.low,
        CI_high = conf.high,
        P_value = p.value
      )
  }))
}

piecewise_results <- run_piecewise()
piecewise_publication <- piecewise_results |>
  transmute(
    `Time window after the D3 landmark` = Time_window,
    Comparison,
    `Analysis sample, n` = N,
    `Events, n` = Events,
    `HR (95% CI), P-value` = format_hr(HR, CI_low, CI_high, P_value)
  )
write_csv(piecewise_publication, file.path(settings$out_dir, "Table_S21_heparin_piecewise_Cox.csv"), na = "")

extract_rmst_difference <- function(result, outcome, comparison, adjustment, tau, n, events) {
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
  find_column <- function(patterns) {
    hits <- unique(unlist(lapply(patterns, function(pattern) grep(pattern, names(row), ignore.case = TRUE, value = TRUE))))
    if (!length(hits)) NA_character_ else hits[1]
  }
  estimate_column <- find_column(c("^Est\\.$", "^Est$", "Estimate"))
  lower_column <- find_column("lower")
  upper_column <- find_column("upper")
  p_column <- find_column(c("^p$", "p.value", "p_value"))
  tibble(
    Outcome = outcome,
    Comparison = comparison,
    Adjustment = adjustment,
    Tau = tau,
    N = n,
    Events = events,
    RMST_difference_days = as_num(row[[estimate_column]][1]),
    CI_low = as_num(row[[lower_column]][1]),
    CI_high = as_num(row[[upper_column]][1]),
    P_value = as_num(row[[p_column]][1])
  )
}

run_rmst <- function(outcome, time_var, event_var, tau, comparison_level) {
  covars <- c(settings$model3_covars, "cohort")
  analysis_data <- d |>
    filter(heparin_timing %in% c("No_by_D3", comparison_level)) |>
    droplevels()
  needed <- c(time_var, event_var, "heparin_timing", covars)
  analysis_data <- analysis_data[complete.cases(analysis_data[, needed, drop = FALSE]), , drop = FALSE]
  analysis_data$arm <- as.integer(analysis_data$heparin_timing == comparison_level)
  comparison <- paste0(unname(settings$heparin_labels[comparison_level]), " vs No heparin")
  unadjusted_fit <- rmst2(
    time = analysis_data[[time_var]],
    status = analysis_data[[event_var]],
    arm = analysis_data$arm,
    tau = tau
  )
  covariate_matrix <- model.matrix(reformulate(covars), data = analysis_data)
  covariate_matrix <- covariate_matrix[, colnames(covariate_matrix) != "(Intercept)", drop = FALSE]
  adjusted_fit <- rmst2(
    time = analysis_data[[time_var]],
    status = analysis_data[[event_var]],
    arm = analysis_data$arm,
    tau = tau,
    covariates = covariate_matrix
  )
  bind_rows(
    extract_rmst_difference(
      unadjusted_fit$unadjusted.result,
      outcome,
      comparison,
      "Unadjusted",
      tau,
      nrow(analysis_data),
      sum(analysis_data[[event_var]])
    ),
    extract_rmst_difference(
      adjusted_fit$adjusted.result,
      outcome,
      comparison,
      "Model III + cohort source",
      tau,
      nrow(analysis_data),
      sum(analysis_data[[event_var]])
    )
  )
}

rmst_results <- bind_rows(lapply(names(outcomes), function(outcome) {
  info <- outcomes[[outcome]]
  bind_rows(
    run_rmst(outcome, info$time, info$event, info$tau, "D1_0_24h"),
    run_rmst(outcome, info$time, info$event, info$tau, "D2D3_24_72h")
  )
}))
rmst_publication <- rmst_results |>
  transmute(
    Outcome,
    Comparison,
    Adjustment,
    `RMST difference, days (95% CI), P value` = format_rmst(
      RMST_difference_days,
      CI_low,
      CI_high,
      P_value
    )
  )
write_csv(rmst_publication, file.path(settings$out_dir, "Table_S20_heparin_RMST.csv"), na = "")

km_data <- d |>
  transmute(
    time = time_28day_landmark,
    event = event_28day_landmark,
    heparin_timing_km
  )
survival_fit <- survfit(Surv(time, event) ~ heparin_timing_km, data = km_data)
logrank <- survdiff(Surv(time, event) ~ heparin_timing_km, data = km_data)
logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1L, lower.tail = FALSE)
km_plot <- ggsurvplot(
  survival_fit,
  data = km_data,
  risk.table = TRUE,
  pval = ifelse(logrank_p < 0.001, "Log-rank p < 0.001", sprintf("Log-rank p = %.3f", logrank_p)),
  conf.int = FALSE,
  censor = TRUE,
  palette = unname(settings$heparin_colors[settings$heparin_km_levels]),
  xlab = "Days after D3 landmark",
  ylab = "Survival probability",
  xlim = c(0, settings$tau28),
  break.time.by = 5,
  legend = "bottom",
  legend.title = NULL,
  legend.labs = unname(settings$heparin_labels[settings$heparin_km_levels]),
  ggtheme = theme_bw(base_size = 12),
  tables.theme = theme_cleantable(base_size = 11),
  risk.table.height = 0.28,
  risk.table.y.text = FALSE,
  risk.table.y.text.col = TRUE
)
km_plot$plot <- km_plot$plot +
  coord_cartesian(ylim = c(0.75, 1)) +
  scale_y_continuous(limits = c(0.75, 1), breaks = seq(0.75, 1, 0.05), labels = percent_format(accuracy = 1))
pdf(file.path(settings$out_dir, "figures", "Figure_S8_heparin_KM_28day.pdf"), width = 9.5, height = 7.6)
print(km_plot, newpage = FALSE)
dev.off()

forest_plot <- function(data, title, color_variable, colors, output_file, height) {
  plot_data <- data |>
    filter(is.finite(HR), is.finite(CI_low), is.finite(CI_high)) |>
    mutate(
      Row = factor(paste(Comparison, color_variable, sep = " | "), levels = rev(unique(paste(Comparison, color_variable, sep = " | ")))),
      Label = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high)
    )
  xmax <- max(plot_data$CI_high) * 1.8
  p <- ggplot(plot_data, aes(HR, Row, color = .data[[color_variable]])) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18, linewidth = 0.75) +
    geom_point(shape = 18, size = 3.2) +
    geom_text(aes(x = xmax, label = Label), hjust = 1, color = "black", size = 3.5) +
    scale_x_log10() +
    scale_color_manual(values = colors) +
    coord_cartesian(xlim = c(min(plot_data$CI_low) / 1.2, xmax), clip = "off") +
    labs(x = "Hazard ratio (log scale)", y = NULL, color = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  ggsave(file.path(settings$out_dir, "figures", output_file), p, width = 10.5, height = height)
}

for (outcome in names(outcomes)) {
  figure_number <- ifelse(outcome == "14-day mortality", "4A", "4B")
  overall_data <- main_results |>
    filter(
      Outcome == outcome,
      Model %in% c("Model III", "Model III + trajectory", "ATO weighted")
    ) |>
    mutate(Model = droplevels(Model))
  forest_plot(
    overall_data,
    paste0(outcome, ": pooled heparin timing associations"),
    "Model",
    c(`Model III` = "#D62728", `Model III + trajectory` = "#F28E2B", `ATO weighted` = "#1F77B4"),
    paste0("Figure_", figure_number, "_heparin_overall_forest.pdf"),
    5.9
  )
}

for (outcome in names(outcomes)) {
  figure_number <- ifelse(outcome == "14-day mortality", "4C", "4D")
  trajectory_data <- interaction_results |>
    filter(Outcome == outcome) |>
    mutate(
      Trajectory = factor(Trajectory, levels = settings$trajectory_levels),
      Comparison = factor(
        Comparison,
        levels = c("0-24 h initiation vs No heparin", "24-72 h initiation vs No heparin")
      )
    )
  forest_plot(
    trajectory_data,
    paste0(outcome, ": associations within SIC trajectories"),
    "Trajectory",
    settings$trajectory_colors,
    paste0("Figure_", figure_number, "_heparin_trajectory_forest.pdf"),
    8
  )
}

forest_plot(
  piecewise_results,
  "Heparin timing and 28-day mortality by time window",
  "Time_window",
  c(`0-7 days` = "#D62728", `7-25 days` = "#1F77B4"),
  "Figure_S9_heparin_piecewise_forest.pdf",
  5.6
)
