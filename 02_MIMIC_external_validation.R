#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

settings <- list(
  input_file = file.path("data", "MIMIC_D1D3_SIC_lcmm.xlsx"),
  input_sheet = "Sheet1",
  cdic_dir = file.path("results", "01_CDIC_trajectory_modeling"),
  out_dir = file.path("results", "02_MIMIC_external_validation"),
  id_var = "patient_id_cleaned",
  sic_vars = c("sic_d1", "sic_d2", "sic_d3"),
  missing_tokens = c("", "NA", "N/A", ".", "NULL", "NaN"),
  missing_codes = c("20", "35"),
  min_observed = 2L,
  final_ng = 5L,
  raw_to_trajectory = c(`1` = "T4", `2` = "T2", `3` = "T1", `4` = "T3", `5` = "T5"),
  trajectory_names = c(
    T1 = "Low-stable",
    T2 = "Early-rising",
    T3 = "Moderate-decreasing",
    T4 = "High-decreasing",
    T5 = "Persistent-high"
  ),
  colors = c(T1 = "#1F77B4", T2 = "#FF7F0E", T3 = "#2CA02C", T4 = "#9467BD", T5 = "#D62728")
)

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (!startsWith(arg, "--")) next
  kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
  key <- gsub("-", "_", kv[1])
  value <- if (length(kv) > 1L) paste(kv[-1], collapse = "=") else "TRUE"
  if (key %in% c("input_file", "input_sheet", "cdic_dir", "out_dir")) settings[[key]] <- value
}

required_packages <- c("lcmm", "readxl", "dplyr", "tidyr", "tibble", "ggplot2", "readr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing R packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(lcmm)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  R_MAX_NUM_THREADS = "1"
)

dir.create(settings$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(settings$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

clean_sic <- function(x, variable) {
  z <- trimws(as.character(x))
  missing <- is.na(z) | toupper(z) %in% toupper(settings$missing_tokens) | z %in% settings$missing_codes
  value <- as_num(z)
  invalid_text <- !missing & is.na(value)
  invalid_numeric <- !missing & !is.na(value) &
    (value < 0 | value > 6 | abs(value - round(value)) > sqrt(.Machine$double.eps))
  if (any(invalid_text | invalid_numeric)) {
    stop("Unexpected SIC value in ", variable, ": ", paste(head(unique(z[invalid_text | invalid_numeric]), 10), collapse = ", "))
  }
  value[missing] <- NA_real_
  value
}

repair_model <- function(model) {
  model$call$fixed <- as.formula("sic_score ~ time0", env = .GlobalEnv)
  model$call$mixture <- as.formula("~ time0", env = .GlobalEnv)
  model$call$random <- as.formula("~ 1", env = .GlobalEnv)
  model$call$classmb <- as.formula("~ 1", env = .GlobalEnv)
  model$call$subject <- "subject_id_num"
  model$call$ng <- settings$final_ng
  model$call$nwg <- FALSE
  model$call$na.action <- 1
  model$call$var.time <- "time0"
  model$call$data <- quote(newdata)
  model
}

standardize_prediction <- function(x) {
  dat <- as.data.frame(x, check.names = FALSE)
  names_dat <- names(dat)
  subject_candidates <- c("subject_id_num", "subject", "id")
  subject_col <- subject_candidates[subject_candidates %in% names_dat][1]
  if (is.na(subject_col)) subject_col <- names_dat[1]
  class_candidates <- c("class", "raw_class", "pred_class", "predicted_class")
  class_col <- class_candidates[class_candidates %in% names_dat][1]
  prob_cols <- grep("^prob[._]?[0-9]+$", names_dat, value = TRUE, ignore.case = TRUE)
  if (length(prob_cols)) {
    suffix <- as.integer(gsub("[^0-9]", "", prob_cols))
    prob_cols <- prob_cols[order(suffix)]
  }
  if (length(prob_cols) < settings$final_ng) {
    remaining <- setdiff(names_dat, c(subject_col, class_col))
    prob_cols <- remaining[vapply(remaining, function(v) {
      z <- as_num(dat[[v]])
      all(is.na(z) | (z >= -1e-8 & z <= 1 + 1e-8))
    }, logical(1))]
  }
  if (length(prob_cols) < settings$final_ng) stop("Could not identify posterior-probability columns.")
  p <- as.data.frame(lapply(dat[, prob_cols[seq_len(settings$final_ng)], drop = FALSE], as_num))
  names(p) <- paste0("prob", seq_len(settings$final_ng))
  bind_cols(
    tibble(
      subject_id_num = as.integer(as_num(dat[[subject_col]])),
      raw_class_reported = if (!is.na(class_col)) as.integer(as_num(dat[[class_col]])) else NA_integer_
    ),
    p
  )
}

raw <- read_excel(
  settings$input_file,
  sheet = settings$input_sheet,
  col_types = "text",
  trim_ws = TRUE
) |>
  as.data.frame(check.names = FALSE)

required_columns <- c(settings$id_var, settings$sic_vars)
if (length(setdiff(required_columns, names(raw)))) {
  stop("Missing columns: ", paste(setdiff(required_columns, names(raw)), collapse = ", "))
}

raw[[settings$id_var]] <- trimws(as.character(raw[[settings$id_var]]))
if (any(is.na(raw[[settings$id_var]]) | raw[[settings$id_var]] == "")) stop("Missing patient_id_cleaned.")
if (anyDuplicated(raw[[settings$id_var]])) stop("patient_id_cleaned must be unique.")

wide <- raw
for (variable in settings$sic_vars) wide[[variable]] <- clean_sic(raw[[variable]], variable)
wide$n_sic_observed <- rowSums(!is.na(wide[, settings$sic_vars, drop = FALSE]))
if (any(wide$n_sic_observed < settings$min_observed)) {
  stop("Every patient must have at least ", settings$min_observed, " valid SIC observations.")
}
wide$subject_id_num <- seq_len(nrow(wide))

long <- wide |>
  select(subject_id_num, all_of(settings$id_var), n_sic_observed, all_of(settings$sic_vars)) |>
  pivot_longer(all_of(settings$sic_vars), names_to = "time_point", values_to = "sic_score") |>
  mutate(
    time0 = match(time_point, settings$sic_vars) - 1L,
    time_day = time0 + 1L
  ) |>
  filter(!is.na(sic_score)) |>
  arrange(subject_id_num, time0)

model_path <- file.path(settings$cdic_dir, "CDIC_locked_ng5_model.rds")
cdic_assignment_path <- file.path(settings$cdic_dir, "CDIC_posterior_assignments.csv")
cdic_profile_path <- file.path(settings$cdic_dir, "CDIC_observed_trajectory_profiles.csv")
cdic_fitted_path <- file.path(settings$cdic_dir, "CDIC_fitted_trajectory_profiles.csv")
if (!all(file.exists(c(model_path, cdic_assignment_path, cdic_profile_path, cdic_fitted_path)))) {
  stop("Run 01_CDIC_trajectory_modeling.R before external validation.")
}

model <- repair_model(readRDS(model_path))
if (is.null(model$conv) || model$conv != 1L) stop("The locked CDIC model is not strictly converged.")

external_raw <- local({
  ng <- settings$final_ng
  fixed_formula <- sic_score ~ time0
  mixture_formula <- ~ time0
  random_formula <- ~ 1
  predictClass(model, newdata = long, subject = "subject_id_num")
})

posterior <- standardize_prediction(external_raw) |>
  arrange(subject_id_num)

probability_columns <- paste0("prob", seq_len(settings$final_ng))
probability_matrix <- as.matrix(posterior[, probability_columns, drop = FALSE])
if (any(!is.finite(probability_matrix))) stop("Non-finite posterior probability.")
if (any(probability_matrix < -1e-8 | probability_matrix > 1 + 1e-8)) stop("Posterior probability outside [0,1].")
if (max(abs(rowSums(probability_matrix) - 1)) > 1e-4) stop("Posterior probabilities do not sum to one.")
probability_matrix <- probability_matrix / rowSums(probability_matrix)
posterior[probability_columns] <- as.data.frame(probability_matrix)
posterior$raw_class <- max.col(probability_matrix, ties.method = "first")
posterior$max_postprob <- apply(probability_matrix, 1, max)
posterior$assigned_prob <- probability_matrix[cbind(seq_len(nrow(probability_matrix)), posterior$raw_class)]

if (nrow(posterior) != nrow(wide) || !setequal(posterior$subject_id_num, wide$subject_id_num)) {
  stop("External-classification subjects do not match the eligible MIMIC cohort.")
}

posterior <- posterior |>
  left_join(
    wide |> select(subject_id_num, all_of(settings$id_var), n_sic_observed),
    by = "subject_id_num"
  ) |>
  mutate(
    trajectory_group = unname(settings$raw_to_trajectory[as.character(raw_class)]),
    trajectory_group = factor(trajectory_group, levels = paste0("T", 1:5)),
    trajectory_label = unname(settings$trajectory_names[as.character(trajectory_group)])
  )

for (group in paste0("T", 1:5)) {
  raw_class <- as.integer(names(settings$raw_to_trajectory)[settings$raw_to_trajectory == group])
  posterior[[paste0("prob_", group)]] <- posterior[[paste0("prob", raw_class)]]
}

class_summary <- posterior |>
  group_by(trajectory_group, trajectory_label, raw_class) |>
  summarise(
    n = n(),
    proportion = n() / nrow(posterior),
    APP = mean(assigned_prob),
    .groups = "drop"
  )

mimic_profile <- long |>
  left_join(
    posterior |> select(subject_id_num, trajectory_group, trajectory_label, raw_class),
    by = "subject_id_num"
  ) |>
  group_by(trajectory_group, trajectory_label, raw_class, time0) |>
  summarise(
    n_subjects = n_distinct(subject_id_num),
    mean_sic = mean(sic_score),
    sd_sic = ifelse(n() > 1L, sd(sic_score), NA_real_),
    .groups = "drop"
  ) |>
  mutate(time_day = time0 + 1L)

cdic_profile <- read_csv(cdic_profile_path, show_col_types = FALSE) |>
  mutate(trajectory_group = factor(trajectory_group, levels = paste0("T", 1:5)))
cdic_fitted <- read_csv(cdic_fitted_path, show_col_types = FALSE) |>
  mutate(trajectory_group = factor(trajectory_group, levels = paste0("T", 1:5)))
cdic_assignment <- read_csv(cdic_assignment_path, show_col_types = FALSE) |>
  mutate(trajectory_group = factor(trajectory_group, levels = paste0("T", 1:5)))

comparison <- cdic_profile |>
  select(trajectory_group, time0, CDIC_mean_SIC = mean_sic) |>
  inner_join(
    mimic_profile |> select(trajectory_group, time0, MIMIC_mean_SIC = mean_sic),
    by = c("trajectory_group", "time0")
  ) |>
  mutate(
    Day = time0 + 1L,
    Difference_MIMIC_minus_CDIC = MIMIC_mean_SIC - CDIC_mean_SIC,
    Absolute_difference = abs(Difference_MIMIC_minus_CDIC)
  ) |>
  arrange(trajectory_group, time0)

comparison_publication <- comparison |>
  transmute(
    Trajectory = paste0(
      trajectory_group,
      " ",
      unname(settings$trajectory_names[as.character(trajectory_group)])
    ),
    `Time point` = paste0("D", Day),
    `CDIC mean SIC score` = CDIC_mean_SIC,
    `MIMIC-IV mean SIC score` = MIMIC_mean_SIC,
    `Difference in SIC score (MIMIC-IV minus CDIC)` = Difference_MIMIC_minus_CDIC,
    `Absolute difference` = Absolute_difference
  )

write_csv(class_summary, file.path(settings$out_dir, "MIMIC_classification_summary.csv"), na = "")
write_csv(mimic_profile, file.path(settings$out_dir, "MIMIC_observed_trajectory_profiles.csv"), na = "")
write_csv(comparison_publication, file.path(settings$out_dir, "Table_S8_CDIC_MIMIC_pointwise_profiles.csv"), na = "")
write_csv(posterior, file.path(settings$out_dir, "MIMIC_posterior_assignments.csv"), na = "")

legend_labels <- setNames(
  paste0(names(settings$trajectory_names), ": ", unname(settings$trajectory_names)),
  names(settings$trajectory_names)
)

profile_plot_data <- bind_rows(
  cdic_profile |>
    transmute(trajectory_group, time_day, mean_sic, Cohort = "CDIC"),
  mimic_profile |>
    transmute(trajectory_group, time_day, mean_sic, Cohort = "MIMIC-IV")
)

p_profile <- ggplot(
  profile_plot_data,
  aes(time_day, mean_sic, color = trajectory_group, linetype = Cohort, group = interaction(trajectory_group, Cohort))
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_manual(values = settings$colors, labels = legend_labels) +
  scale_linetype_manual(values = c(CDIC = "solid", `MIMIC-IV` = "3313")) +
  scale_x_continuous(breaks = 1:3, labels = c("D1", "D2", "D3")) +
  scale_y_continuous(breaks = 0:6) +
  coord_cartesian(xlim = c(1, 3), ylim = c(0, 6)) +
  labs(x = NULL, y = "Mean SIC Score", color = NULL, linetype = NULL, title = "Observed SIC trajectories in the derivation and validation cohorts") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right", panel.grid.minor = element_blank())

app_data <- bind_rows(
  cdic_assignment |>
    group_by(trajectory_group) |>
    summarise(APP = mean(assigned_prob), .groups = "drop") |>
    mutate(Cohort = "CDIC"),
  class_summary |>
    transmute(trajectory_group, APP, Cohort = "MIMIC-IV")
)

p_app <- ggplot(app_data, aes(trajectory_group, APP, fill = Cohort)) +
  geom_col(width = 0.66) +
  geom_hline(yintercept = 0.70, linetype = "dashed", color = "grey35") +
  facet_wrap(~Cohort, nrow = 1) +
  scale_fill_manual(values = c(CDIC = "#4C78A8", `MIMIC-IV` = "#F58518"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(x = NULL, y = "Average posterior probability") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

proportion_data <- bind_rows(
  cdic_assignment |>
    count(trajectory_group, name = "n") |>
    mutate(Proportion = n / sum(n), Cohort = "CDIC"),
  class_summary |>
    transmute(trajectory_group, n, Proportion = proportion, Cohort = "MIMIC-IV")
)

p_proportion <- ggplot(proportion_data, aes(trajectory_group, Proportion, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  scale_fill_manual(values = c(CDIC = "#4C78A8", `MIMIC-IV` = "#F58518")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Class proportion", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

supplementary_profile_data <- bind_rows(
  cdic_profile |>
    transmute(
      trajectory_group,
      time_day,
      SIC = mean_sic,
      Panel = "A. CDIC observed class-specific means",
      Source = "Observed"
    ),
  mimic_profile |>
    transmute(
      trajectory_group,
      time_day,
      SIC = mean_sic,
      Panel = "B. MIMIC-IV observed class-specific means",
      Source = "Observed"
    ),
  cdic_fitted |>
    transmute(
      trajectory_group,
      time_day,
      SIC = fitted_sic,
      Panel = "C. Locked LCMM fit and MIMIC-IV observed means",
      Source = "CDIC locked LCMM fit"
    ),
  mimic_profile |>
    transmute(
      trajectory_group,
      time_day,
      SIC = mean_sic,
      Panel = "C. Locked LCMM fit and MIMIC-IV observed means",
      Source = "MIMIC-IV observed"
    )
) |>
  mutate(
    Panel = factor(
      Panel,
      levels = c(
        "A. CDIC observed class-specific means",
        "B. MIMIC-IV observed class-specific means",
        "C. Locked LCMM fit and MIMIC-IV observed means"
      )
    )
  )

p_s4 <- ggplot(
  supplementary_profile_data,
  aes(
    time_day,
    SIC,
    color = trajectory_group,
    linetype = Source,
    group = interaction(trajectory_group, Source)
  )
) +
  geom_line(linewidth = 1.05) +
  geom_point(
    data = supplementary_profile_data |> filter(Source != "CDIC locked LCMM fit"),
    size = 1.8
  ) +
  facet_wrap(~Panel, nrow = 1) +
  scale_color_manual(values = settings$colors, labels = legend_labels) +
  scale_linetype_manual(values = c(Observed = "solid", `CDIC locked LCMM fit` = "solid", `MIMIC-IV observed` = "3313")) +
  scale_x_continuous(breaks = 1:3, labels = c("D1", "D2", "D3")) +
  scale_y_continuous(breaks = 0:6) +
  coord_cartesian(xlim = c(1, 3), ylim = c(0, 6)) +
  labs(x = NULL, y = "Mean SIC Score", color = NULL, linetype = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(settings$out_dir, "figures", "Figure_2B_CDIC_MIMIC_observed_trajectories.pdf"), p_profile, width = 9.2, height = 5.8)
ggsave(file.path(settings$out_dir, "figures", "Figure_S3_average_posterior_probability.pdf"), p_app, width = 7.5, height = 5.3)
ggsave(file.path(settings$out_dir, "figures", "Figure_S4_CDIC_MIMIC_trajectory_profiles.pdf"), p_s4, width = 15, height = 5.6)
ggsave(file.path(settings$out_dir, "figures", "Figure_S5_class_proportions.pdf"), p_proportion, width = 7.5, height = 5.3)
