#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

settings <- list(
  input_file = file.path("data", "CDIC_D1D3_SIC_lcmm.xlsx"),
  input_sheet = "Sheet1",
  out_dir = file.path("results", "01_CDIC_trajectory_modeling"),
  id_var = "patient_id_cleaned",
  sic_vars = c("sic_d1", "sic_d2", "sic_d3"),
  missing_tokens = c("", "NA", "N/A", ".", "NULL", "NaN"),
  missing_codes = c("20", "35"),
  min_observed = 2L,
  candidate_ng = 1:6,
  selected_ng = 5L,
  seed = 20260526L,
  gridsearch_reps = 120L,
  gridsearch_maxiter = 50L,
  hlme_maxiter = 1000L,
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
  if (key %in% c("input_file", "input_sheet", "out_dir")) settings[[key]] <- value
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

set.seed(settings$seed)
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
safe_min <- function(x) if (length(x) == 0L || all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)

clean_sic <- function(x, variable) {
  z <- trimws(as.character(x))
  missing <- is.na(z) | toupper(z) %in% toupper(settings$missing_tokens) | z %in% settings$missing_codes
  value <- as_num(z)
  invalid_text <- !missing & is.na(value)
  if (any(invalid_text)) {
    stop("Unexpected nonnumeric value in ", variable, ": ", paste(head(unique(z[invalid_text]), 10), collapse = ", "))
  }
  valid <- !missing & is.finite(value) & value >= 0 & value <= 6 &
    abs(value - round(value)) <= sqrt(.Machine$double.eps)
  value[!valid] <- NA_real_
  value
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

probability_columns <- function(x) {
  out <- grep("^prob[0-9]+$", names(x), value = TRUE)
  if (!length(out)) out <- grep("^prob", names(x), value = TRUE)
  out
}

repair_call <- function(fit, ng) {
  fit$call$fixed <- as.formula("sic_score ~ time0", env = .GlobalEnv)
  fit$call$random <- as.formula("~ 1", env = .GlobalEnv)
  if (ng > 1L) {
    fit$call$mixture <- as.formula("~ time0", env = .GlobalEnv)
    fit$call$classmb <- as.formula("~ 1", env = .GlobalEnv)
  }
  fit
}

valid_fit <- function(fit, ng, n_subjects) {
  if (is.null(fit) || is.null(fit$conv) || fit$conv != 1L) return(FALSE)
  if (any(!is.finite(c(fit$loglik, fit$AIC, fit$BIC, fit$best)))) return(FALSE)
  if (ng == 1L) return(TRUE)
  pp <- as.data.frame(fit$pprob)
  pc <- probability_columns(pp)
  if (length(pc) < ng) return(FALSE)
  p <- as.matrix(data.frame(lapply(pp[, pc[seq_len(ng)], drop = FALSE], as_num)))
  nrow(p) == n_subjects && all(is.finite(p)) &&
    all(p >= -1e-8 & p <= 1 + 1e-8) &&
    max(abs(rowSums(p) - 1)) <= 1e-4
}

extract_posterior <- function(fit, ng) {
  if (ng == 1L) {
    return(tibble(subject_id_num = wide$subject_id_num, raw_class = 1L, prob1 = 1, max_postprob = 1))
  }
  pp <- as.data.frame(fit$pprob)
  pc <- probability_columns(pp)[seq_len(ng)]
  p <- as.data.frame(lapply(pp[, pc, drop = FALSE], as_num))
  names(p) <- paste0("prob", seq_len(ng))
  subject <- if ("subject_id_num" %in% names(pp)) pp$subject_id_num else if ("subject" %in% names(pp)) pp$subject else pp[[1]]
  cls <- max.col(as.matrix(p), ties.method = "first")
  bind_cols(tibble(subject_id_num = as.integer(subject), raw_class = cls), p) |>
    mutate(
      max_postprob = apply(as.matrix(p), 1, max),
      assigned_prob = as.matrix(p)[cbind(seq_len(nrow(p)), cls)]
    )
}

class_metrics <- function(pp, ng) {
  bind_rows(lapply(seq_len(ng), function(k) {
    selected <- pp$raw_class == k
    tibble(
      raw_class = k,
      class_n = sum(selected),
      class_prop = mean(selected),
      APP = mean(pp[[paste0("prob", k)]][selected]),
      mean_max_postprob = mean(pp$max_postprob[selected])
    )
  }))
}

entropy_value <- function(pp, ng) {
  if (ng == 1L) return(1)
  p <- as.matrix(pp[, paste0("prob", seq_len(ng)), drop = FALSE])
  p <- pmax(pmin(p, 1), 1e-12)
  1 - sum(-p * log(p)) / (nrow(p) * log(ng))
}

fit_ng1 <- hlme(
  fixed = sic_score ~ time0,
  random = ~ 1,
  subject = "subject_id_num",
  ng = 1,
  data = long,
  na.action = 1,
  maxiter = settings$hlme_maxiter,
  var.time = "time0",
  verbose = FALSE
)
fit_ng1 <- repair_call(fit_ng1, 1L)
if (!valid_fit(fit_ng1, 1L, nrow(wide))) stop("The one-class model did not converge.")

fits <- list(`1` = fit_ng1)
candidate_rows <- list()

for (ng in settings$candidate_ng) {
  fit <- if (ng == 1L) {
    fit_ng1
  } else {
    grid_fit <- tryCatch(
      gridsearch(
        m = hlme(
          fixed = sic_score ~ time0,
          mixture = ~ time0,
          random = ~ 1,
          subject = "subject_id_num",
          ng = ng,
          nwg = FALSE,
          data = long,
          na.action = 1,
          maxiter = settings$hlme_maxiter,
          var.time = "time0",
          verbose = FALSE
        ),
        rep = settings$gridsearch_reps,
        maxiter = settings$gridsearch_maxiter,
        minit = fit_ng1
      ),
      error = function(e) NULL
    )
    grid_fit <- if (is.null(grid_fit)) NULL else repair_call(grid_fit, ng)
    if (!valid_fit(grid_fit, ng, nrow(wide))) {
      direct_fit <- tryCatch(
        hlme(
          fixed = sic_score ~ time0,
          mixture = ~ time0,
          random = ~ 1,
          subject = "subject_id_num",
          ng = ng,
          nwg = FALSE,
          data = long,
          B = fit_ng1,
          na.action = 1,
          maxiter = settings$hlme_maxiter,
          var.time = "time0",
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      direct_fit <- if (is.null(direct_fit)) NULL else repair_call(direct_fit, ng)
      candidates <- Filter(function(x) valid_fit(x, ng, nrow(wide)), list(grid_fit, direct_fit))
      if (!length(candidates)) NULL else candidates[[which.min(vapply(candidates, function(x) x$BIC, numeric(1)))]]
    } else {
      grid_fit
    }
  }
  if (is.null(fit) || !valid_fit(fit, ng, nrow(wide))) stop("Candidate ng=", ng, " failed strict validation.")
  fits[[as.character(ng)]] <- fit
  pp <- extract_posterior(fit, ng)
  cm <- class_metrics(pp, ng)
  candidate_rows[[as.character(ng)]] <- tibble(
    Number_of_classes = ng,
    Log_likelihood = fit$loglik,
    AIC = fit$AIC,
    BIC = fit$BIC,
    Entropy = entropy_value(pp, ng),
    Minimum_APP = safe_min(cm$APP),
    Mean_APP = mean(cm$APP),
    Smallest_class_percent = 100 * safe_min(cm$class_prop),
    Eligible = safe_min(cm$APP) >= 0.70 && safe_min(cm$class_prop) >= 0.05,
    Decision = ifelse(ng == settings$selected_ng, "Selected", "")
  )
}

candidate_table <- bind_rows(candidate_rows)
candidate_table <- candidate_table |>
  transmute(
    `Number of classes` = Number_of_classes,
    `Log likelihood` = Log_likelihood,
    AIC,
    BIC,
    Entropy,
    `Minimum APP` = Minimum_APP,
    `Mean APP` = Mean_APP,
    `Smallest class size` = sprintf("%.1f%%", Smallest_class_percent),
    Eligible = ifelse(Eligible, "Yes", "No"),
    Decision
  )
write_csv(candidate_table, file.path(settings$out_dir, "Table_S4_candidate_LCMM_models.csv"), na = "")

fit <- fits[[as.character(settings$selected_ng)]]
posterior <- extract_posterior(fit, settings$selected_ng) |>
  left_join(wide |> select(subject_id_num, all_of(settings$id_var)), by = "subject_id_num") |>
  mutate(
    trajectory_group = unname(settings$raw_to_trajectory[as.character(raw_class)]),
    trajectory_group = factor(trajectory_group, levels = paste0("T", 1:5))
  )

for (group in paste0("T", 1:5)) {
  raw_class <- as.integer(names(settings$raw_to_trajectory)[settings$raw_to_trajectory == group])
  posterior[[paste0("prob_", group)]] <- posterior[[paste0("prob", raw_class)]]
}

class_summary <- posterior |>
  group_by(trajectory_group, raw_class) |>
  summarise(
    n = n(),
    proportion = n() / nrow(posterior),
    APP = mean(assigned_prob),
    .groups = "drop"
  ) |>
  mutate(trajectory_label = unname(settings$trajectory_names[as.character(trajectory_group)]))

observed_profile <- long |>
  left_join(posterior |> select(subject_id_num, trajectory_group), by = "subject_id_num") |>
  group_by(trajectory_group, time0) |>
  summarise(
    n_subjects = n_distinct(subject_id_num),
    mean_sic = mean(sic_score),
    sd_sic = ifelse(n() > 1L, sd(sic_score), NA_real_),
    .groups = "drop"
  ) |>
  mutate(time_day = time0 + 1L)

fit <- repair_call(fit, settings$selected_ng)
grid <- data.frame(time0 = seq(0, 2, length.out = 201))
predicted <- as.matrix(predictY(fit, newdata = grid, var.time = "time0", draws = FALSE)$pred)
fitted_profile <- bind_rows(lapply(seq_len(settings$selected_ng), function(k) {
  tibble(
    raw_class = k,
    time0 = grid$time0,
    time_day = grid$time0 + 1,
    fitted_sic = predicted[, k],
    trajectory_group = unname(settings$raw_to_trajectory[as.character(k)])
  )
})) |>
  mutate(trajectory_group = factor(trajectory_group, levels = paste0("T", 1:5))) |>
  arrange(trajectory_group, time0)

best <- fit$best
covariance <- matrix(0, length(best), length(best))
covariance[upper.tri(covariance, diag = TRUE)] <- fit$V
covariance <- covariance + t(covariance) - diag(diag(covariance))
standard_error <- sqrt(diag(covariance))
intercept_index <- tail(which(grepl("^intercept class", names(best))), settings$selected_ng)
slope_index <- which(grepl("^time0 class", names(best)))

fixed_effects <- bind_rows(lapply(seq_len(settings$selected_ng), function(k) {
  group <- unname(settings$raw_to_trajectory[as.character(k)])
  bind_rows(
    tibble(Model_term = "Intercept", Trajectory = group, Estimate = best[intercept_index[k]], Standard_error = standard_error[intercept_index[k]]),
    tibble(Model_term = "time0", Trajectory = group, Estimate = best[slope_index[k]], Standard_error = standard_error[slope_index[k]])
  )
})) |>
  mutate(
    z_statistic = Estimate / Standard_error,
    P_value = 2 * pnorm(-abs(z_statistic)),
    Trajectory = factor(Trajectory, levels = paste0("T", 1:5))
  ) |>
  arrange(Trajectory, factor(Model_term, levels = c("Intercept", "time0"))) |>
  transmute(
    `Model term` = Model_term,
    `Trajectory label` = paste0(
      Trajectory,
      " ",
      unname(settings$trajectory_names[as.character(Trajectory)])
    ),
    Coefficient = Estimate,
    `Standard error` = Standard_error,
    `Wald z statistic` = z_statistic,
    `P-value` = P_value
  )

write_csv(fixed_effects, file.path(settings$out_dir, "Table_S5_selected_ng5_fixed_effects.csv"), na = "")
write_csv(class_summary, file.path(settings$out_dir, "CDIC_classification_summary.csv"), na = "")
write_csv(observed_profile, file.path(settings$out_dir, "CDIC_observed_trajectory_profiles.csv"), na = "")
write_csv(fitted_profile, file.path(settings$out_dir, "CDIC_fitted_trajectory_profiles.csv"), na = "")
write_csv(posterior, file.path(settings$out_dir, "CDIC_posterior_assignments.csv"), na = "")
saveRDS(fit, file.path(settings$out_dir, "CDIC_locked_ng5_model.rds"))

legend_labels <- setNames(
  paste0(names(settings$trajectory_names), ": ", unname(settings$trajectory_names)),
  names(settings$trajectory_names)
)

p_trajectory <- ggplot(
  fitted_profile,
  aes(time_day, fitted_sic, color = trajectory_group, group = trajectory_group)
) +
  geom_line(linewidth = 1.2, lineend = "round") +
  scale_color_manual(values = settings$colors, labels = legend_labels) +
  scale_x_continuous(breaks = 1:3, labels = c("D1", "D2", "D3")) +
  scale_y_continuous(breaks = 2:6) +
  coord_cartesian(xlim = c(1, 3), ylim = c(2, 6)) +
  labs(x = NULL, y = "SIC Score", color = NULL, title = "Five SIC trajectory groups", subtitle = "LCMM-fitted trajectories") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "right", panel.grid.minor = element_blank())

ggsave(file.path(settings$out_dir, "figures", "Figure_2A_CDIC_fitted_trajectories.pdf"), p_trajectory, width = 8.5, height = 5.6)
