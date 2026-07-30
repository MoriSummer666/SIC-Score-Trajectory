#!/usr/bin/env Rscript

options(
  sic.survival.settings = list(
    cohort = "MIMIC",
    input_file = file.path(
      "results",
      "03_random_forest_imputation",
      "MIMIC_D1D3_covariate_missForest_completed_dataset_survival_used.xlsx"
    ),
    input_sheet = 1,
    posterior_file = file.path("results", "02_MIMIC_external_validation", "MIMIC_posterior_assignments.csv"),
    classification_file = file.path("results", "02_MIMIC_external_validation", "MIMIC_classification_summary.csv"),
    profile_file = file.path("results", "02_MIMIC_external_validation", "MIMIC_observed_trajectory_profiles.csv"),
    out_dir = file.path("results", "05_MIMIC_survival_analysis")
  )
)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_directory <- dirname(normalizePath(sub("^--file=", "", script_argument)))
source(file.path(script_directory, "04_CDIC_survival_analysis.R"), chdir = FALSE)
