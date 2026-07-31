# Run this file to fit or load both Si slope models and create the final plots
# Edit the local paths and run switches in the marked sections below
# Edit model settings and predictors in rf_trend/specification.R
# Edit plot labels and appearance in rf_trend/figure2.R

# ---- Packages ----

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "randomForest",
  "fastshap", "cowplot", "data.table"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install_command <- paste0(
    "install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
  stop(
    "Install the missing R packages, then run this script again:\n",
    install_command,
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(randomForest)
  library(fastshap)
  library(cowplot)
})

# ---- Run controls ----

# Use 10 for a quick preliminary run or 100 for the full stability analysis
# Each choice writes to its own output folder so results are never replaced
stability_bootstraps <- 100L
# TRUE refits and replaces saved results for this bootstrap setting
# FALSE reuses matching saved results and only recreates summaries and plots
refit_models <- FALSE
# Leave TRUE to create the three final plots under data_root
make_plots <- TRUE

# ---- Paths ----

script_argument <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
project_dir <- if (length(script_argument) > 0) {
  script_path <- sub("^--file=", "", script_argument[[1]])
  normalizePath(dirname(script_path), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

# Change this to the local long-term-change folder
# Its data subfolder must contain the slopes, spatial data, and clusters
data_root <- paste0(
  "/Users/sidneybush/Library/CloudStorage/Box-Box/",
  "Sidney_Bush/SiSyn/long-term-change"
)
# Change this to the folder containing both annual chemistry input files
master_data_root <- paste0(
  "/Users/sidneybush/Library/CloudStorage/Box-Box/",
  "Sidney_Bush/SiSyn/spatial-data-extractions/master-datasets"
)
# Keep these output paths unchanged to store results outside the repository
run_output_root <- file.path(
  data_root, "rf_outputs", paste0(stability_bootstraps, "_bootstraps")
)
output_dir <- file.path(run_output_root, "model_results")
plot_dir <- file.path(run_output_root, "plots_only")
figure2_path <- file.path(plot_dir, "Figure2_Si_slope_RF_SHAP.png")
sensitivity_plot_path <- file.path(
  plot_dir, "Split_sensitivity_80-20_vs_60-40.png"
)
comparison_plot_path <- file.path(
  plot_dir, "Predictor_set_comparison.png"
)

# ---- Workflow modules ----

module_dir <- file.path(project_dir, "rf_trend")
source(file.path(module_dir, "specification.R"))
source(file.path(module_dir, "data-prep.R"))
source(file.path(module_dir, "resampling.R"))
source(file.path(module_dir, "modeling.R"))
source(file.path(module_dir, "sensitivity.R"))
source(file.path(module_dir, "figure2.R"))

# ---- Analysis settings ----

# Change predictor years, repeats, tuning, or split settings in specification.R
settings <- default_rf_settings(stability_bootstraps)

# ---- Inputs and model definitions ----

input_files <- rf_input_files(data_root, master_data_root)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

predictor_spec <- predictor_specification()
# Keep these values for the current Si-only workflow
target_chemical <- "DSi"
target_label <- "Si"
model_specs <- slope_model_specification(target_chemical, target_label)
predictor_sets <- c(
  "average_only",
  "average_plus_trends"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Model-ready data ----

site_predictors <- build_site_predictors(
  input_files,
  predictor_spec,
  settings
)
prepared_models <- lapply(seq_len(nrow(model_specs)), function(i) {
  spec <- model_specs[i, ]
  slopes <- load_slope_response(
    input_files[[spec$data_type]],
    spec$chemical
  )
  prepare_model_data(
    slopes,
    site_predictors,
    spec,
    predictor_spec
  )
})
names(prepared_models) <- model_specs$model_id

# ---- Primary random forest models ----

comparison_results <- lapply(seq_len(nrow(model_specs)), function(i) {
  spec <- model_specs[i, ]
  setNames(lapply(predictor_sets, function(predictor_set) {
    fit_slope_model(
      prepared_models[[spec$model_id]],
      spec,
      predictor_set,
      settings,
      model_index = i,
      output_dir = output_dir,
      calculate_shap = predictor_set == "average_plus_trends",
      refit_models = refit_models
    )
  }), predictor_sets)
})
names(comparison_results) <- model_specs$model_id
results <- lapply(comparison_results, `[[`, "average_plus_trends")

# ---- Split sensitivity ----

sensitivity_results <- lapply(seq_len(nrow(model_specs)), function(i) {
  spec <- model_specs[i, ]
  run_split_sensitivity(
    prepared_models[[spec$model_id]],
    spec,
    settings,
    model_index = i,
    output_dir = output_dir,
    primary_result = results[[spec$model_id]],
    refit_models = refit_models
  )
})
names(sensitivity_results) <- model_specs$model_id

# ---- Plots and console summary ----

model_summary <- bind_rows(lapply(comparison_results, function(model_sets) {
  bind_rows(lapply(model_sets, `[[`, "summary"))
}))
sensitivity_summary <- bind_rows(lapply(sensitivity_results, `[[`, "summary"))

if (make_plots) {
  create_figure2(results, figure2_path)
  create_predictor_comparison_figure(
    comparison_results,
    comparison_plot_path
  )
  create_split_sensitivity_figure(
    sensitivity_results,
    sensitivity_plot_path
  )
}

print(
  model_summary %>%
    select(
      model_id, predictor_set, model_sites, incomplete_sites_removed,
      significant_sites,
      testing_R2_median, testing_R2_IQR, final_feature_count
    )
)
print(sensitivity_summary)
message("Model results: ", normalizePath(output_dir))
if (make_plots) message("Plots: ", normalizePath(plot_dir))
