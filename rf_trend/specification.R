# ---- Model settings ----

# Change the predictor period, repeats, and tuning controls in this list
# Set refit_models to TRUE in the main script after changing these values
default_rf_settings <- function(stability_iterations = 100L) {
  stability_iterations <- as.integer(stability_iterations)
  if (length(stability_iterations) != 1L ||
      !is.finite(stability_iterations) || stability_iterations < 1L) {
    stop("stability_iterations must be one positive integer", call. = FALSE)
  }

  list(
    analysis_version = "si_slopes_v14_full_stability",
    seed = 666L,
    # Annual predictors are summarized over this fixed period
    analysis_start_year = 2002L,
    analysis_end_year = 2022L,
    minimum_driver_observations = 5L,
    # Primary training/testing design
    primary_train_proportion = 0.80,
    primary_strata = "final_cluster",
    outer_repeats = 10L,
    response_bins = 3L,
    # Predictor stability is estimated within each training set only
    inner_stability_iterations = stability_iterations,
    importance_quantile = 0.75,
    stability_frequency_threshold = 0.95,
    fallback_feature_count = 5L,
    # Tree count and mtry are selected automatically
    tree_grid = seq(100L, 2000L, by = 100L),
    mtry_step_factor = 1.5,
    mtry_improvement = 0.01,
    shap_nsim = 20L,
    # Alternative training/testing designs use the same full RF procedure
    sensitivity_train_proportions = c(0.80, 0.60),
    sensitivity_strata = c("final_cluster", "environment_cluster"),
    sensitivity_design_version = "paired_outer_splits_v2",
    sensitivity_repeats = 10L
  )
}

# ---- Input files ----

# Change a filename here only if a shared input file is renamed
rf_input_files <- function(data_root, master_data_root) {
  c(
    concentration = file.path(data_root, "data", "conc_slopes_export.csv"),
    yield = file.path(data_root, "data", "yield_slopes_export.csv"),
    spatial = file.path(
      data_root, "data", "all-data_si-extract_3_20260629.csv"
    ),
    chemistry = file.path(
      master_data_root, "Full_Results_WRTDS_kalman_annual.csv"
    ),
    raw_chemistry = file.path(
      master_data_root, "20260105_masterdata_chem.csv"
    ),
    environment_clusters = file.path(
      data_root, "data", "Si_sites_clusters_six_names.csv"
    )
  )
}

# ---- Predictor definitions ----

# Change the selected predictor lists in this function
# Keep average and trend names aligned with data-prep.R
predictor_specification <- function() {
  land <- c(
    "land_Bare", "land_Cropland", "land_Forest",
    "land_Grassland_Shrubland", "land_Ice_Snow", "land_Impervious",
    "land_Salt_Water", "land_Tidal_Wetland", "land_Water",
    "land_Wetland_Marsh"
  )
  rock <- c(
    "rocks_volcanic", "rocks_sedimentary",
    "rocks_carbonate_evaporite", "rocks_metamorphic", "rocks_plutonic"
  )
  selected_land <- c(
    "land_Cropland", "land_Forest", "land_Grassland_Shrubland",
    "land_Wetland_Marsh"
  )
  selected_rock <- c(
    "rocks_volcanic", "rocks_carbonate_evaporite",
    "rocks_metamorphic", "rocks_plutonic"
  )
  dynamic_average <- c(
    "median_N", "median_P", "npp", "evapotrans", "precip", "temp",
    "snow_cover", "qnorm"
  )
  static <- c(
    "permafrost", "elevation", "basin_slope", "RBI", "recession_slope",
    selected_land, selected_rock
  )
  trend_stubs <- c(
    "N", "P", "npp", "evapotrans", "precip", "temp", "snow_cover",
    "qnorm"
  )
  trends <- paste0("trend_", trend_stubs)

  list(
    land = land,
    rock = rock,
    selected_land = selected_land,
    selected_rock = selected_rock,
    static = static,
    dynamic_average = dynamic_average,
    average = c(static, dynamic_average),
    trends = trends,
    candidates = c(static, dynamic_average, trends)
  )
}

# ---- Response definitions ----

# The main script supplies DSi and Si for the current workflow
slope_model_specification <- function(target_chemical = "DSi",
                                      target_label = "Si") {
  model_stub <- gsub("[^A-Za-z0-9]+", "_", target_label)
  tibble(
    model_id = paste(c("concentration", "yield"), model_stub, sep = "_"),
    data_type = c("concentration", "yield"),
    chemical = target_chemical,
    target_label = target_label,
    plot_label = paste(c("Concentration", "Yield"), "slope")
  )
}
