# ---- Final all-site model ----

# Fit one documented model for reuse after testing performance is established
fit_final_model <- function(model_data, predictors, settings, model_index) {
  x <- model_data[predictors]
  y <- model_data$response
  assert_complete_predictors(x, "Final model data")
  fitted <- fit_stable_rf(
    x,
    y,
    settings,
    settings$seed + model_index * 50000000L
  )
  importance <- rf_importance(fitted$rf2$model)

  list(
    model = fitted$rf2$model,
    features = fitted$features,
    feature_values = x[, fitted$features, drop = FALSE],
    rf1_tuning = fitted$rf1$settings,
    tuning = fitted$rf2$settings,
    stability = tibble(
      variable = predictors,
      inner_frequency = as.numeric(fitted$stability$frequency[predictors]),
      baseline_importance = as.numeric(
        fitted$stability$baseline_importance[predictors]
      ),
      selected = predictors %in% fitted$features
    ),
    feature_source = fitted$stability$feature_source,
    importance = tibble(
      variable = names(importance),
      permutation_importance = as.numeric(importance)
    ) %>%
      arrange(desc(permutation_importance)),
    oob_r2 = tail(fitted$rf2$model$rsq, 1)
  )
}

# ---- Model summary ----

summarize_model_result <- function(spec, predictor_set, audit, validation,
                                   final_model, settings) {
  metrics <- validation$outer_metrics
  strata_label <- if (settings$primary_strata == "final_cluster") {
    "GRL lithology"
  } else if (settings$primary_strata == "environment_cluster") {
    "environmental clusters"
  } else {
    settings$primary_strata
  }

  audit %>%
    mutate(
      analysis_version = settings$analysis_version,
      predictor_set = predictor_set,
      validation_method = paste(
        settings$outer_repeats,
        "repeated site splits with",
        round(100 * settings$primary_train_proportion),
        "percent training,", strata_label, "and response balancing,",
        "and training-only bootstrap stability selection"
      ),
      testing_R2_mean = mean(metrics$R2),
      testing_R2_median = median(metrics$R2),
      testing_R2_sd = sd(metrics$R2),
      testing_R2_IQR = IQR(metrics$R2),
      testing_Pearson_R2_median = median(metrics$Pearson_R2),
      testing_RMSE_mean = mean(metrics$RMSE),
      testing_RMSE_median = median(metrics$RMSE),
      testing_RMSE_IQR = IQR(metrics$RMSE),
      testing_site_coverage = validation$testing_site_coverage,
      median_stable_features = median(metrics$stable_features),
      final_feature_count = length(final_model$features),
      final_features = paste(final_model$features, collapse = "; "),
      final_oob_R2 = final_model$oob_r2,
      final_num_trees = final_model$tuning$num.trees,
      final_mtry = final_model$tuning$mtry
    )
}

# ---- Cached model run ----

fit_slope_model <- function(prepared, spec, predictor_set, settings,
                            model_index, output_dir, calculate_shap,
                            refit_models = FALSE) {
  model_dir <- file.path(output_dir, spec$model_id, predictor_set)
  result_path <- file.path(model_dir, "model_result.rds")

  if (file.exists(result_path) && !refit_models) {
    result <- readRDS(result_path)
    if (identical(result$analysis_version, settings$analysis_version)) {
      message("Using saved model: ", spec$model_id, " / ", predictor_set)
      return(result)
    }
  }

  message("Fitting model: ", spec$model_id, " / ", predictor_set)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  predictors <- prepared$predictor_sets[[predictor_set]]
  validation <- run_repeated_validation(
    prepared$data,
    predictors,
    settings,
    model_index,
    train_proportion = settings$primary_train_proportion,
    strata_variable = settings$primary_strata,
    repeats = settings$outer_repeats,
    calculate_shap = calculate_shap
  )
  message("  Fitting final all-site model")
  final_model <- fit_final_model(
    prepared$data,
    predictors,
    settings,
    model_index
  )
  summary <- summarize_model_result(
    spec,
    predictor_set,
    prepared$audit,
    validation,
    final_model,
    settings
  )

  result <- list(
    analysis_version = settings$analysis_version,
    specification = spec,
    predictor_set = predictor_set,
    settings = settings,
    audit = prepared$audit,
    missingness = prepared$missingness,
    excluded_sites = prepared$excluded_sites,
    trend_coverage = prepared$trend_coverage,
    candidate_predictors = predictors,
    site_data = prepared$data %>%
      select(
        stream_key, Stream_Name, response, p_value, significant,
        final_cluster, environment_cluster, environment_cluster_name,
        predictor_start_year, predictor_end_year, predictor_years
      ),
    validation = validation,
    final = final_model,
    summary = summary
  )
  saveRDS(result, result_path, compress = TRUE)
  result
}
