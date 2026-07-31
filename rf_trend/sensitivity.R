# ---- Split-design sensitivity ----

# Apply the same tuning and stability pipeline to all four split designs
run_split_sensitivity <- function(prepared, spec, settings, model_index,
                                  output_dir, primary_result,
                                  refit_models = FALSE) {
  predictor_set <- "average_plus_trends"
  sensitivity_dir <- file.path(output_dir, "split-sensitivity")
  result_path <- file.path(
    sensitivity_dir,
    paste0(spec$model_id, "_split_sensitivity.rds")
  )

  if (file.exists(result_path) && !refit_models) {
    result <- readRDS(result_path)
    if (
      identical(result$analysis_version, settings$analysis_version) &&
      identical(result$stratification_methods, settings$sensitivity_strata) &&
      identical(
        result$sensitivity_design_version,
        settings$sensitivity_design_version
      )
    ) {
      message("Using saved split sensitivity: ", spec$model_id)
      return(result)
    }
  }

  dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)
  model_data <- prepared$data
  predictors <- prepared$predictor_sets[[predictor_set]]
  metrics_list <- list()
  feature_list <- list()
  counter <- 0L
  message("Running split sensitivity: ", spec$model_id)

  for (strata_variable in settings$sensitivity_strata) {
    stratification <- if (
      identical(strata_variable, "final_cluster")
    ) {
      "Dominant lithology"
    } else {
      "Environmental clusters"
    }

    for (train_proportion in settings$sensitivity_train_proportions) {
      counter <- counter + 1L
      split_label <- sprintf(
        "%d/%d",
        round(100 * train_proportion),
        round(100 * (1 - train_proportion))
      )
      is_primary <- identical(strata_variable, settings$primary_strata) &&
        isTRUE(all.equal(train_proportion, settings$primary_train_proportion))

      if (is_primary) {
        message("  Reusing primary ", split_label, " ", stratification)
        validation <- primary_result$validation
      } else {
        message("  Fitting ", split_label, " ", stratification)
        validation <- run_repeated_validation(
          model_data,
          predictors,
          settings,
          model_index,
          train_proportion = train_proportion,
          strata_variable = strata_variable,
          repeats = settings$sensitivity_repeats,
          calculate_shap = FALSE
        )
      }

      metrics_list[[counter]] <- validation$outer_metrics %>%
        mutate(
          model_id = spec$model_id,
          predictor_set = predictor_set,
          stratification = stratification,
          strata_variable = strata_variable,
          split = split_label,
          .before = 1
        )
      feature_list[[counter]] <- validation$feature_summary %>%
        mutate(
          model_id = spec$model_id,
          stratification = stratification,
          strata_variable = strata_variable,
          split = split_label,
          .before = 1
        )
    }
  }

  metrics <- bind_rows(metrics_list)
  result <- list(
    analysis_version = settings$analysis_version,
    specification = spec,
    predictor_set = predictor_set,
    stratification_methods = settings$sensitivity_strata,
    sensitivity_design_version = settings$sensitivity_design_version,
    metrics = metrics,
    feature_stability = bind_rows(feature_list),
    summary = metrics %>%
      group_by(
        model_id, stratification, strata_variable, split, train_proportion
      ) %>%
      summarise(
        repeats = n(),
        train_sites_median = median(training_sites),
        test_sites_median = median(testing_sites),
        R2_median = median(R2),
        R2_IQR = IQR(R2),
        R2_mean = mean(R2),
        RMSE_median = median(RMSE),
        RMSE_IQR = IQR(RMSE),
        RMSE_mean = mean(RMSE),
        stable_features_median = median(stable_features),
        .groups = "drop"
      )
  )
  saveRDS(result, result_path, compress = TRUE)
  result
}
