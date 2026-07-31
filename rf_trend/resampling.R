# ---- Metrics ----

regression_metrics <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]

  if (length(observed) < 3) {
    return(tibble(
      n = length(observed), R2 = NA_real_, RMSE = NA_real_,
      MAE = NA_real_, Pearson_R2 = NA_real_
    ))
  }

  sst <- sum((observed - mean(observed))^2)
  tibble(
    n = length(observed),
    R2 = if (sst > 0) 1 - sum((observed - predicted)^2) / sst else NA_real_,
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted)),
    Pearson_R2 = suppressWarnings(cor(observed, predicted)^2)
  )
}

# ---- Complete-case check ----

assert_complete_predictors <- function(x, context) {
  incomplete <- names(x)[vapply(x, function(column) {
    any(!is.finite(column))
  }, logical(1))]
  if (length(incomplete) > 0) {
    stop(
      context, " contains missing predictors: ",
      paste(incomplete, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(x)
}

# ---- Automatic RF tuning ----

# Tree count and mtry are selected from training data only
tune_rf <- function(x, y, settings, seed_value) {
  assert_complete_predictors(x, "RF tuning data")
  tree_grid <- sort(unique(as.integer(settings$tree_grid)))
  if (length(tree_grid) == 0 || any(tree_grid < 1)) {
    stop("tree_grid must contain positive integers", call. = FALSE)
  }

  set.seed(seed_value)
  tree_scan <- randomForest::randomForest(
    x = x,
    y = y,
    ntree = max(tree_grid),
    importance = TRUE
  )
  tree_mse <- tree_scan$mse[tree_grid]
  selected_trees <- tree_grid[[which.min(tree_mse)]]

  if (ncol(x) == 1) {
    selected_mtry <- 1L
    mtry_grid <- tibble(mtry = 1L, oob_error = NA_real_)
  } else {
    set.seed(seed_value + 1L)
    mtry_scan <- NULL
    mtry_start <- max(1L, floor(ncol(x) / 3))
    invisible(capture.output({
      mtry_scan <- randomForest::tuneRF(
        x = x,
        y = y,
        mtryStart = mtry_start,
        ntreeTry = selected_trees,
        stepFactor = settings$mtry_step_factor,
        improve = settings$mtry_improvement,
        trace = FALSE,
        plot = FALSE,
        doBest = FALSE
      )
    }))
    selected_mtry <- as.integer(mtry_scan[which.min(mtry_scan[, 2]), 1])
    selected_mtry <- max(1L, min(ncol(x), selected_mtry))
    mtry_grid <- as.data.frame(mtry_scan) %>%
      transmute(
        mtry = as.integer(mtry),
        oob_error = as.numeric(OOBError)
      )
  }

  set.seed(seed_value + 2L)
  model <- randomForest::randomForest(
    x = x,
    y = y,
    ntree = selected_trees,
    mtry = selected_mtry,
    importance = TRUE
  )

  list(
    model = model,
    settings = tibble(
      num.trees = selected_trees,
      mtry = selected_mtry,
      oob_mse = tail(model$mse, 1)
    ),
    tree_grid = tibble(num.trees = tree_grid, oob_mse = tree_mse),
    mtry_grid = mtry_grid
  )
}

rf_importance <- function(model) {
  values <- randomForest::importance(model, type = 1, scale = TRUE)
  if (is.matrix(values)) {
    output <- as.numeric(values[, 1])
    names(output) <- rownames(values)
  } else {
    output <- as.numeric(values)
    names(output) <- names(values)
  }
  output[!is.finite(output)] <- 0
  output
}

# ---- Stratified site splits ----

# Build one balanced site order so paired split sizes use the same sites first
make_stratified_order <- function(model_data, strata_variable, response_bins,
                                  seed_value) {
  if (anyDuplicated(model_data$stream_key)) {
    stop("Each site must appear once before splitting", call. = FALSE)
  }
  group <- model_data[[strata_variable]]
  group_rows <- c(
    list(which(is.na(group))),
    lapply(sort(unique(group[!is.na(group)])), function(value) {
      which(!is.na(group) & group == value)
    })
  )
  group_rows <- Filter(length, group_rows)
  response_bin <- integer(nrow(model_data))
  cells <- list()

  for (rows in group_rows) {
    ordered <- rows[order(model_data$response[rows], model_data$stream_key[rows])]
    response_bin[ordered] <- pmin(
      response_bins,
      floor((seq_along(ordered) - 1L) * response_bins / length(ordered)) + 1L
    )
    cells <- c(cells, unname(split(rows, response_bin[rows])))
  }

  set.seed(seed_value)
  cells <- lapply(cells, function(rows) rows[sample.int(length(rows))])
  cell_size <- lengths(cells)
  selected <- integer(length(cells))
  next_site <- rep(1L, length(cells))
  tie_order <- sample.int(length(cells))
  site_order <- integer(nrow(model_data))

  for (position in seq_len(nrow(model_data))) {
    deficit <- position * cell_size / nrow(model_data) - selected
    deficit[selected == cell_size] <- -Inf
    candidates <- which(abs(deficit - max(deficit)) < 1e-12)
    chosen <- tie_order[which(tie_order %in% candidates)[1]]
    site_order[[position]] <- cells[[chosen]][next_site[[chosen]]]
    next_site[[chosen]] <- next_site[[chosen]] + 1L
    selected[[chosen]] <- selected[[chosen]] + 1L
  }

  list(site_order = site_order, response_bin = response_bin)
}

# Select an exact testing prefix from the shared balanced site order
make_site_split <- function(model_data, train_proportion, seed_value,
                            strata_variable, response_bins) {
  n_sites <- nrow(model_data)
  n_train <- round(n_sites * train_proportion)
  n_train <- max(1L, min(n_sites - 1L, n_train))
  n_test <- n_sites - n_train
  ordered <- make_stratified_order(
    model_data, strata_variable, response_bins, seed_value
  )
  testing_rows <- ordered$site_order[seq_len(n_test)]
  split_group <- as.character(model_data[[strata_variable]])
  stratum_group <- ifelse(is.na(split_group), "Unassigned", split_group)

  list(
    train = setdiff(seq_len(n_sites), testing_rows),
    test = sort(testing_rows),
    assignment = tibble(
      row_id = seq_len(n_sites),
      split_group = split_group,
      response_bin = ordered$response_bin,
      stratum = paste(stratum_group, ordered$response_bin, sep = "__"),
      subset = if_else(
        seq_len(n_sites) %in% testing_rows,
        "testing",
        "training"
      )
    )
  )
}

# ---- Stability selection ----

# Bootstrap only the outer training sites and retain repeatable predictors
bootstrap_stability <- function(x, y, rf1, settings, seed_value) {
  baseline_importance <- rf_importance(rf1)
  importance_cutoff <- as.numeric(quantile(
    baseline_importance,
    settings$importance_quantile,
    na.rm = TRUE
  ))
  iterations <- as.integer(settings$inner_stability_iterations)
  selected <- matrix(
    FALSE,
    nrow = iterations,
    ncol = ncol(x),
    dimnames = list(NULL, colnames(x))
  )

  for (i in seq_len(iterations)) {
    set.seed(seed_value + i)
    sampled_rows <- sample(seq_len(nrow(x)), nrow(x), replace = TRUE)
    model <- randomForest::randomForest(
      x = x[sampled_rows, , drop = FALSE],
      y = y[sampled_rows],
      ntree = rf1$ntree,
      mtry = rf1$mtry,
      importance = TRUE
    )
    importance <- rf_importance(model)[colnames(x)]
    selected[i, ] <- importance > importance_cutoff
  }

  frequency <- colMeans(selected)
  stable_features <- names(frequency)[
    frequency >= settings$stability_frequency_threshold
  ]
  feature_source <- "stability_selection"
  if (length(stable_features) == 0) {
    stable_features <- names(sort(
      baseline_importance,
      decreasing = TRUE
    ))[seq_len(min(settings$fallback_feature_count, ncol(x)))]
    feature_source <- "importance_fallback"
  }

  list(
    features = stable_features,
    frequency = frequency,
    baseline_importance = baseline_importance,
    importance_cutoff = importance_cutoff,
    feature_source = feature_source
  )
}

# Tune RF1, select stable predictors, then tune and fit RF2
fit_stable_rf <- function(x, y, settings, seed_value) {
  rf1 <- tune_rf(x, y, settings, seed_value)
  stability <- bootstrap_stability(
    x,
    y,
    rf1$model,
    settings,
    seed_value + 100000L
  )
  selected_x <- x[, stability$features, drop = FALSE]
  rf2 <- tune_rf(selected_x, y, settings, seed_value + 200000L)

  list(
    rf1 = rf1,
    stability = stability,
    rf2 = rf2,
    features = stability$features
  )
}

# ---- Testing SHAP ----

predict_rf <- function(object, newdata) {
  as.numeric(predict(object, newdata = as.data.frame(newdata)))
}

compute_testing_shap <- function(model, training_data, testing_data, nsim,
                                 seed_value) {
  set.seed(seed_value)
  as.data.frame(fastshap::explain(
    object = model,
    X = as.data.frame(training_data),
    newdata = as.data.frame(testing_data),
    pred_wrapper = predict_rf,
    nsim = nsim
  ))
}

# ---- Repeated outer validation ----

# Testing sites validate each fitted model and never enter tuning or selection

run_repeated_validation <- function(model_data, predictors, settings,
                                    model_index, train_proportion,
                                    strata_variable, repeats,
                                    calculate_shap = TRUE) {
  outer_metrics <- vector("list", repeats)
  prediction_list <- vector("list", repeats)
  assignment_list <- vector("list", repeats)
  feature_list <- vector("list", repeats)
  shap_list <- if (calculate_shap) vector("list", repeats) else NULL

  for (outer_split in seq_len(repeats)) {
    message("  Outer split ", outer_split, " of ", repeats)
    split_seed <- settings$seed + model_index * 1000000L + outer_split * 1000L
    partition <- make_site_split(
      model_data,
      train_proportion,
      split_seed,
      strata_variable,
      settings$response_bins
    )
    x_train <- model_data[partition$train, predictors, drop = FALSE]
    x_test <- model_data[partition$test, predictors, drop = FALSE]
    y_train <- model_data$response[partition$train]
    y_test <- model_data$response[partition$test]
    assert_complete_predictors(x_train, "Training data")
    assert_complete_predictors(x_test, "Testing data")

    fitted <- fit_stable_rf(
      x_train,
      y_train,
      settings,
      split_seed + 300000L
    )
    features <- fitted$features
    testing_prediction <- predict_rf(
      fitted$rf2$model,
      x_test[, features, drop = FALSE]
    )
    training_prediction <- predict_rf(
      fitted$rf2$model,
      x_train[, features, drop = FALSE]
    )
    testing_metrics <- regression_metrics(y_test, testing_prediction)
    training_metrics <- regression_metrics(y_train, training_prediction)

    outer_metrics[[outer_split]] <- testing_metrics %>%
      mutate(
        outer_split = outer_split,
        train_proportion = train_proportion,
        training_sites = length(partition$train),
        testing_sites = length(partition$test),
        training_R2 = training_metrics$R2,
        training_RMSE = training_metrics$RMSE,
        rf1_ntree = fitted$rf1$settings$num.trees,
        rf1_mtry = fitted$rf1$settings$mtry,
        importance_cutoff = fitted$stability$importance_cutoff,
        stable_features = length(features),
        feature_source = fitted$stability$feature_source,
        rf2_ntree = fitted$rf2$settings$num.trees,
        rf2_mtry = fitted$rf2$settings$mtry,
        .before = 1
      )
    prediction_list[[outer_split]] <- tibble(
      outer_split = outer_split,
      stream_key = model_data$stream_key[partition$test],
      Stream_Name = model_data$Stream_Name[partition$test],
      observed = y_test,
      predicted = testing_prediction
    )
    assignment_list[[outer_split]] <- partition$assignment %>%
      mutate(
        outer_split = outer_split,
        stream_key = model_data$stream_key,
        Stream_Name = model_data$Stream_Name,
        .before = 1
      )
    feature_list[[outer_split]] <- tibble(
      outer_split = outer_split,
      variable = predictors,
      inner_frequency = as.numeric(
        fitted$stability$frequency[predictors]
      ),
      baseline_importance = as.numeric(
        fitted$stability$baseline_importance[predictors]
      ),
      selected = predictors %in% features
    )

    if (!calculate_shap) next
    shap_values <- compute_testing_shap(
      fitted$rf2$model,
      x_train[, features, drop = FALSE],
      x_test[, features, drop = FALSE],
      settings$shap_nsim,
      split_seed + 400000L
    )
    identifiers <- tibble(
      outer_split = outer_split,
      stream_key = model_data$stream_key[partition$test],
      Stream_Name = model_data$Stream_Name[partition$test]
    )
    shap_list[[outer_split]] <- bind_cols(identifiers, shap_values) %>%
      pivot_longer(
        cols = all_of(features),
        names_to = "variable",
        values_to = "shap"
      ) %>%
      left_join(
        bind_cols(identifiers, x_test[, features, drop = FALSE]) %>%
          pivot_longer(
            cols = all_of(features),
            names_to = "variable",
            values_to = "predictor_value"
          ),
        by = c("outer_split", "stream_key", "Stream_Name", "variable")
      )
  }

  predictions <- bind_rows(prediction_list)
  outer_metrics <- bind_rows(outer_metrics)
  feature_detail <- bind_rows(feature_list)
  feature_summary <- feature_detail %>%
    group_by(variable) %>%
    summarise(
      outer_selection_frequency = mean(selected),
      median_inner_frequency = median(inner_frequency),
      median_baseline_importance = median(baseline_importance),
      .groups = "drop"
    ) %>%
    arrange(desc(outer_selection_frequency), desc(median_inner_frequency))
  site_predictions <- predictions %>%
    group_by(stream_key, Stream_Name, observed) %>%
    summarise(
      testing_appearances = n(),
      predicted = mean(predicted),
      predicted_sd = if_else(n() > 1, sd(predicted), NA_real_),
      .groups = "drop"
    )

  output <- list(
    outer_metrics = outer_metrics,
    predictions = predictions,
    assignments = bind_rows(assignment_list),
    feature_detail = feature_detail,
    feature_summary = feature_summary,
    site_predictions = site_predictions,
    site_metrics = regression_metrics(
      site_predictions$observed,
      site_predictions$predicted
    ),
    testing_site_coverage = nrow(site_predictions) / nrow(model_data)
  )
  if (!calculate_shap) return(output)

  testing_shap <- bind_rows(shap_list)
  output$testing_shap <- testing_shap
  output$testing_shap_site <- testing_shap %>%
    group_by(stream_key, Stream_Name, variable) %>%
    summarise(
      shap = mean(shap),
      predictor_value = mean(predictor_value),
      .groups = "drop"
    )
  output$shap_importance <- testing_shap %>%
    group_by(outer_split, variable) %>%
    summarise(split_mean_abs_shap = mean(abs(shap)), .groups = "drop") %>%
    complete(
      outer_split = seq_len(repeats),
      variable = predictors,
      fill = list(split_mean_abs_shap = 0)
    ) %>%
    group_by(variable) %>%
    summarise(
      mean_abs_shap = mean(split_mean_abs_shap),
      shap_outer_splits = sum(split_mean_abs_shap > 0),
      .groups = "drop"
    ) %>%
    left_join(feature_summary, by = "variable") %>%
    arrange(desc(mean_abs_shap))
  output
}
