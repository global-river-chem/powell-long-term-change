# ---- Plot labels ----

# Change the names displayed in plots in this vector
predictor_labels <- c(
  median_N = "Median annual N", median_P = "Median annual P", npp = "NPP",
  evapotrans = "ET", precip = "Precip",
  temp = "Temp", snow_cover = "Snow cover", qnorm = "Mean specific discharge",
  permafrost = "Permafrost",
  elevation = "Elevation", basin_slope = "Basin slope",
  RBI = "RBI", recession_slope = "RCS",
  land_Cropland = "Cropland cover", land_Forest = "Forest cover",
  land_Grassland_Shrubland = "Grass and shrub cover",
  land_Wetland_Marsh = "Wetland cover",
  rocks_volcanic = "Volcanic rock",
  rocks_carbonate_evaporite = "Carbonate-evaporite rock",
  rocks_metamorphic = "Metamorphic rock",
  rocks_plutonic = "Plutonic rock",
  trend_N = "N Sen slope", trend_P = "P Sen slope",
  trend_npp = "NPP Sen slope", trend_evapotrans = "ET Sen slope",
  trend_precip = "Precip Sen slope", trend_temp = "Temp Sen slope",
  trend_snow_cover = "Snow-cover Sen slope",
  trend_qnorm = "Specific discharge Sen slope"
)

pretty_predictor <- function(x) {
  labels <- unname(predictor_labels[x])
  labels[is.na(labels)] <- x[is.na(labels)]
  labels
}

# ---- Plot data ----

scale_zero_one <- function(x) {
  limits <- range(x, finite = TRUE)
  if (!all(is.finite(limits)) || diff(limits) == 0) return(rep(0.5, length(x)))
  (x - limits[[1]]) / diff(limits)
}

shap_plot_data <- function(result) {
  top_variables <- result$validation$shap_importance %>%
    slice_max(mean_abs_shap, n = 8, with_ties = FALSE) %>%
    pull(variable)

  result$validation$testing_shap_site %>%
    filter(variable %in% top_variables) %>%
    group_by(variable) %>%
    mutate(scaled_value = scale_zero_one(predictor_value)) %>%
    ungroup() %>%
    mutate(label = pretty_predictor(variable))
}

# ---- Performance panels ----

performance_panel <- function(result, tag) {
  predictions <- result$validation$site_predictions %>%
    left_join(result$site_data, by = c("stream_key", "Stream_Name", "observed" = "response"))
  metrics <- result$validation$outer_metrics
  limits <- range(c(predictions$observed, predictions$predicted), finite = TRUE)
  padding <- 0.04 * diff(limits)
  if (!is.finite(padding) || padding == 0) padding <- 1
  limits <- limits + c(-padding, padding)
  annotation <- sprintf(
    paste(
      "Median split predictive R² = %.2f",
      "Median split RMSE = %.3g",
      "Testing sites shown = %d/%d",
      sep = "\n"
    ),
    median(metrics$R2), median(metrics$RMSE),
    nrow(predictions), nrow(result$site_data)
  )

  ggplot(predictions, aes(predicted, observed)) +
    geom_abline(linetype = "dashed", linewidth = 0.7, color = "grey45") +
    geom_point(
      aes(fill = significant),
      shape = 21,
      size = 2.8,
      stroke = 0.6,
      color = "grey25",
      alpha = 0.85
    ) +
    annotate(
      "text", x = limits[[1]], y = limits[[2]], label = annotation,
      hjust = 0, vjust = 1, size = 3.8
    ) +
    scale_fill_manual(
      values = c(`TRUE` = "#525693", `FALSE` = "white"),
      labels = c(
        `TRUE` = "Significant slope",
        `FALSE` = "Not significant"
      ),
      name = NULL
    ) +
    coord_equal(xlim = limits, ylim = limits) +
    labs(
      title = result$specification$plot_label,
      subtitle = "Static + average drivers + 2002–2022 Sen slopes",
      tag = tag,
      x = "Repeated-split predicted Sen slope",
      y = "Observed Sen slope"
    )
}

# ---- SHAP panels ----

importance_panel <- function(result, tag) {
  importance <- result$validation$shap_importance %>%
    slice_max(mean_abs_shap, n = 8, with_ties = FALSE) %>%
    mutate(
      label = pretty_predictor(variable),
      label = reorder(label, mean_abs_shap)
    )

  ggplot(importance, aes(mean_abs_shap, label)) +
    geom_col(fill = "grey25", width = 0.72) +
    labs(
      tag = tag,
      x = "Mean absolute SHAP value",
      y = NULL
    )
}

shap_distribution_panel <- function(result, tag) {
  plot_data <- shap_plot_data(result)
  order <- result$validation$shap_importance %>%
    slice_max(mean_abs_shap, n = 8, with_ties = FALSE) %>%
    arrange(mean_abs_shap) %>%
    mutate(label = pretty_predictor(variable)) %>%
    pull(label)
  plot_data$label <- factor(plot_data$label, levels = order)

  set.seed(42)
  ggplot(plot_data, aes(shap, label)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.7) +
    geom_jitter(
      aes(fill = scaled_value),
      shape = 21,
      color = "grey55",
      height = 0.18,
      size = 2.3,
      alpha = 0.85
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      name = "Scaled predictor value"
    ) +
    labs(tag = tag, x = "SHAP value", y = NULL)
}

# ---- Figure assembly ----

create_figure2 <- function(results, output_path) {
  theme_set(
    theme_classic(base_size = 14) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        plot.tag = element_text(size = 18),
        plot.title = element_text(size = 16)
      )
  )

  data_types <- vapply(
    results,
    function(result) result$specification$data_type,
    character(1)
  )
  concentration <- results[[which(data_types == "concentration")]]
  yield <- results[[which(data_types == "yield")]]
  panel_a <- performance_panel(concentration, "a)")
  panel_b <- performance_panel(yield, "b)")
  panel_c <- importance_panel(concentration, "c)")
  panel_d <- importance_panel(yield, "d)")
  panel_e <- shap_distribution_panel(concentration, "e)")
  panel_f <- shap_distribution_panel(yield, "f)")

  significance_legend <- cowplot::get_legend(
    panel_a + theme(legend.position = "bottom")
  )
  shap_legend <- cowplot::get_legend(
    panel_e + theme(legend.position = "bottom")
  )
  row_performance <- cowplot::plot_grid(
    panel_a + theme(legend.position = "none"),
    panel_b + theme(legend.position = "none"),
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )
  row_importance <- cowplot::plot_grid(
    panel_c,
    panel_d,
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )
  row_shap <- cowplot::plot_grid(
    panel_e + theme(legend.position = "none"),
    panel_f + theme(legend.position = "none"),
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )
  figure <- cowplot::plot_grid(
    row_performance,
    significance_legend,
    row_importance,
    row_shap,
    shap_legend,
    ncol = 1,
    rel_heights = c(1.05, 0.10, 1.05, 1.10, 0.10)
  )
  footnote <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste(
        "Models use static + average drivers + 2002–2022 Sen slopes; SHAP",
        "explains testing predictions only; splits balance lithology and Si slopes"
      ),
      x = 0,
      hjust = 0,
      size = 10
    )
  figure <- cowplot::plot_grid(
    figure,
    footnote,
    ncol = 1,
    rel_heights = c(1, 0.025)
  )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    output_path,
    figure,
    width = 15,
    height = 18,
    dpi = 300,
    bg = "white"
  )
  invisible(figure)
}

# ---- Predictor-set comparison ----

create_predictor_comparison_figure <- function(comparison_results,
                                               output_path) {
  plot_data <- bind_rows(lapply(comparison_results, function(model_sets) {
    bind_rows(lapply(model_sets, function(result) {
      result$validation$outer_metrics %>%
        mutate(
          plot_label = result$specification$plot_label,
          predictor_set = result$predictor_set
        )
    }))
  })) %>%
    select(plot_label, predictor_set, outer_split, R2, RMSE) %>%
    pivot_longer(c(R2, RMSE), names_to = "metric", values_to = "value") %>%
    mutate(
      predictor_set = recode(
        predictor_set,
        average_only = "Static + average drivers",
        average_plus_trends = "Static + average drivers + Sen slopes"
      ),
      predictor_set = factor(
        predictor_set,
        levels = c(
          "Static + average drivers",
          "Static + average drivers + Sen slopes"
        )
      ),
      metric = recode(
        metric,
        R2 = "Testing predictive R²",
        RMSE = "Testing RMSE"
      ),
      panel = paste(plot_label, metric, sep = " — ")
    )
  zero_panels <- plot_data %>%
    filter(metric == "Testing predictive R²") %>%
    distinct(panel)

  plot <- ggplot(
    plot_data,
    aes(predictor_set, value, group = interaction(plot_label, outer_split))
  ) +
    geom_hline(
      data = zero_panels,
      aes(yintercept = 0),
      color = "grey60",
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_line(color = "grey70", linewidth = 0.6) +
    geom_point(aes(fill = predictor_set), shape = 21, size = 3) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("white", "#525693")) +
    labs(
      title = "Do driver Sen slopes improve prediction beyond site summaries?",
      subtitle = "Predictor sets use identical sites in every repeated split",
      x = NULL,
      y = NULL,
      fill = NULL,
      caption = "Each line is one paired 80/20 training/testing split"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 15, hjust = 1),
      plot.title.position = "plot"
    )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    output_path,
    plot,
    width = 11,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  invisible(plot)
}

# ---- Split sensitivity figure ----

create_split_sensitivity_figure <- function(results, output_path) {
  plot_data <- bind_rows(lapply(results, function(result) {
    result$metrics %>%
      mutate(plot_label = result$specification$plot_label)
  })) %>%
    select(plot_label, stratification, outer_split, split, R2, RMSE) %>%
    pivot_longer(c(R2, RMSE), names_to = "metric", values_to = "value") %>%
    mutate(
      split = factor(split, levels = c("60/40", "80/20")),
      metric = recode(
        metric,
        R2 = "Testing predictive R²",
        RMSE = "Testing RMSE"
      ),
      panel = paste(plot_label, metric, sep = " — ")
    )
  zero_panels <- plot_data %>%
    filter(metric == "Testing predictive R²") %>%
    distinct(panel)

  plot <- ggplot(plot_data, aes(split, value, fill = stratification)) +
    geom_hline(
      data = zero_panels,
      aes(yintercept = 0),
      color = "grey60",
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.75,
      position = position_dodge(width = 0.75)
    ) +
    geom_point(
      size = 1.6,
      alpha = 0.60,
      position = position_jitterdodge(
        jitter.width = 0.10,
        dodge.width = 0.75
      )
    ) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    scale_fill_manual(
      values = c(
        `Dominant lithology` = "#525693",
        `Environmental clusters` = "#b9d7ef"
      )
    ) +
    labs(
      title = "Sensitivity of testing performance to the site split",
      subtitle = paste(
        "Static + average drivers + Sen slopes; both methods balance observed",
        "slope distribution within their stratification groups"
      ),
      x = "Training/testing proportion",
      y = NULL,
      fill = "Split stratification",
      caption = "Each point is one repeated split; predictive R² is 1 − SSE/SST"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      plot.title.position = "plot"
    )

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    output_path,
    plot,
    width = 11,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  invisible(plot)
}
