# ---- Basic helpers ----

# Routine path, predictor, and model changes do not require edits here
normalize_stream_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[[:space:]]+", " ", x)
}

numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}

mean_or_na <- function(x) {
  x <- numeric_or_na(x)
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

median_or_na <- function(x) {
  x <- numeric_or_na(x)
  if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
}

mode_or_na <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

harmonize_permafrost <- function(x) {
  x <- numeric_or_na(x)
  percentages <- is.finite(x) & x > 1 & x <= 100
  x[percentages] <- x[percentages] / 100
  100 * pmin(pmax(x, 0), 1)
}

# ---- Lithology groups ----

grl_lithology_cluster <- function(major_rock, sedimentary_fraction) {
  rock <- tolower(trimws(as.character(major_rock)))
  consolidated <- case_when(
    rock %in% c(
      "sedimentary", "sedimentary; metamorphic",
      "sedimentary; carbonate_evaporite",
      "sedimentary; volcanic; carbonate_evaporite",
      "sedimentary; plutonic; metamorphic; carbonate_evaporite"
    ) ~ "Sedimentary",
    rock %in% c("volcanic", "plutonic; volcanic") ~ "Volcanic",
    rock %in% c(
      "plutonic", "plutonic; metamorphic",
      "plutonic; volcanic; metamorphic"
    ) ~ "Plutonic",
    rock %in% c(
      "metamorphic", "metamorphic; carbonate_evaporite"
    ) ~ "Metamorphic",
    rock %in% c(
      "carbonate_evaporite", "volcanic; carbonate_evaporite"
    ) ~ "Carbonate Evaporite",
    TRUE ~ NA_character_
  )
  final <- case_when(
    consolidated == "Sedimentary" &
      is.finite(sedimentary_fraction) & sedimentary_fraction >= 70 ~
        "Sedimentary",
    consolidated == "Sedimentary" ~ "Mixed Sedimentary",
    TRUE ~ consolidated
  )
  final
}

# ---- Site-average spatial predictors ----

build_spatial_predictors <- function(path, predictor_spec) {
  spatial <- read.csv(path, check.names = FALSE)
  composition <- c(predictor_spec$land, predictor_spec$rock)

  for (variable in composition) {
    if (!variable %in% names(spatial)) spatial[[variable]] <- NA_real_
    spatial[[variable]] <- numeric_or_na(spatial[[variable]])
    spatial[[variable]][is.na(spatial[[variable]])] <- 0
  }

  site_rows <- spatial %>%
    transmute(
      stream_key = normalize_stream_key(Stream_Name),
      Stream_Name_spatial = Stream_Name,
      LTER_spatial = LTER,
      permafrost = harmonize_permafrost(permafrost_mean_m),
      elevation = numeric_or_na(elevation_mean_m),
      basin_slope = numeric_or_na(basin_slope_mean_degree),
      RBI = numeric_or_na(RBI),
      recession_slope = numeric_or_na(RCS),
      major_rock = as.character(major_rock),
      final_cluster = grl_lithology_cluster(
        major_rock, numeric_or_na(rocks_sedimentary)
      ),
      across(all_of(composition), numeric_or_na)
    )

  spatial_predictors <- predictor_spec$static
  site_rows %>%
    group_by(stream_key) %>%
    summarise(
      Stream_Name_spatial = mode_or_na(Stream_Name_spatial),
      LTER_spatial = paste(sort(unique(LTER_spatial)), collapse = ";"),
      major_rock = mode_or_na(major_rock),
      final_cluster = mode_or_na(final_cluster),
      across(all_of(spatial_predictors), mean_or_na),
      .groups = "drop"
    )
}

build_environment_clusters <- function(path) {
  read.csv(path, check.names = TRUE) %>%
    transmute(
      stream_key = normalize_stream_key(Stream_Name),
      environment_cluster = as.character(cluster),
      environment_cluster_name = as.character(cluster_name)
    ) %>%
    group_by(stream_key) %>%
    summarise(
      environment_cluster = mode_or_na(environment_cluster),
      environment_cluster_name = mode_or_na(environment_cluster_name),
      .groups = "drop"
    )
}

# ---- Annual driver table ----

wide_driver_series <- function(data, prefix, suffix, variable) {
  pattern <- paste0("^", prefix, "([0-9]{4})", suffix, "$")
  columns <- grep(pattern, names(data), value = TRUE)
  if (length(columns) == 0) stop("No annual columns found for ", variable)

  output <- data %>%
    transmute(
      stream_key = normalize_stream_key(Stream_Name),
      across(all_of(columns))
    ) %>%
    pivot_longer(
      cols = all_of(columns), names_to = "source_column", values_to = "value"
    ) %>%
    mutate(Year = as.integer(sub(pattern, "\\1", source_column)))

  output$value <- numeric_or_na(output$value)

  output %>%
    group_by(stream_key, Year) %>%
    summarise(value = mean_or_na(value), .groups = "drop") %>%
    rename(!!variable := value)
}

build_spatial_annual_drivers <- function(path) {
  spatial <- read.csv(path, check.names = FALSE)
  specifications <- list(
    c("npp_", "_kgC_m2_year", "npp"),
    c("evapotrans_", "_kg_m2", "evapotrans"),
    c("precip_", "_mm_per_day", "precip"),
    c("temp_", "_degC", "temp"),
    c("snow_", "_max_prop_area", "snow_cover")
  )
  tables <- lapply(specifications, function(specification) {
    wide_driver_series(
      spatial,
      prefix = specification[[1]],
      suffix = specification[[2]],
      variable = specification[[3]]
    )
  })
  Reduce(
    function(x, y) full_join(x, y, by = c("stream_key", "Year")),
    tables
  )
}

build_raw_annual_nutrients <- function(path) {
  data.table::fread(
    path,
    select = c("Stream_Name", "date", "variable", "value"),
    showProgress = FALSE
  ) %>%
    as_tibble() %>%
    filter(
      variable %in% c("NO3", "NOx", "SRP", "PO4"),
      is.finite(value),
      value > 0
    ) %>%
    mutate(
      stream_key = normalize_stream_key(Stream_Name),
      Year = as.integer(substr(date, 1, 4)),
      nutrient = if_else(variable %in% c("NO3", "NOx"), "N", "P")
    ) %>%
    filter(is.finite(Year)) %>%
    group_by(stream_key, Year, nutrient) %>%
    summarise(value = median(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = nutrient, values_from = value) %>%
    rename(N_raw = N, P_raw = P)
}

build_annual_chemistry_drivers <- function(path, raw_path) {
  chemistry <- read.csv(path, check.names = FALSE) %>%
    mutate(stream_key = normalize_stream_key(Stream_Name))

  nutrients <- chemistry %>%
    filter(
      chemical %in% c("NO3", "NOx", "P"),
      is.finite(GenConc_uM),
      GenConc_uM > 0
    ) %>%
    mutate(nutrient = if_else(chemical %in% c("NO3", "NOx"), "N", "P")) %>%
    group_by(stream_key, Year, nutrient) %>%
    summarise(value = median(GenConc_uM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = nutrient, values_from = value) %>%
    rename(N_wrtds = N, P_wrtds = P)
  raw_nutrients <- build_raw_annual_nutrients(raw_path)
  nutrients <- full_join(
    nutrients, raw_nutrients, by = c("stream_key", "Year")
  ) %>%
    transmute(
      stream_key,
      Year,
      N = coalesce(N_wrtds, N_raw),
      P = coalesce(P_wrtds, P_raw)
    )

  discharge <- chemistry %>%
    filter(
      chemical == "DSi",
      is.finite(Discharge_cms),
      is.finite(drainSqKm),
      drainSqKm > 0
    ) %>%
    group_by(stream_key, Year) %>%
    summarise(
      qnorm = median(Discharge_cms / drainSqKm, na.rm = TRUE),
      .groups = "drop"
    )

  full_join(nutrients, discharge, by = c("stream_key", "Year"))
}

build_annual_driver_table <- function(input_files) {
  spatial <- build_spatial_annual_drivers(input_files[["spatial"]])
  chemistry <- build_annual_chemistry_drivers(
    input_files[["chemistry"]], input_files[["raw_chemistry"]]
  )
  full_join(spatial, chemistry, by = c("stream_key", "Year"))
}

# ---- Driver trends ----

sen_slope_by_year <- function(year, value, minimum_years) {
  keep <- is.finite(year) & is.finite(value)
  year <- numeric_or_na(year[keep])
  value <- numeric_or_na(value[keep])
  unique_years <- length(unique(year))
  if (unique_years < minimum_years) return(NA_real_)

  pairs <- utils::combn(seq_along(year), 2)
  year_difference <- year[pairs[2, ]] - year[pairs[1, ]]
  valid <- year_difference != 0
  median(
    (value[pairs[2, valid]] - value[pairs[1, valid]]) /
      year_difference[valid],
    na.rm = TRUE
  )
}

summarise_driver_trends <- function(aligned, driver_variables,
                                    minimum_years, prefix) {
  aligned %>%
    group_by(stream_key) %>%
    summarise(
      across(
        all_of(driver_variables),
        ~ sen_slope_by_year(Year, .x, minimum_years),
        .names = paste0(prefix, "{.col}")
      ),
      .groups = "drop"
    )
}

build_driver_summaries <- function(annual_drivers, site_keys, settings) {
  driver_variables <- setdiff(
    names(annual_drivers), c("stream_key", "Year")
  )
  analysis_years <- tidyr::crossing(
    stream_key = site_keys,
    Year = seq.int(settings$analysis_start_year, settings$analysis_end_year)
  )
  aligned <- analysis_years %>%
    left_join(annual_drivers, by = c("stream_key", "Year"))
  averages <- aligned %>%
    group_by(stream_key) %>%
    summarise(
      median_N = median_or_na(N),
      median_P = median_or_na(P),
      npp = mean_or_na(npp),
      evapotrans = mean_or_na(evapotrans),
      precip = mean_or_na(precip),
      temp = mean_or_na(temp),
      snow_cover = mean_or_na(snow_cover),
      qnorm = mean_or_na(qnorm),
      .groups = "drop"
    )
  trends <- summarise_driver_trends(
    aligned,
    driver_variables,
    settings$minimum_driver_observations,
    "trend_"
  )
  coverage <- aligned %>%
    select(stream_key, Year, all_of(driver_variables)) %>%
    pivot_longer(
      cols = all_of(driver_variables),
      names_to = "driver",
      values_to = "value"
    ) %>%
    group_by(stream_key, driver) %>%
    summarise(years = n_distinct(Year[is.finite(value)]), .groups = "drop")

  list(
    values = averages %>%
      full_join(trends, by = "stream_key") %>%
      mutate(
        predictor_start_year = settings$analysis_start_year,
        predictor_end_year = settings$analysis_end_year,
        predictor_years = predictor_end_year - predictor_start_year + 1L
      ),
    coverage = coverage
  )
}

# ---- Site-level predictor table ----

build_site_predictors <- function(input_files, predictor_spec, settings) {
  spatial <- build_spatial_predictors(
    input_files[["spatial"]], predictor_spec
  )
  environment_clusters <- build_environment_clusters(
    input_files[["environment_clusters"]]
  )
  annual_drivers <- build_annual_driver_table(input_files)
  driver_summaries <- build_driver_summaries(
    annual_drivers, spatial$stream_key, settings
  )

  data <- spatial %>%
    left_join(environment_clusters, by = "stream_key") %>%
    left_join(driver_summaries$values, by = "stream_key") %>%
    select(
      stream_key, Stream_Name_spatial, LTER_spatial,
      major_rock, final_cluster, environment_cluster,
      environment_cluster_name,
      predictor_start_year, predictor_end_year, predictor_years,
      all_of(predictor_spec$candidates)
  )

  list(data = data, trend_coverage = driver_summaries$coverage)
}

# ---- Response data ----

load_slope_response <- function(path, target_chemical) {
  read.csv(path, check.names = FALSE) %>%
    filter(chemical == target_chemical, is.finite(estimates)) %>%
    transmute(
      stream_key = normalize_stream_key(Stream_Name),
      Stream_Name,
      response = numeric_or_na(estimates),
      p_value = numeric_or_na(p.value),
      significant = is.finite(p_value) & p_value <= 0.05
    )
}

prepare_model_data <- function(slopes, site_predictors, spec,
                               predictor_spec) {
  matched_data <- slopes %>%
    inner_join(site_predictors$data, by = "stream_key")
  candidates <- predictor_spec$candidates
  missing_fraction <- vapply(
    matched_data[candidates],
    function(x) mean(!is.finite(x)),
    numeric(1)
  )
  predictor_sets <- list(
    average_only = predictor_spec$average,
    average_plus_trends = c(predictor_spec$average, predictor_spec$trends)
  )
  # Use one complete cohort so the predictor-set comparison stays paired
  missing_matrix <- !is.finite(
    as.matrix(matched_data[predictor_sets$average_plus_trends])
  )
  complete_rows <- rowSums(missing_matrix) == 0
  excluded_sites <- matched_data[!complete_rows, ] %>%
    transmute(
      stream_key,
      Stream_Name,
      missing_predictors = apply(
        missing_matrix[!complete_rows, , drop = FALSE],
        1,
        function(row) paste(names(row)[row], collapse = ", ")
      )
    )
  model_data <- matched_data[complete_rows, ]
  varies <- vapply(model_data[candidates], function(x) {
    length(unique(x)) > 1
  }, logical(1))
  if (any(!varies)) {
    stop(
      "Selected predictors have no variation: ",
      paste(names(varies)[!varies], collapse = ", "),
      call. = FALSE
    )
  }

  audit <- tibble(
    model_id = spec$model_id,
    data_type = spec$data_type,
    response = paste(spec$target_label, "Sen slope"),
    available_response_sites = nrow(slopes),
    matched_spatial_sites = nrow(matched_data),
    model_sites = nrow(model_data),
    unmatched_spatial_sites = nrow(slopes) - nrow(matched_data),
    incomplete_sites_removed = nrow(excluded_sites),
    significant_sites = sum(model_data$significant),
    nonsignificant_sites = sum(!model_data$significant),
    average_predictors = length(predictor_sets$average_only),
    trend_predictors = length(setdiff(
      predictor_sets$average_plus_trends, predictor_sets$average_only
    )),
    significance_filtered = FALSE,
    complete_case_filtered = TRUE
  )

  list(
    data = model_data,
    predictor_sets = predictor_sets,
    audit = audit,
    missingness = tibble(
      variable = names(missing_fraction),
      missing_fraction = missing_fraction
    ),
    excluded_sites = excluded_sites,
    trend_coverage = site_predictors$trend_coverage %>%
      filter(stream_key %in% model_data$stream_key)
  )
}
