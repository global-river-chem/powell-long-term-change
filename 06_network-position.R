suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

# Calculate each site's downstream position in the HydroRIVERS network


# Settings ----------------------------------------------------------------

project_folder <- "/Users/sidneybush/Library/CloudStorage/Box-Box/Sidney_Bush/SiSyn/long-term-change"

spatial_file <- file.path(
  project_folder,
  "data",
  "all-data_si-extract_3_20260629.csv"
)

site_reference_file <- file.path(
  project_folder,
  "data",
  "Site_Reference_Table - WRTDS_Reference_Table_LTER_V3.csv"
)

hydrorivers_folder <- file.path(project_folder, "data", "hydrorivers")
output_folder <- file.path(project_folder, "outputs", "network-position")

distance_cutoffs_km <- c(10, 25, 50, 100)
basin_coverage_cutoffs <- c(0.80, 0.90, 0.95)
snap_distance_to_review_km <- 10

if (!file.exists(spatial_file) || !file.exists(site_reference_file)) {
  stop("The spatial or site-reference input file is missing", call. = FALSE)
}

dir.create(hydrorivers_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
sf::sf_use_s2(TRUE)


# Clean names and coordinates ---------------------------------------------

# Remove extra spaces and turn blank text into missing values
clean_text <- function(x) {
  x <- trimws(gsub("\u00A0", " ", as.character(x), fixed = TRUE))
  x[x == ""] <- NA_character_
  x
}

# Use one project name when the input files use different labels
clean_lter_name <- function(x) {
  x <- clean_text(x)
  x[tolower(x) %in% c("swedish goverment", "swedish government", "sweden")] <- "Sweden"
  x[tolower(x) == "cameroon"] <- "Congo Basin"
  x[tolower(x) == "eastriversfa"] <- "Coal Creek"
  x
}

# Apply known stream-name differences between the two input files
clean_stream_name <- function(x) {
  x <- clean_text(x)
  x[tolower(x) %in% c("east fork", "eastfork")] <- "east fork"
  x[tolower(x) %in% c("west fork", "westfork")] <- "west fork"

  dplyr::recode(
    x,
    "Amazon River at Obidos" = "Obidos",
    "MGWEIR" = "MG_WEIR",
    "ORlow" = "OR_low",
    "OR_WEIR" = "OR_low",
    "coal_11" = "Coal Creek",
    .default = x
  )
}

# Ignore capitalization, spaces, accents, and punctuation when matching names
match_key <- function(x) {
  x <- iconv(clean_text(x), to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("[^a-z0-9]+", "", x)
}

parse_coordinate <- function(x) {
  x <- gsub("\u2212", "-", as.character(x), fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

# Use each site location to choose a regional HydroRIVERS file
hydrorivers_region <- function(latitude, longitude) {
  dplyr::case_when(
    !is.finite(latitude) | !is.finite(longitude) ~ NA_character_,
    !between(latitude, -90, 90) | !between(longitude, -180, 180) ~ NA_character_,
    latitude < -60 ~ NA_character_,
    latitude >= 58 & longitude >= -75 & longitude <= -5 ~ "gr",
    latitude >= 60 & longitude < -50 ~ "ar",
    longitude >= -170 & longitude < -25 & latitude >= 15 ~ "na",
    longitude >= -170 & longitude < -25 & latitude < 15 ~ "sa",
    longitude >= -20 & longitude < 60 & latitude < 30 ~ "af",
    longitude >= -20 & longitude < 60 & latitude >= 30 ~ "eu",
    longitude >= 60 & longitude < 180 & latitude < -12 ~ "au",
    longitude >= 60 & longitude < 180 & latitude >= 50 ~ "si",
    longitude >= 60 & longitude < 180 ~ "as",
    longitude < -170 ~ "na",
    TRUE ~ NA_character_
  )
}


# Read the site data -------------------------------------------------------

spatial_raw <- read.csv(spatial_file, check.names = TRUE)
site_reference_raw <- read.csv(site_reference_file, check.names = TRUE)

sites <- spatial_raw %>%
  transmute(
    LTER = clean_lter_name(LTER),
    Stream_Name = clean_stream_name(Stream_Name),
    Discharge_File_Name = clean_text(Discharge_File_Name),
    reference_drainage_area_km2 = suppressWarnings(as.numeric(drainage_area)),
    lter_key = match_key(LTER),
    stream_key = match_key(Stream_Name),
    discharge_key = match_key(Discharge_File_Name)
  )

site_coordinates <- site_reference_raw %>%
  transmute(
    LTER = clean_lter_name(LTER),
    Stream_Name = clean_stream_name(Stream_Name),
    Discharge_File_Name = clean_text(Discharge_File_Name),
    Latitude = parse_coordinate(Latitude),
    Longitude = parse_coordinate(Longitude),
    lter_key = match_key(LTER),
    stream_key = match_key(Stream_Name),
    discharge_key = match_key(Discharge_File_Name)
  ) %>%
  mutate(
    # Correct source-table coordinates confirmed from site records
    Latitude = case_when(
      LTER == "PIE" & Stream_Name == "Aberjona" ~ 42.4474568,
      LTER == "HYBAM" & tolower(Stream_Name) == "labrea" ~ -7.254421,
      LTER == "Sweden" & Stream_Name == "Raan Helsingborg" ~ 55.99957,
      TRUE ~ Latitude
    ),
    Longitude = case_when(
      LTER == "PIE" & Stream_Name == "Aberjona" ~ -71.1380816,
      LTER == "HYBAM" & tolower(Stream_Name) == "labrea" ~ -64.801287,
      LTER == "Sweden" & Stream_Name == "Raan Helsingborg" ~ 12.77950,
      TRUE ~ Longitude
    ),
    region_as_supplied = hydrorivers_region(Latitude, Longitude),
    region_if_swapped = hydrorivers_region(Longitude, Latitude),
    swap_coordinates =
      (abs(Latitude) > 90 & abs(Longitude) <= 90) |
      (is.na(region_as_supplied) & !is.na(region_if_swapped)),
    original_latitude = Latitude,
    Latitude = ifelse(swap_coordinates, Longitude, Latitude),
    Longitude = ifelse(swap_coordinates, original_latitude, Longitude)
  )

invalid_coordinates <- site_coordinates %>%
  filter(
    is.finite(Latitude),
    is.finite(Longitude),
    !between(Latitude, -90, 90) | !between(Longitude, -180, 180)
  )

valid_coordinates <- site_coordinates %>%
  filter(
    is.finite(Latitude),
    is.finite(Longitude),
    between(Latitude, -90, 90),
    between(Longitude, -180, 180)
  )


# Add latitude and longitude to each site ---------------------------------

# Use an exact name first when the reference table contains similar site names
exact_stream_coordinates <- valid_coordinates %>%
  filter(!is.na(lter_key), !is.na(Stream_Name)) %>%
  mutate(coordinate = paste(Latitude, Longitude)) %>%
  group_by(lter_key, Stream_Name) %>%
  filter(n_distinct(coordinate) == 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    lter_key,
    Stream_Name,
    exact_stream_latitude = Latitude,
    exact_stream_longitude = Longitude
  )

exact_discharge_coordinates <- valid_coordinates %>%
  filter(!is.na(lter_key), !is.na(Discharge_File_Name)) %>%
  mutate(coordinate = paste(Latitude, Longitude)) %>%
  group_by(lter_key, Discharge_File_Name) %>%
  filter(n_distinct(coordinate) == 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    lter_key,
    Discharge_File_Name,
    exact_discharge_latitude = Latitude,
    exact_discharge_longitude = Longitude
  )

# Match remaining rows after standardizing capitalization and punctuation
stream_coordinates <- valid_coordinates %>%
  filter(!is.na(lter_key), !is.na(stream_key)) %>%
  mutate(coordinate = paste(Latitude, Longitude)) %>%
  group_by(lter_key, stream_key) %>%
  filter(n_distinct(coordinate) == 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    lter_key,
    stream_key,
    stream_latitude = Latitude,
    stream_longitude = Longitude
  )

discharge_coordinates <- valid_coordinates %>%
  filter(!is.na(lter_key), !is.na(discharge_key)) %>%
  mutate(coordinate = paste(Latitude, Longitude)) %>%
  group_by(lter_key, discharge_key) %>%
  filter(n_distinct(coordinate) == 1) %>%
  slice(1) %>%
  ungroup() %>%
  select(
    lter_key,
    discharge_key,
    discharge_latitude = Latitude,
    discharge_longitude = Longitude
  )

sites <- sites %>%
  left_join(
    exact_stream_coordinates,
    by = c("lter_key", "Stream_Name"),
    na_matches = "never"
  ) %>%
  left_join(
    exact_discharge_coordinates,
    by = c("lter_key", "Discharge_File_Name"),
    na_matches = "never"
  ) %>%
  left_join(
    stream_coordinates,
    by = c("lter_key", "stream_key"),
    na_matches = "never"
  ) %>%
  left_join(
    discharge_coordinates,
    by = c("lter_key", "discharge_key"),
    na_matches = "never"
  ) %>%
  mutate(
    Latitude = coalesce(
      exact_stream_latitude,
      exact_discharge_latitude,
      stream_latitude,
      discharge_latitude
    ),
    Longitude = coalesce(
      exact_stream_longitude,
      exact_discharge_longitude,
      stream_longitude,
      discharge_longitude
    ),
    hydrorivers_region = hydrorivers_region(Latitude, Longitude)
  ) %>%
  select(
    -lter_key,
    -stream_key,
    -discharge_key,
    -ends_with("_latitude"),
    -ends_with("_longitude")
  )

missing_coordinates <- sites %>%
  filter(!is.finite(Latitude) | !is.finite(Longitude))

outside_hydrorivers_coverage <- sites %>%
  filter(is.finite(Latitude), is.finite(Longitude), is.na(hydrorivers_region))

sites_to_match <- sites %>%
  filter(!is.na(hydrorivers_region))


# Match the sites to the nearest HydroRIVERS segment ----------------------

# HydroRIVERS is distributed in regional files, so each region is read once
region_results <- list()

for (region in sort(unique(sites_to_match$hydrorivers_region))) {
  region_name <- paste0("HydroRIVERS_v10_", region)
  zip_file <- file.path(hydrorivers_folder, paste0(region_name, "_shp.zip"))
  unzip_folder <- file.path(hydrorivers_folder, paste0(region_name, "_shp"))
  shapefile <- file.path(
    unzip_folder,
    paste0(region_name, "_shp"),
    paste0(region_name, ".shp")
  )

  if (!file.exists(shapefile)) {
    if (!file.exists(zip_file)) {
      options(timeout = 3600)
      download.file(
        paste0(
          "https://data.hydrosheds.org/file/HydroRIVERS/",
          basename(zip_file)
        ),
        zip_file,
        mode = "wb"
      )
    }

    dir.create(unzip_folder, recursive = TRUE, showWarnings = FALSE)
    unzip(zip_file, exdir = unzip_folder)
  }

  rivers <- st_read(shapefile, quiet = TRUE) %>%
    select(
      HYRIV_ID,
      MAIN_RIV,
      DIST_DN_KM,
      UPLAND_SKM,
      ENDORHEIC,
      DIS_AV_CMS,
      ORD_STRA
    )

  region_sites <- sites_to_match %>%
    filter(hydrorivers_region == region) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)

  nearest_river <- st_nearest_feature(region_sites, rivers)
  matched_rivers <- rivers[nearest_river, ]

  outlet_row <- match(
    as.character(matched_rivers$MAIN_RIV),
    as.character(rivers$HYRIV_ID)
  )
  outlet_upland_skm <- as.numeric(rivers$UPLAND_SKM[outlet_row])

  snap_distance_km <- as.numeric(
    st_distance(region_sites, matched_rivers, by_element = TRUE)
  ) / 1000

  region_results[[region]] <- bind_cols(
    st_drop_geometry(region_sites),
    st_drop_geometry(matched_rivers)
  ) %>%
    mutate(
      downstream_to_outlet_km = as.numeric(DIST_DN_KM),
      downstream_to_ocean_km = ifelse(
        ENDORHEIC == 0,
        as.numeric(DIST_DN_KM),
        NA_real_
      ),
      upland_skm = as.numeric(UPLAND_SKM),
      outlet_upland_skm = outlet_upland_skm,
      basin_coverage_fraction = ifelse(
        outlet_upland_skm > 0,
        upland_skm / outlet_upland_skm,
        NA_real_
      ),
      reference_basin_coverage_fraction = ifelse(
        outlet_upland_skm > 0,
        reference_drainage_area_km2 / outlet_upland_skm,
        NA_real_
      ),
      discharge_est_cms = as.numeric(DIS_AV_CMS),
      drains_to_ocean = ENDORHEIC == 0,
      snap_distance_km = snap_distance_km
    )

  rm(rivers, region_sites, matched_rivers)
  gc(verbose = FALSE)
}

site_metrics <- bind_rows(region_results) %>%
  arrange(LTER, Stream_Name, Discharge_File_Name) %>%
  transmute(
    LTER,
    Stream_Name,
    Latitude,
    Longitude,
    HYRIV_ID,
    downstream_to_outlet_km,
    downstream_to_ocean_km,
    upland_skm,
    outlet_upland_skm,
    basin_coverage_fraction,
    reference_drainage_area_km2,
    reference_basin_coverage_fraction,
    discharge_est_cms,
    drains_to_ocean,
    ORD_STRA,
    snap_distance_km
  )


# Compare possible definitions of a near-mouth site -----------------------

has_distance_and_coverage <- site_metrics$drains_to_ocean &
  is.finite(site_metrics$downstream_to_ocean_km) &
  is.finite(site_metrics$basin_coverage_fraction)

cutoff_comparison <- expand.grid(
  distance_cutoff_km = distance_cutoffs_km,
  basin_coverage_cutoff = basin_coverage_cutoffs
) %>%
  rowwise() %>%
  mutate(
    n_all_matched_sites = nrow(site_metrics),
    n_ocean_connected_sites = sum(site_metrics$drains_to_ocean),
    n_sites_with_distance_and_coverage = sum(has_distance_and_coverage),
    n_near_mouth_high_coverage_sites = sum(
      has_distance_and_coverage &
        site_metrics$downstream_to_ocean_km <= distance_cutoff_km &
        site_metrics$basin_coverage_fraction >= basin_coverage_cutoff
    ),
    near_mouth_percent_of_sites_with_metrics =
      100 * n_near_mouth_high_coverage_sites /
      n_sites_with_distance_and_coverage
  ) %>%
  ungroup() %>%
  arrange(distance_cutoff_km, basin_coverage_cutoff)


# Save sites that need a coordinate or matching review --------------------

large_snap_distance <- site_metrics %>%
  filter(snap_distance_km > snap_distance_to_review_km)

sites_needing_review <- bind_rows(
  missing_coordinates %>%
    transmute(
      review_issue = "missing_site_coordinates",
      LTER,
      Stream_Name,
      Latitude,
      Longitude,
      snap_distance_km = NA_real_
    ),
  invalid_coordinates %>%
    transmute(
      review_issue = "invalid_reference_coordinates",
      LTER,
      Stream_Name,
      Latitude,
      Longitude,
      snap_distance_km = NA_real_
    ),
  outside_hydrorivers_coverage %>%
    transmute(
      review_issue = "outside_hydrorivers_coverage",
      LTER,
      Stream_Name,
      Latitude,
      Longitude,
      snap_distance_km = NA_real_
    ),
  large_snap_distance %>%
    transmute(
      review_issue = "large_snap_distance",
      LTER,
      Stream_Name,
      Latitude,
      Longitude,
      snap_distance_km
    )
) %>%
  distinct()


# Save the output files ----------------------------------------------------

write.csv(
  site_metrics,
  file.path(output_folder, "network_position_site_metrics.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  cutoff_comparison,
  file.path(output_folder, "network_position_cutoff_comparison.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  sites_needing_review,
  file.path(output_folder, "network_position_sites_needing_review.csv"),
  row.names = FALSE,
  na = ""
)

message("Network-position analysis finished: ", output_folder)
