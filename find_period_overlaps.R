

library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

# ============================================================
# SETTINGS
# ============================================================
window_length     <- 15  # length of the period, in years
min_years_sampled <- 12  # minimum sampled years required within that period
max_gap_years     <- 2   # longest allowed gap between consecutive sampled
# years within the window (e.g. 3 = no more than
# 3 years in a row with no sample)


# ============================================================
# STEP 1: Get distinct sampled years per site
# Replace `raw` with your real long-format data (site_id, date)
# ============================================================

raw <- data.frame(
  site_id = c(rep("A", 18), rep("B", 15), rep("C", 10), rep("D", 20)),
  date = as.Date(c(
    seq(as.Date("2000-01-01"), as.Date("2019-01-01"), by = "year")[1:18],
    seq(as.Date("1998-01-01"), as.Date("2015-01-01"), by = "year")[1:15],
    seq(as.Date("1985-01-01"), as.Date("1996-01-01"), by = "year")[1:10],
    seq(as.Date("2001-01-01"), as.Date("2022-01-01"), by = "year")[1:20]
  ))
)

raw <- annual_dat %>% 
  # used to select streams that had 20 year records for Si, N and P
  #filter(Stream_Name %in% long_dat$Stream_Name) %>% 
  mutate(chemical = case_when(chemical == "NO3" ~ "NO3",
                              chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical == "P") %>%
  select(Stream_Name,Year,Conc_mgL) %>% 
  rename(site_id = Stream_Name,
         year = Year)

site_years <- raw %>%
  #mutate(year = year(date)) %>%
  distinct(site_id, year)


# ============================================================
# STEP 2: Build a site x year matrix (1 = sampled, 0 = not)
# Columns cover the FULL year range across the whole dataset,
# so windows can be evaluated consistently
# ============================================================

full_year_range <- seq(min(site_years$year), max(site_years$year))

year_matrix <- site_years %>%
  mutate(sampled = 1L) %>%
  pivot_wider(names_from = year, values_from = sampled, values_fill = 0L)

# Ensure every year in the full range has a column, even if no site
# has any samples that year (fills with 0)
missing_yr_cols <- setdiff(as.character(full_year_range), names(year_matrix))
for (yr in missing_yr_cols) year_matrix[[yr]] <- 0L

sites <- year_matrix$site_id
mat   <- as.matrix(year_matrix[, as.character(full_year_range)])
rownames(mat) <- sites


# ============================================================
# STEP 3: Slide a window of length `window_length` across the
# full year range. For each window, a site qualifies if it has
# >= min_years_sampled sampled years WITHIN that window, AND no
# gap between consecutive sampled years longer than max_gap_years.
# ============================================================

# Helper: given a window's start/end years and a site's sampled years
# within that window, return the longest internal gap, the leading gap
# (window start -> first sample), and the trailing gap (last sample ->
# window end). All three must be checked against max_gap_years.
gap_stats <- function(window_start, window_end, years) {
  if (length(years) == 0) {
    return(list(internal_gap = Inf, leading_gap = Inf, trailing_gap = Inf))
  }
  years <- sort(years)
  internal_gap <- if (length(years) > 1) max(diff(years) - 1) else 0
  leading_gap  <- years[1] - window_start
  trailing_gap <- window_end - years[length(years)]
  list(internal_gap = internal_gap, leading_gap = leading_gap, trailing_gap = trailing_gap)
}

possible_starts <- full_year_range[1]:(tail(full_year_range, 1) - window_length + 1)

window_results <- map_dfr(possible_starts, function(start_yr) {
  end_yr <- start_yr + window_length - 1
  window_cols <- as.character(start_yr:end_yr)
  window_cols <- window_cols[window_cols %in% colnames(mat)]
  window_yrs_numeric <- as.numeric(window_cols)
  
  # per-site count of sampled years, and gap checks (internal + edges)
  site_stats <- map_dfr(sites, function(s) {
    sampled_yrs_in_window <- window_yrs_numeric[mat[s, window_cols] == 1]
    g <- gap_stats(start_yr, end_yr, sampled_yrs_in_window)
    tibble(
      site_id = s,
      n_sampled     = length(sampled_yrs_in_window),
      internal_gap  = g$internal_gap,
      leading_gap   = g$leading_gap,
      trailing_gap  = g$trailing_gap
    )
  })
  
  qualifying <- site_stats %>%
    filter(n_sampled >= min_years_sampled,
           internal_gap <= max_gap_years,
           leading_gap  <= max_gap_years,
           trailing_gap <= max_gap_years) %>%
    pull(site_id)
  
  tibble(
    window_start     = start_yr,
    window_end       = end_yr,
    n_qualifying     = length(qualifying),
    qualifying_sites = list(qualifying)
  )
})

print(window_results %>% select(window_start, window_end, n_qualifying))


# ============================================================
# STEP 4: Identify the best window(s) -- the one(s) with the
# most qualifying sites. (There may be ties across adjacent or
# distant windows -- all are shown.)
# ============================================================

best_n <- max(window_results$n_qualifying)
best_windows <- window_results %>% filter(n_qualifying == best_n)

cat("Best window(s) found, with", best_n, "qualifying sites:\n")
print(best_windows %>% select(window_start, window_end, n_qualifying))

# View the qualifying sites for the (first) best window:
best_windows$qualifying_sites[[1]]


# ============================================================
# OPTIONAL: full detail for the best window -- how many years
# each qualifying site was actually sampled within it
# ============================================================

best_start <- best_windows$window_start[1]
best_end   <- best_windows$window_end[1]
best_cols  <- as.character(best_start:best_end)
best_cols  <- best_cols[best_cols %in% colnames(mat)]
best_cols_numeric <- as.numeric(best_cols)

site_detail <- map_dfr(sites, function(s) {
  sampled_yrs <- best_cols_numeric[mat[s, best_cols] == 1]
  g <- gap_stats(best_start, best_end, sampled_yrs)
  tibble(
    site_id = s,
    years_sampled_in_window = length(sampled_yrs),
    internal_gap = g$internal_gap,
    leading_gap  = g$leading_gap,
    trailing_gap = g$trailing_gap
  )
}) %>%
  filter(years_sampled_in_window >= min_years_sampled,
         internal_gap <= max_gap_years,
         leading_gap  <= max_gap_years,
         trailing_gap <= max_gap_years) %>%
  arrange(desc(years_sampled_in_window), internal_gap)

print(site_detail)