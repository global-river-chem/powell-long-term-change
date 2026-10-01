## evaluate overlapping site years - draft code

library(tidyverse)
library(trend)
library(viridis)
library(googledrive)
library(ggridges)
library(cowplot)
library(viridis)
library(sf)
library(mapview)
library(purrr)
library(lubridate)
#install.packages("igraph") # if needed
library(igraph)

# directory for downloading files from google drive
(path <- scicomptools::wd_loc(local = FALSE, remote_path = file.path('/', "home","jankowski","data")))

# directory for accessing WRTDS results
(path_2 <- scicomptools::wd_loc(local = FALSE, remote_path = file.path('/', "home", "shares", "lter-si", "WRTDS","WRTDS Results_2025")))

# Data import -------------------------------------------------------------

## Reference table
ref_link<-"https://drive.google.com/file/d/1u4yNGOFOVKfQoEak7z9ovkFhtAm9XRE_/view?usp=drive_link"
ref_folder = drive_get(as_id(ref_link))
ref <- drive_download(file = ref_folder$id, path = file.path(path,"Site_Reference_Table.csv"), overwrite = T)
ref_dat<-read_csv(ref$local_path)

# Annual data - generalized flow normalization 
annual_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_annual.csv"))
kalman_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_kalman_annual.csv"))


# Evaluate overlapping records --------------------------------------------

## OPTION A - one row with just start/end date. Data has no gaps.
# Function: years of overlap between two date ranges (0 if none)

# Set up dataframe - select for one or more chemicals
df <- duration_dat %>% filter(chemical=="DSi")
colnames(df) <- c("Stream_Name","chemical","LTER","start_date","end_date","duration")

overlap_years <- function(start1, end1, start2, end2) {
  ov_start <- pmax(start1, start2)
  ov_end   <- pmin(end1, end2)
  ov_years  <- as.numeric(ov_end - ov_start)
  ifelse(ov_years <= 0, 0, ov_years)
}

# All pairwise combinations of sites
pairs <- combn(df$Stream_Name, 2, simplify = FALSE)
glimpse(pairs)

# look across all pairs to see which have overlapping years
overlap_df <- map_dfr(pairs, function(p) {
  s1 <- filter(df, Stream_Name == p[1])
  s2 <- filter(df, Stream_Name == p[2])
  tibble(
    site1 = p[1],
    site2 = p[2],
    overlap_years = overlap_years(s1$start_date, s1$end_date,
                                  s2$start_date, s2$end_date)
  )
})

overlap_df$overlap_years_2 <- overlap_df$overlap_years*365

# Keep only pairs with at least 15 years of overlap
overlap_15 <- overlap_df %>% filter(overlap_years_2 >= 15)
print(overlap_15)

# Which individual sites have at least one qualifying partner?
sites_with_overlap <- unique(c(overlap_15$site1, overlap_15$site2))
df_screened <- df %>% filter(site_id %in% sites_with_overlap)
print(df_screened)


## OPTION B - raw records, overlap based on actual sampled years, not just full time span
raw <- annual_dat %>% filter(chemical=="DSi") %>% select(LTER,Stream_Name,Year,chemical)

#Get the set of sampled years per site
site_years <- raw %>%
  group_by(Stream_Name) %>%
  summarise(years = list(unique(Year)), .groups = "drop")

# Pairwise overlap = number of years both sites were sampled
pairs_raw <- combn(site_years$Stream_Name, 2, simplify = FALSE)

overlap_raw_df <- map_dfr(pairs_raw, function(p) {
  y1 <- site_years$years[site_years$Stream_Name == p[1]][[1]]
  y2 <- site_years$years[site_years$Stream_Name == p[2]][[1]]
  shared <- sort(intersect(y1, y2))
  tibble(
    site1 = p[1],
    site2 = p[2],
    overlap_years = length(shared),
    shared_years  = list(shared)   # the actual years both sites were sampled
  )
})

overlap_raw_15 <- overlap_raw_df %>% filter(overlap_years >= 15)
print(overlap_raw_15)

# shared_years is a list-column (each cell holds a vector of years).
# To view the actual years for one pair:
overlap_raw_15$shared_years[[1]]

# To "unnest" into a long format (one row per site-pair-year):
overlap_raw_long <- overlap_raw_15 %>%
  select(site1, site2, shared_years) %>%
  unnest(shared_years) %>%
  rename(year = shared_years)
print(overlap_raw_long)





# now filter out set of sites that all  have overlapping time period
g <- graph_from_data_frame(overlap_raw_15[, c("site1", "site2")], directed = FALSE)

# first need to check for "clique density"
cat("Number of sites with >=1 qualifying overlap:", vcount(g), "\n")
cat("Number of qualifying pairs (edges):", ecount(g), "\n")
cat("Graph density (0 = sparse, 1 = fully connected):",
    round(edge_density(g), 3), "\n")

# Use max_cliques() rather than cliques():
# - cliques() returns EVERY clique, including small ones nested inside
#   bigger ones -- this is what explodes in dense graphs
# - max_cliques() returns only MAXIMAL cliques (not contained in any
#   larger qualifying group), which is almost always what you actually
#   want, and is a far smaller, faster result
max_cliques_15 <- max_cliques(g, min = 3)
cat("Number of maximal cliques (size 3+) found:", length(max_cliques_15), "\n")
print(max_cliques_15)

# For each clique, find the years shared by ALL member sites
# (requires site_years from Option B -- the actual sampled years per site)
clique_shared_years <- map_dfr(seq_along(max_cliques_15), function(i) {
  clique_sites <- names(max_cliques_15[[i]])  # site IDs in this clique
  
  # pull each site's set of sampled years
  years_list <- site_years$years[match(clique_sites, site_years$Stream_Name)]
  
  # intersect across ALL sites in the clique, not just pairs
  shared <- sort(Reduce(intersect, years_list))
  
  tibble(
    clique_id    = i,
    n_sites      = length(clique_sites),
    sites        = paste(clique_sites, collapse = ", "),
    n_shared_yrs = length(shared),
    shared_years = list(shared)
  )
})

# Keep only cliques where the FULL group still shares >=15 years
# (note: pairwise overlap >=15yrs does NOT guarantee the whole clique does --
# adding a 3rd+ site can only shrink or maintain the shared-year count)
clique_shared_years_15 <- clique_shared_years %>% filter(n_shared_yrs >= 16)
print(clique_shared_years_15)


#View the actual shared years for one clique:
clique_shared_years_15$shared_years[[10]]

# Long format: one row per clique-year
clique_years_long <- clique_shared_years_15 %>%
  select(clique_id, sites, shared_years) %>%
  unnest(shared_years) %>%
  rename(year = shared_years)
print(clique_years_long)

# If density is high (say, > 0.3-0.4) and max_cliques() is still slow,
# a practical alternative is to just find the SINGLE largest group:
largest_group <- largest_cliques(g)
print(largest_group)

