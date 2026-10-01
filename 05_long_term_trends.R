## looking at long term trends

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

## Monthly data
month_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_monthly.csv"))
month_kalman_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_kalman_monthly.csv"))

# Trend data
trend_link <- "https://drive.google.com/file/d/1dHCrqcF1scgBnk3lSzVOg3RRZSLBlqSr/view?usp=drive_link"
trend_folder <- drive_get(as_id(trend_link))
trend <- drive_download(file=trend_folder$id, path = file.path(path, "Full_Results_WRTDS_trends.csv"), overwrite=T)
trend_dat <- read_csv(trend$local_path)

## spatial data
spatial_link<- "https://drive.google.com/file/d/1mwpcPdxjPVSyempGG9-9yECSdHxdxhoh/view?usp=drive_link"
spatial_folder <- drive_get(as_id(spatial_link))
spatial <- drive_download(file=spatial_folder$id, path=file.path(path,"spatial_predictors_glass_climate_models_glc_20260919.csv"),overwrite=TRUE)
spatial_dat <- read_csv(spatial$local_path)

spatial_dat_check <- spatial_dat %>% 
  select(contains("land"))

## Coastal sites
coastal_link <- "https://drive.google.com/file/d/1Xf8ErPzJJJBwR230NSTQrH5tAH4U0Ovb/view?usp=drive_link"
coastal_folder <- drive_get(as_id(coastal_link))
coastal <- drive_download(file=coastal_folder$id, path = file.path(path,"coastal_sites_list.csv"),overwrite = TRUE)
coastal_list <- read_csv(coastal$local_path)

# subset annual data to only coastal rivers
coastal_rivers <- unique(coastal_list$Stream_Name)
annual_dat_coast <- annual_dat %>% 
  filter(Stream_Name %in% coastal_rivers | LTER == "Sweden"|LTER == "Seine"|Stream_Name == "COLUMBIA RIVER AT PORT WESTWARD")


## cluster data
# spatial data
cluster_link<-"https://drive.google.com/file/d/13_cncYkrEw4ZgSY15TGMjMVivs_sr9d5/view?usp=drive_link"
cluster_folder = drive_get(as_id(cluster_link))
cluster <- drive_download(file = cluster_folder$id, path = file.path(path,"Si_sites_cluster_eight.csv"), overwrite = T)
cluster_dat<-read_csv(cluster$local_path)

cluster_dat$cluster_name <- as.factor(cluster_dat$cluster_name)
glimpse(cluster_dat)


# sites select ------------------------------------------------------------

# summarizing record lengths by stream and chemical
duration_dat <- annual_dat %>% 
  # just a check to be sure all are still in reference table
  filter(Stream_Name %in% ref_dat$Stream_Name) %>% 
  # only select needed columns
  select(LTER,Stream_Name,Year,chemical,Conc_mgL) %>% 
  group_by(Stream_Name,chemical,LTER) %>% 
  summarise(min_year = min(Year),
            max_year = max(Year),
            duration = (max_year - min_year)+1)

# check for long records for main solutes (ratios regenerated below)
chems <- c("NO3","NOx","P","DSi")
n_yrs <- 15

# select sites that have at least a 15 year record
duration_dat_wide <- duration_dat %>% 
  filter(chemical %in% chems) %>% 
  select(-min_year,-max_year) %>%
  select(Stream_Name,chemical, duration,LTER) %>% 
  # combine NOx and NO3 as single variable (NO3)
  mutate(chemical = case_when(chemical == "NOx"|chemical == "NO3"~ "NO3",
                              .default = chemical)) %>% 
  # These had duplicate nitrate data (NOx and NO3) as identified by code below. 
  filter(!(Stream_Name == "MURRAY RIVER DOWNSTREAM YARRAWONGA WEIR" & chemical == "NO3"& duration == 17)) %>% 
  filter(!(Stream_Name == "BARWON RIVER AT MUNGINDI" & chemical == "NO3"& duration == 4)) %>% 
  filter(!(Stream_Name == "DARLING RIVER AT MENINDEE UPSTREAM WEIR 32" & chemical == "NO3"& duration == 13)) %>% 
  pivot_wider(names_from = "chemical", values_from="duration") %>% 
  mutate(long_si = case_when(DSi>=n_yrs ~ "yes",
                             DSi<n_yrs ~ "no", 
                             .default = NA),
         long_si_n = case_when(DSi >= n_yrs & NO3 >= n_yrs ~ "yes", 
                              DSi < n_yrs | NO3 < n_yrs ~ "no", 
                              .default = NA),
         long_si_p = case_when(DSi >= n_yrs &  P >=n_yrs ~ "yes", 
                              DSi < n_yrs |  P < n_yrs ~ "no", 
                              .default = NA),
         long_all = case_when(DSi >= n_yrs & NO3 >= n_yrs & P >=n_yrs ~ "yes", 
                              DSi < n_yrs | NO3 < n_yrs | P < n_yrs ~ "no", 
                              .default = NA))

## check for "hidden duplicates" 
duration_dat_wide %>%
  # Group by the columns you are using to reshape the data
  group_by(Stream_Name) %>% 
  # Filter for groups that have more than one row
  filter(n() > 1) %>% 
  ungroup()

## look at results - number of streams with chemical records
nrow(filter(duration_dat_wide, long_all == "yes"))
nrow(filter(duration_dat_wide, long_si == "yes"))
nrow(filter(duration_dat_wide, long_si_n == "yes"))
nrow(filter(duration_dat_wide, long_si_p == "yes"))




# list of sites with at least a 10 year record
annual_dat_coast_long <- annual_dat_coast %>% 
  group_by(Stream_Name,chemical) %>% 
  summarise(max_year = max(Year),
            min_year = min(Year),
            duration = (max_year - min_year)+1) %>% 
  filter(duration >= 15) %>% 
  distinct(Stream_Name)


## Select data records
# annual data with long-term coastal records selected
annual_dat_coast_v2 <- annual_dat_coast %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN" & chemical != "Si:P") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(Stream_Name %in% annual_dat_coast_long$Stream_Name)

# all long term data
annual_dat_long <- annual_dat %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN" & chemical != "Si:P") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>%  
  filter(Stream_Name %in% long$Stream_Name)


# set land cover groupings ------------------------------------------------



# checking land cover groupings
annual_dat_long_land <- annual_dat_long %>% 
  filter(chemical == "DSi") %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA)) %>% 
  distinct(Stream_Name,land_group) %>% 
  group_by(land_group) %>% 
  summarise(n=n())
  



# Sen slope estimates --------------------------------------------------------

modified_sens.slope <- function(x, ...) {
  result <- sens.slope(x, ...)
  tibble(
    p.value = result$p.value,
    statistic = result$statistic,
    estimates = result$estimates[1],
    low.conf = result$conf.int[1],
    high.conf = result$conf.int[2])
}

modified_mk_test <- function(x, ...) {
  result <- mk.test(x, ...)
  tibble(
    p.value = result$p.value,
    statistic = result$statistic,
    estimates = result$estimates[1])
}

# Analyze trends  ----------------------------------------------------------

## CONCENTRATION
conc_slope <- annual_dat_coast_v2 %>% 
  filter(!is.na(FNConc_uM)) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM))

conc_mk <- annual_dat_coast_v2 %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNConc_uM)) %>% 
  mutate(variable = rep("FNConc"))

## YIELD
yield_slope <- annual_dat_coast_v2 %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNYield_10_6kmol_yr_km2)) 

yield_mk <- annual_dat_coast_v2 %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNYield_10_6kmol_yr_km2)) %>% 
  mutate(variable = rep("FNYield"))


## DISCHARGE
dis_slope <- annual_dat_coast_v2 %>% 
  filter(chemical == "DSi") %>% 
  filter(!is.na(Discharge_cms)) %>% 
  group_by(Stream_Name) %>%
  group_modify(~ modified_sens.slope(.x$Discharge_cms)) %>% 
  mutate(variable = rep("Discharge"))


## RATIOS 
## need to adjust Si:N to be Si:NO3 because of DIN issue

# Calculating ratios
ratio_trend_dat <- annual_dat_coast_v2 %>% 
  select(LTER,Stream_Name,Year,chemical,FNConc_uM) %>% 
  filter(chemical == "NO3"|chemical == "DSi"|chemical == "NOx"|chemical == "P"|chemical == "Si:P") %>% 
  ## Calculate DIN (DIN = NOx <or> NO3 + NH4)
  # Handle "duplicate" values for sites that break across a year so have two values for one year
  dplyr::group_by(LTER, Stream_Name, chemical, Year) %>%
  dplyr::summarize(response_values = mean(FNConc_uM, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(chemical = dplyr::case_when(
    ### NOx is preferred for calculating DIN because it is NO3 + NOx
    chemical == "NOx" | chemical == "NO3" ~ "NO3x",
    .default = chemical)) %>% 
  tidyr::pivot_wider(names_from = chemical,
                     values_from = response_values,
                     values_fn = mean) %>% 
  ## Calculate ratios
  dplyr::mutate(Si_NO3x = ifelse(test = (!is.na(DSi) & !is.na(NO3x)),
                                 yes = (DSi / NO3x), no = NA)) %>% 
  dplyr::mutate(Si_P = ifelse(test = (!is.na(DSi) & !is.na(P)),
                                 yes = (DSi / P), no = NA)) %>% 
  pivot_longer(cols = DSi:Si_P,names_to = "chemical", values_to = "FNConc_uM")


# look at trends in ratios
conc_ratio_slope <- ratio_trend_dat %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical == "Si_NO3x"|chemical == "Si_P"|chemical == "Si:P") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM))

conc_ratio_mk <- ratio_trend_dat %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNConc_uM))

# Ratios of YIELD
ratio_trend_dat_yield <- annual_dat_coast_v2 %>% 
  select(LTER,Stream_Name,Year,chemical,FNYield_10_6kmol_yr_km2) %>% 
  filter(chemical == "NO3"|chemical == "DSi"|chemical == "NOx"|chemical == "P") %>% 
  ## Calculate DIN (DIN = NOx <or> NO3 + NH4)
  # Handle "duplicate" values for sites that break across a year so have two values for one year
  dplyr::group_by(LTER, Stream_Name, chemical, Year) %>%
  dplyr::summarize(response_values = mean(FNYield_10_6kmol_yr_km2, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(chemical = dplyr::case_when(
    ### NOx is preferred for calculating DIN because it is NO3 + NOx
    chemical == "NOx" | chemical == "NO3" ~ "NO3x",
    .default = chemical)) %>% 
  tidyr::pivot_wider(names_from = chemical,
                     values_from = response_values,
                     values_fn = mean) %>% 
  ## Calculate ratios
  dplyr::mutate(Si_NO3x = ifelse(test = (!is.na(DSi) & !is.na(NO3x)),
                                 yes = (DSi / NO3x), no = NA)) %>% 
  dplyr::mutate(Si_P = ifelse(test = (!is.na(DSi) & !is.na(P)),
                              yes = (DSi / P), no = NA)) %>% 
  pivot_longer(cols = DSi:Si_P,names_to = "chemical", values_to = "FNYield_10_6kmol_yr_km2")

yield_ratio_slope <- ratio_trend_dat_yield %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  filter(chemical == "Si_NO3x"|chemical == "Si_P") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNYield_10_6kmol_yr_km2))

yield_ratio_mk <- ratio_trend_dat_yield %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNYield_10_6kmol_yr_km2))


# export slope values -----------------------------------------------------

# bind datasets together for export
# concentration slopes

conc_slopes_export <- conc_slope %>% 
  filter(Stream_Name %in% annual_dat_long$Stream_Name) %>% 
  bind_rows(conc_ratio_slope) %>% 
  filter(chemical != "Si:P")

yield_slopes_export <- yield_slope %>% 
  filter(Stream_Name %in% annual_dat_coast_long$Stream_Name) %>% 
  bind_rows(yield_ratio_slope) %>% 
  filter(chemical != "Si:P")


write.csv(x = conc_slopes_export, row.names = F, na = "",
          file = file.path(path, 
                           "conc_slopes_export.csv"))

tidy_dest <- googledrive::as_id("https://drive.google.com/drive/u/1/folders/10HMZLr9TAf2asrprZMyyuiYlIesm4Jug")

## Export to it
googledrive::drive_upload(path = tidy_dest, overwrite = T,
                          media = file.path(path,
                                            "conc_slopes_export.csv"))

write.csv(x = yield_slopes_export, row.names = F, na = "",
          file = file.path(path, 
                           "yield_slopes_export.csv"))

tidy_dest <- googledrive::as_id("https://drive.google.com/drive/u/1/folders/10HMZLr9TAf2asrprZMyyuiYlIesm4Jug")

## Export to it
googledrive::drive_upload(path = tidy_dest, overwrite = T,
                          media = file.path(path,
                                            "yield_slopes_export.csv"))


# summarize direction of trends --------------------------------------------------------

# concentration of Si, N, and P - identified human impact vs not
conc_change <- conc_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))

# Yield of Si, N, and P
yield_change <- yield_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))


# Discharge
dis_change <- dis_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))


# Ratio of concentrations
conc_ratio_change <- conc_ratio_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))


# Ratio of yields
yield_ratio_change <- yield_ratio_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change"))  %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))



## Concentration 

conc_change_all <- bind_rows(conc_change,conc_ratio_change)

conc.agr = aggregate(Stream_Name~chemical*change*land_group, FUN=length, data=conc_change_all)
conc.agr1 = aggregate(Stream_Name~chemical*land_group, FUN=length, data=conc_change_all)
conc.agr2 <- conc.agr %>% 
  left_join(conc.agr1, by=c("chemical","land_group"))
conc.agr2$prop.streams <-  conc.agr2$Stream_Name.x/conc.agr2$Stream_Name.y

conc.agr2


conc.agr2 %>% 
  filter(chemical != "NH4") %>% 
  ggplot(aes(x = land_group, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Concentration")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  facet_wrap(~chemical, nrow=1)

## Yield 
yield_change_all <- bind_rows(yield_change,yield_ratio_change)


yield.agr = aggregate(Stream_Name~chemical*change*land_group, FUN=length, data=yield_change_all)
yield.agr1 = aggregate(Stream_Name~chemical*land_group, FUN=length, data=yield_change_all)
yield.agr2 <- yield.agr %>% 
  left_join(yield.agr1, by=c("chemical","land_group"))
yield.agr2$prop.streams <-  yield.agr2$Stream_Name.x/yield.agr2$Stream_Name.y


yield.agr2 %>% 
  filter(chemical != "NH4") %>% 
  ggplot(aes(x = land_group, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  facet_wrap(~chemical, nrow=1)

## Discharge
dis.agr = aggregate(Stream_Name~cluster_name*change, FUN=length, data=dis_cluster)
dis.agr1 = aggregate(Stream_Name~cluster_name, FUN=length, data=dis_cluster)
dis.agr2 <- dis.agr %>% 
  left_join(dis.agr1, by=c("cluster_name"))

dis.agr2$prop.streams <-  dis.agr2$Stream_Name.x/dis.agr2$Stream_Name.y

dis.agr2 %>% 
  ggplot(aes(x = cluster_name, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  ggtitle("Discharge")


# EGRET trends ------------------------------------------------------------

trend_dat_v1 <- trend_dat %>% 
  # coastal rivers
  filter(stream %in% coastal_rivers | LTER == "Sweden"|LTER == "Seine"|stream == "COLUMBIA RIVER AT PORT WESTWARD")
  #filter(stream %in% long$Stream_Name)

trend_conc <- trend_dat_v1 %>% 
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NOx", 
                              .default = chemical)) %>% 
  rename(Stream_Name = stream) %>% 
  filter(Metric == "Concentration") %>% 
  filter(chemical != "NH4") %>% 
  filter(!is.na(change_mg_L)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(Stream_Name,chemical) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=15) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA)) %>%  
  mutate(trend_direction = case_when(slope_percent_yr>0.1 ~ "increase",
                                     slope_percent_yr< -0.1 ~ "decrease",
                                     slope_percent_yr <= 0.1 & slope_percent_yr >= -0.1 ~ "no change"))

trend_conc %>% 
  group_by(chemical,trend_direction) %>% 
  summarise(n=n())

trend_conc %>%
  filter(!is.na(land_group)) %>% 
  filter(chemical != "NH4") %>% 
  # congo percentage is way off
  #filter(slope_percent_yr <20) %>% 
  ggplot(aes(land_group, slope_percent_yr,col=chemical))+
  geom_boxplot()+
  geom_jitter(aes(land_group, slope_percent_yr, fill=chemical), width = 0.2)+
  geom_hline(yintercept = 0, lty=2)+
  facet_wrap(~chemical, scales="free")+
  ylim(-2,2)

# Yield
trend_yield <- trend_dat_v1 %>% 
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NOx", 
                              .default = chemical)) %>% 
  rename(Stream_Name = stream) %>% 
  filter(Metric == "Flux") %>% 
  filter(chemical != "NH4") %>% 
  #filter(!is.na(change_mg_L)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(Stream_Name,chemical) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=15) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA)) %>%  
  mutate(trend_direction = case_when(slope_percent_yr>0.1 ~ "increase",
                                     slope_percent_yr< -0.1 ~ "decrease",
                                     slope_percent_yr <= 0.1 & slope_percent_yr >= -0.1 ~ "no change"))

trend_yield %>% 
  group_by(chemical,trend_direction) %>% 
  summarise(n=n())

trend_yield %>%
  filter(!is.na(land_group)) %>% 
  filter(slope_10_3kg_yr_yr < 2000& slope_10_3kg_yr_yr > -2000) %>% 
  group_by(chemical) %>% 
  summarise(mean = median(slope_10_3kg_yr_yr, na.rm=TRUE))

trend_yield %>% 
  filter(!is.na(land_group)) %>% 
  filter(slope_10_3kg_yr_yr < 2000& slope_10_3kg_yr_yr > -2000) %>% 
  ggplot(aes(land_group, slope_10_3kg_yr_yr))+
  geom_boxplot() +
  facet_wrap(~chemical, scales="free")

# grouping sites 
map_sites <- ref_dat %>%
  filter(!is.na(Latitude)) %>% 
  filter(LTER == "Sweden"|LTER == "NIVA"|LTER == "Finnish Environmental Institute")

si_sf <- st_as_sf(map_sites, coords = c("Longitude", "Latitude"), crs = 4326)
mapview(si_sf)


# assign ocean basins to sites
arctic <- c("Lena","Ob","Kolyma","Yenisey","Mackenzie","FINEALT", "FINETAN","FINEPAS")
n_atl <- c("NOREVEF","STREORK","STRENID","MROEDRI",
           "ST. LAWRENCE","CONNECTICUT RIVER AT THOMPSONVILLE","HUDSON RIVER","SUSQUEHANNA RIVER","POTOMAC RIVER",
           "EDISTO RIVER","ALTAMAHA RIVER","SOPCHOPPY RIVER","APALACHICOLA RIVER","Mississippi River near St. Francisville",
           "Brazos River near Rosharon")
s_atl <- c("Obidos","Saut Maripa","Langa Tabiki","Ciudad Bolivar", "Congo a Beach Brazzaville")
n_pac <- c("YUKON RIVER","Yukon","COLUMBIA RIVER AT PORT WESTWARD","SACRAMENTO R","SANTA ANA","SKEENA RIVER AT USK")
s_pac <- c("HUNTER RIVER AT SINGLETON")
indian <- c("Lock 5")
north_sea_lter <- c("UK","Seine")
north_sea_streams <- c("SFJENAU","HOREVOS","ROGEVIK","ROGEORR","ROGEBJE","VAGEOTR","AAGEVEG", "TELESKI",
                       "VESENUM","OSTEGLO","BUSEDRA","OSLEALN","Enningdalsalv N.Bullaren","Orekilsalven Munkedal","Gota Alv Trollhattan",
                       "Viskan Asbro","Atran Falkenberg","Nissan Halmstad","Lagan Laholm","Ronnean Klippan")
baltic_streams <- c("Helgean Hammarsjon","Morrumsan Morrum","Lyckebyan Lyckeby","Alsteran Getebro","Eman Emsfors",
                    "Motala Strom Norrkoping","Nykopingsan Spanga","Dalalven Alvkarleby","Gavlean Gavle stadspark",
                    "Ljusnan Ljusne Strommar","Delangersan Iggesund","Indalsalven Bergeforsen",
                    "Angermanalven Solleftea","Gide alv Gideabacka",'Ore alv Torrbole',"Ume alv Stornorrfors",
                    "Skellefte alv Kvistforsen","Pite alv Bolebyn","Lulealven","Raan Helsingborg","Kalixalven","Torne alv Mattila")
baltic_lter <- "Finnish Environmental Institute"


trend_conc <- trend_dat_v1 %>% 
  filter(Metric == "Concentration") %>% 
  filter(!is.na(change_mg_L)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(stream,chemical) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=10) %>% 
  mutate(trend_direction = case_when(slope_percent_yr>0.5 ~ "increase",
                                     slope_percent_yr< -0.5 ~ "decrease",
                                     slope_percent_yr <= 0.5 & slope_percent_yr >= -0.5 ~ "no change")) %>% 
  mutate(ocean_basin = case_when(stream %in% arctic ~"arctic",
                                 stream %in% n_atl ~ "north atlantic",
                                 stream %in% north_sea_streams | LTER %in% north_sea_lter ~ "north sea",
                                 stream %in% baltic_streams| LTER %in% baltic_lter ~ "baltic",
                                 stream %in% s_atl ~ "south atlantic",
                                 stream %in% indian ~ "indian",
                                 stream %in% n_pac | stream %in% s_pac ~ "pacific",
                                 .default = NA
  ))

trend_flux <- trend_dat_v1 %>% 
  filter(Metric == "Flux") %>% 
  filter(!is.na(change_10_3kg_yr)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(stream,chemical) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=10) %>% 
  mutate(trend_direction = case_when(slope_10_3kg_yr_yr > 5 ~ "increase",
                                     slope_10_3kg_yr_yr < -5 ~ "decrease",
                                     slope_10_3kg_yr_yr  <= 5 & slope_percent_yr >= -5 ~ "no change")) %>% 
  mutate(ocean_basin = case_when(stream %in% arctic ~"arctic",
                                 stream %in% n_atl ~ "north atlantic",
                                 stream %in% north_sea_streams | LTER %in% north_sea_lter ~ "north sea",
                                 stream %in% baltic_streams| LTER %in% baltic_lter ~ "baltic",
                                 stream %in% s_atl ~ "south atlantic",
                                 stream %in% indian ~ "indian",
                                 stream %in% n_pac | stream %in% s_pac ~ "pacific",
                                 .default = NA
  ))


## CONCENTRATION
# no significance associated with these trends, so used arbitrary change/no change threshold for now
fig1 <- trend_conc %>% 
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NOx", 
                              .default = chemical)) %>% 
  group_by(chemical, trend_direction) %>%
  summarise(n=n()) %>% 
  ggplot(aes(chemical,n,fill=trend_direction))+
  geom_col(position="stack")

fig1

## FLUX
# no "percent change" estimate for flux, so used slope value and arbitrary thresholds for increase/decrease
fig2 <- trend_flux %>% 
  filter(!is.na(trend_direction)) %>% 
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NOx", 
                              .default = chemical)) %>% 
  group_by(chemical, trend_direction) %>%
  summarise(n=n()) %>% 
  ggplot(aes(chemical,n,fill=trend_direction))+
  geom_col(position="stack")

fig2

plot_grid(fig1,fig2,nrow=1)


fig3 <- trend_conc %>%
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NOx", 
                              .default = chemical)) %>% 
  # congo percentage is way off
  filter(slope_percent_yr <20) %>% 
  ggplot(aes(ocean_basin, slope_percent_yr))+
  geom_boxplot()+
  geom_jitter(width = 0.2, color = "blue")+
  geom_hline(yintercept = 0, lty=2)+
  facet_wrap(~chemical, nrow=4)

fig4 <- trend_flux %>%
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NOx", 
                              .default = chemical)) %>% 
  # congo percentage is way off
  filter(slope_10_3kg_yr_yr < 10000) %>% 
  ggplot(aes(ocean_basin, slope_10_3kg_yr_yr))+
  geom_boxplot()+
  geom_jitter(width = 0.2, color = "blue")+
  geom_hline(yintercept = 0, lty=2)+
  facet_wrap(~chemical, nrow=4,scales="free_y")


plot_grid(fig3,fig4,ncol=2)


# summarize magnitude of trends -------------------------------------------

si_trends <- trend_dat_v1 %>% 
  rename(Stream_Name = stream) %>% 
  filter(chemical == "DSi") %>% 
  filter(Metric == "Concentration") %>% 
  filter(!is.na(change_mg_L)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(Stream_Name) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=10) %>% 
  arrange(slope_percent_yr) %>% 
  left_join(spatial_dat, by="Stream_Name")

si_trends %>% 
  ggplot(aes(major_land, change_percent))+
  geom_boxplot()


p_trends <- trend_dat_v1 %>% 
  rename(Stream_Name = stream) %>% 
  filter(chemical == "P") %>% 
  filter(Metric == "Concentration") %>% 
  filter(!is.na(change_mg_L)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(Stream_Name) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=10) %>% 
  arrange(slope_percent_yr) %>% 
  left_join(spatial_dat, by="Stream_Name")

p_trends %>% 
  ggplot(aes(major_land, slope_percent_yr))+
  geom_boxplot()+
  geom_jitter()

n_trends <- trend_dat_v1 %>% 
  rename(Stream_Name = stream) %>% 
  mutate(chemical = case_when(chemical == "NO3"|chemical == "NOx" ~ "NO3",
                              .default = chemical)) %>% 
  filter(chemical == "NO3") %>% 
  filter(Metric == "Concentration") %>% 
  filter(!is.na(change_mg_L)) %>% 
  mutate(duration = (Year2-Year1)+1) %>% 
  group_by(Stream_Name) %>%
  filter(duration == max(duration, na.rm = TRUE)) %>%
  ungroup() %>% 
  filter(duration >=10) %>% 
  arrange(slope_percent_yr) %>% 
  left_join(spatial_dat, by="Stream_Name")

n_trends %>% 
  ggplot(aes(major_land, slope_percent_yr))+
  geom_boxplot()+
  geom_jitter()



# Analyze Seasonal trends ---------------------------------------------------------


# Si:N ratios of concentration
ratio_trend_dat <- month_dat %>% 
  select(LTER,Stream_Name,Year,Month,chemical,FNConc_uM) %>% 
  filter(Stream_Name %in% long$Stream_Name) %>% 
  filter(chemical == "NO3"|chemical == "DSi"|chemical == "NOx"|chemical == "P"|chemical == "Si:P") %>% 
  ## Calculate DIN (DIN = NOx <or> NO3 + NH4)
  # Handle "duplicate" values for sites that break across a year so have two values for one year
  dplyr::group_by(LTER, Stream_Name, chemical, Month,Year) %>%
  dplyr::summarize(response_values = mean(FNConc_uM, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(chemical = dplyr::case_when(
    ### NOx is preferred for calculating DIN because it is NO3 + NOx
    chemical == "NOx" | chemical == "NO3" ~ "NO3x",
    .default = chemical)) %>% 
  tidyr::pivot_wider(names_from = chemical,
                     values_from = response_values,
                     values_fn = mean) %>% 
  ## Calculate ratios
  dplyr::mutate(Si_NO3x = ifelse(test = (!is.na(DSi) & !is.na(NO3x)),
                                 yes = (DSi / NO3x), no = NA)) %>% 
  dplyr::mutate(Si_P = ifelse(test = (!is.na(DSi) & !is.na(P)),
                              yes = (DSi / P), no = NA)) %>% 
  pivot_longer(cols = DSi:Si_P,names_to = "chemical", values_to = "FNConc_uM")

## Concentration
conc_slope <- month_dat %>% 
  filter(Stream_Name %in% long$Stream_Name) %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical,Month) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM)) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))

conc_ratio_slope <- ratio_trend_dat %>% 
  filter(Stream_Name %in% long$Stream_Name) %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical,Month) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM)) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))  

yield_slope <- month_dat %>% 
  filter(Stream_Name %in% long$Stream_Name) %>% 
  filter(!is.na(FNYield_kmol_day_km2)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical,Month) %>%
  group_modify(~ modified_sens.slope(.x$FNYield_kmol_day_km2)) %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  mutate(land_group = case_when(land_Impervious < 3 & land_Cropland <65 ~ "less impacted",
                                land_Impervious > 3 & land_Cropland <65 ~ "urban", 
                                land_Impervious < 3 & land_Cropland >65 ~ "agricultural",
                                land_Impervious > 3 & land_Cropland > 65 ~ "mixed", 
                                .default = NA))

# plot by month

conc_slope %>%
  filter(estimates < 5 & estimates >-5) %>% 
  filter(p.value <=0.05) %>%
  filter(land_group != "NA") %>% 
  filter(chemical != "NH4" & chemical != "Si:P") %>%
  filter(chemical == "NO3") %>% 
  ggplot(aes(as.factor(Month),estimates))+
  geom_boxplot()+
  geom_hline(yintercept = 0, lty=2)+
  #geom_jitter(alpha=0.5)+
  facet_wrap(~land_group, scales="free")

conc_ratio_slope %>%
  filter(estimates < 0.1 & estimates > -1) %>% 
  filter(p.value <=0.05) %>% 
  filter(land_group != "NA") %>% 
  filter(chemical == "Si_NO3x") %>%
  ggplot(aes(as.factor(Month),estimates))+
  geom_boxplot()+
  geom_hline(yintercept = 0, lty=2)+
  #geom_jitter(alpha=0.5)+
  facet_wrap(~land_group, scale="free")



# OLD BASEMENT ------------------------------------------------------------





# Plot time series and save --------------------------------------------------------
# plot change over time

plot_dat <- kalman_dat %>%
  filter(Stream_Name %in% long_records$Stream_Name) %>% 
  filter(chemical == "DSi")

## save plot 
p = ggplot(data = plot_dat, aes(x = Year, y = Discharge_cms)) + 
  geom_point()+
  geom_smooth(method="loess")

plots = plot_dat %>%
  group_by(Stream_Name) %>%
  do(plots = p %+% . + facet_wrap(~Stream_Name))

setwd("//home/jankowski/data/plots")
pdf()
plots$plots
dev.off()



# group datasets by decade ------------------------------------------------

# plot by year
annual_dat_v0 %>% 
  filter(chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Year,chemical,Stream_Name) %>% 
  summarise(n = n()) %>% 
  #filter(chemical == "DSi") %>% 
  ggplot(aes(Year,n))+
  geom_col()+
  facet_wrap(~chemical, ncol=1)

# Add decade tag
annual_dat_v1 <- annual_dat_v0 %>% 
  mutate(decade = case_when(Year >=1960 & Year < 1970 ~ "1960",
                            Year >=1970 & Year < 1980 ~ "1970",
                            Year >=1980 & Year < 1990 ~ "1980",
                            Year >=1990 & Year < 2000 ~ "1990",
                            Year >=2000 & Year < 2010 ~ "2000",
                            Year >=2010 & Year < 2020 ~ "2010",
                            Year >=2020 & Year < 2030 ~ "2020",
                            .default = NA)) %>% 
  unite(c(Stream_Name,decade,chemical), col="stream_decade_chemical",sep= "__",remove=FALSE)



# count how many years of data for each stream in each decade for each chemical
decades <- annual_dat_v1 %>% 
  group_by(Stream_Name,chemical,decade) %>% 
  summarise(n=n())

# filter to streams with 7 years of data in a decade for each chemical
decades_v2 <- decades %>% 
  filter(n>=7)

si_decades <- decades_v2 %>% 
  filter(chemical == "DSi") %>% 
  select(Stream_Name,decade,chemical,n) %>% 
  unite(Stream_Name:decade, col="stream_ID",sep= "__",remove=FALSE)

all_chem_decades <- decades_v2 %>% 
  select(Stream_Name,decade,chemical,n) %>% 
  unite(Stream_Name:chemical, col="stream_decade_chemical",sep= "__",remove=FALSE)


decades_v3 <- decades_v2 %>% 
  group_by(decade,chemical) %>% 
  summarise(n=n())

decades_v3 %>% 
  filter(chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  ggplot(aes(decade,n))+
  geom_col()+
  facet_wrap(~chemical, ncol=1)


# analyze trends by decade ------------------------------------------------

conc_slope_decade <- annual_dat_v1 %>% 
  filter(chemical != "DIN") %>% 
  filter(stream_decade_chemical %in% all_chem_decades$stream_decade_chemical) %>% 
  filter(!is.na(FNConc_uM)) %>% 
  group_by(Stream_Name,decade,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM))

yield_slope_decade <- annual_dat_v1 %>% 
  filter(chemical != "DIN") %>% 
  filter(stream_decade_chemical %in% all_chem_decades$stream_decade_chemical) %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  group_by(Stream_Name,decade,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNYield_10_6kmol_yr_km2))

dis_slope_decade <- annual_dat_v1 %>% 
  filter(chemical == "DSi") %>% 
  filter(stream_decade_chemical %in% all_chem_decades$stream_decade_chemical) %>% 
  group_by(Stream_Name,decade) %>%
  group_modify(~ modified_sens.slope(.x$Discharge_cms))

## need to adjust Si:N to be Si:NO3 because of DIN issue
ratio_trend_dat <- annual_dat_v1 %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(Stream_Name != "BARWON RIVER AT MUNGINDI") %>% 
  select(LTER,Stream_Name,Year,decade, chemical,FNConc_uM) %>% 
  filter(chemical == "NO3"|chemical == "DSi"|chemical == "NOx") %>% 
  ## Calculate DIN (DIN = NOx <or> NO3 + NH4)
  # Handle "duplicate" values for sites that break across a year so have two values for one year
  dplyr::group_by(LTER, Stream_Name, chemical, Year) %>%
  dplyr::summarize(response_values = mean(FNConc_uM, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(chemical = dplyr::case_when(
    ### NOx is preferred for calculating DIN because it is NO3 + NOx
    chemical == "NOx" | chemical == "NO3" ~ "NO3x",
    .default = chemical)) %>% 
  tidyr::pivot_wider(names_from = chemical,
                     values_from = response_values,
                     values_fn = mean) %>% 
  ## Calculate ratios
  dplyr::mutate(Si_NO3x = ifelse(test = (!is.na(DSi) & !is.na(NO3x)),
                                 yes = (DSi / NO3x), no = NA)) %>% 
  pivot_longer(cols = DSi:Si_NO3x,names_to = "chemical", values_to = "FNConc_uM")


conc_ratio_slope <- ratio_trend_dat %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM))

sig_conc_slope <- conc_ratio_slope %>% 
  filter(p.value <= 0.05)

conc_ratio_mk <- ratio_trend_dat %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNConc_uM))

slopes <- conc_slope %>% 
  bind_rows(conc_ratio_slope)

write.csv(x = slopes, row.names = F, na = '',
          file = "sen_slopes_conc.csv")


# plot slope by decade ----------------------------------------------------

plot_slope <- conc_slope_decade|> 
  left_join(cluster_dat,by="Stream_Name") |> 
  filter(!is.na(cluster_name)) |> 
  filter(chemical != "Si:DIN"  & chemical != "DIN") |> 
  filter(estimates < 1000 & estimates > -1000) |> 
  filter(cluster_name != "grassland Australia") |> 
  filter(p.value<=0.05) |> 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) |>  
  ggplot(aes(cluster_name,estimates,col=cluster_name))+
  geom_boxplot()+
  geom_hline(yintercept=0, lty=1)+
  geom_jitter(position="jitter")+
  facet_grid(chemical~decade,scales="free_y")+
  theme(axis.text.x = element_text(angle=90))

plot_slope

# ridge plot
si_ridges <- conc_slope_decade |>
  left_join(cluster_dat,by="Stream_Name") |> 
  filter(!is.na(cluster_name)) |> 
  #filter(cluster_name != "grassland Australia") %>% 
  filter(chemical != "Si:DIN"  & chemical != "DIN") |> 
  filter(estimates <1000 &estimates>-80) %>% 
  filter(p.value<=0.05) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical == "NO3") %>% 
  # ordered so that productivity variables in column 1 and others in column 2
  #mutate(chemical = fct_relevel(chemical, "DSi","NOx","P","Si:P","Si:DIN")) |> 
  #mutate(cluster=fct_relevel(cluster,"2","5","6","3","1","4")) |> 
  ggplot(aes(estimates,y=decade,fill=decade))+
  geom_density_ridges(
    aes(
      point_color=decade,
      point_fill=decade),
    #point_shape=21,
    #vline_color=cluster),
    alpha=0.2,
    point_alpha=0.5, 
    point_size=1.5,
    jittered_points=TRUE,
    scale=1.5)+
  geom_vline(xintercept=0)+
  #scale_fill_manual(values=yel_green, guide="none")+
  #scale_color_manual(values=yel_green)+
  #scale_discrete_manual(aesthetics = c("point_fill","point_color"), values=yel_green)+
  facet_wrap(~cluster_name, scales="free")+
  theme_minimal()+
  theme(legend.position="none")+
  xlab("Sen Slope")+ylab("")

si_ridges

# summarize change by decade and chemical

conc_slope_decade_v0 <- conc_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  select(Stream_Name,decade,chemical,cluster_name,change,p.value,statistic,estimates,low.conf,high.conf)


conc_slope_decade_v1 <- conc_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>%
  filter(!is.na(cluster)) %>% 
  filter(chemical != "DIN") %>% 
  group_by(chemical, decade, cluster_name,change) %>% 
  summarise(n=n())

conc_slope_decade_v1 %>% 
  #filter(!is.na(cluster)) %>% 
  filter(chemical == "DSi") %>% 
  #filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(decade, n, fill=change))+
  geom_col(position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  facet_grid(chemical~cluster_name)+
  theme(axis.text.x = element_text(angle=90))

# create dataset to show proportions by decade
conc.change = aggregate(Stream_Name~chemical*decade*cluster_name*change, FUN=length, data=conc_slope_decade_v0)
conc.total = aggregate(Stream_Name~chemical*decade*cluster_name, FUN=length, data=conc_slope_decade_v0)

# not sure why getting an error
conc.prop <- conc.change %>% 
  filter(chemical == "NO3"|chemical == "P"|chemical == "DSi") %>% 
  left_join(conc.total, by=c("chemical","decade","cluster_name")) %>% 
  mutate(prop.streams = Stream_Name.x/Stream_Name.y)

conc.prop %>% 
  filter(chemical == "DSi") %>% 
  #filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(x = decade, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Concentration")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right",
        axis.text.x = element_text(angle = 90))+
  facet_grid(chemical~cluster_name)


#######################################
# Yield change by decade

yield_slope_decade_v0 <- yield_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  select(Stream_Name,decade,chemical,cluster_name,change,p.value,statistic,estimates,low.conf,high.conf)


yield_slope_decade_v1 <- yield_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>%
  filter(!is.na(cluster)) %>% 
  filter(chemical != "DIN") %>% 
  group_by(chemical, decade, cluster_name,change) %>% 
  summarise(n=n())

yield_slope_decade_v1 %>% 
  #filter(!is.na(cluster)) %>% 
  filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(decade, n, fill=change))+
  geom_col(position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  facet_grid(chemical~cluster_name)

# create dataset to show proportions by decade
yield.change = aggregate(Stream_Name~chemical*decade*cluster_name*change, FUN=length, data=yield_slope_decade_v0)
yield.total = aggregate(Stream_Name~chemical*decade*cluster_name, FUN=length, data=yield_slope_decade_v0)

# not sure why getting an error
yield.prop <- yield.change %>% 
  filter(chemical == "NO3"|chemical == "P"|chemical == "DSi") %>% 
  left_join(yield.total, by=c("chemical","decade","cluster_name")) %>% 
  mutate(prop.streams = Stream_Name.x/Stream_Name.y)

yield.prop %>% 
  filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(x = decade, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Yield")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  facet_grid(chemical~cluster_name)

########################################
## Discharge change by decade 
dis_slope_decade_v0 <- dis_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  select(Stream_Name,decade,cluster_name,change,p.value,statistic,estimates,low.conf,high.conf)


dis_slope_decade_v1 <- dis_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  filter(!is.na(cluster)) %>% 
  group_by(decade, cluster_name,change) %>% 
  summarise(n=n())

dis_slope_decade_v1 %>% 
  #filter(!is.na(cluster)) %>% 
  #filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(decade, n, fill=change))+
  geom_col(position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  facet_grid(~cluster_name)



