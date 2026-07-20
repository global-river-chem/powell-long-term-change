## looking at long term trends

library(tidyverse)
library(trend)
library(viridis)
library(googledrive)
library(ggridges)
library(cowplot)


# directory for downloading files from google drive
(path <- scicomptools::wd_loc(local = FALSE, remote_path = file.path('/', "home","jankowski","data")))

# directory for accessing WRTDS results
(path_2 <- scicomptools::wd_loc(local = FALSE, remote_path = file.path('/', "home", "shares", "lter-si", "WRTDS","WRTDS Results_2025")))

# Data import -------------------------------------------------------------
# Annual data - generalized flow normalization 
annual_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_annual.csv"))
kalman_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_kalman_annual.csv"))

## spatial data
spatial_link<-"https://drive.google.com/file/d/1abavWUHwPbWpd_PbCHfA-Cf-B8_HYArV/view?usp=drive_link"
spatial_folder = drive_get(as_id(spatial_link))
spatial <- drive_download(file = spatial_folder$id, path = file.path(path,"Spatial_data_average_v3.csv"), overwrite = T)
spatial_dat<-read_csv(spatial$local_path)

spatial_link_2 <- "https://drive.google.com/file/d/1abavWUHwPbWpd_PbCHfA-Cf-B8_HYArV/view?usp=drive_link"
spatial_folder_2 <- drive_get(as_id(spatial_link_2))
spatial_2 <- drive_download(file=spatial_folder_2$id, path=file.path(path,"all_data_si-extract_3_20260629.csv"),overwrite=TRUE)
spatial_dat_2 <- read_csv(spatial_2$local_path)

## cluster data
# spatial data
cluster_link<-"https://drive.google.com/file/d/13_cncYkrEw4ZgSY15TGMjMVivs_sr9d5/view?usp=drive_link"
cluster_folder = drive_get(as_id(cluster_link))
cluster <- drive_download(file = cluster_folder$id, path = file.path(path,"Si_sites_cluster_eight.csv"), overwrite = T)
cluster_dat<-read_csv(cluster$local_path)

cluster_dat$cluster_name <- as.factor(cluster_dat$cluster_name)
glimpse(cluster_dat)


# sites select ------------------------------------------------------------

## Coastal sites
coastal_link <- "https://drive.google.com/file/d/1Xf8ErPzJJJBwR230NSTQrH5tAH4U0Ovb/view?usp=drive_link"
coastal_folder <- drive_get(as_id(coastal_link))
coastal <- drive_download(file=coastal_folder$id, path = file.path(path,"coastal_sites_list.csv"),overwrite = TRUE)
coastal_list <- read_csv(coastal$local_path)

river_list <- unique(coastal_list$Stream_Name)

annual_dat %>% 
  filter(LTER == "Australia"|LTER == "MD") %>% 
  distinct(Stream_Name)

annual_dat_v0 <- annual_dat %>% 
  filter(Stream_Name %in% river_list | LTER == "Sweden"|LTER == "Seine"|Stream_Name == "COLUMBIA RIVER AT PORT WESTWARD")

glimpse(annual_dat_v0)

annual_dat_v0 %>% 
  filter(chemical == "DSi") %>% 
  ggplot(aes(Year,Conc_uM,col=Stream_Name))+
  geom_point()+
  facet_wrap(~LTER)+
  theme(legend.position="none")

# looking at individual LTER data
annual_dat_v0 %>% 
  filter(chemical == "DSi") %>% 
  filter(LTER == "USGS") %>% 
  ggplot(aes(Year,Conc_uM,col=Stream_Name))+
  geom_point()+geom_line()+
  facet_wrap(~LTER)+
  theme(legend.position="bottom")


# Functions needed --------------------------------------------------------

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

# Filter out datasets with long records ----------------------------------

# summarizing record lengths by stream and chemical
duration_dat <- annual_dat_v0 %>% 
  group_by(Stream_Name,chemical,LTER) %>% 
  summarise(min_year = min(Year),
            max_year = max(Year),
            duration = (max_year - min_year)+1)

# create dataset of just long records 
long_records <- duration_dat %>% 
  filter(duration >= 15) 

# plot temporal overlap of chemicals
duration_dat %>% 
  filter(chemical == "DSi"| chemical == "NO3"|chemical == "NOx"|chemical == "P"| chemical == "NH4") %>% 
  #filter(chemical == "Si:P"|chemical == "Si:DIN") %>% 
  ggplot()+
  geom_segment(aes(x=min_year,xend=max_year,
                   y=Stream_Name,xend=Stream_Name))+
  facet_wrap(~chemical)+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank())+
  ylab("Streams")+xlab("Year")

# look at record lengths
long_record_nox <- duration_dat %>% 
  filter(chemical == "NO3"|chemical == "NOx") %>% 
  filter(duration >= 15) %>% 
  select(Stream_Name) 

long_record_nh4 <- duration_dat %>% 
  filter(chemical == "NH4") %>% 
  filter(duration >=15) %>% 
  select(Stream_Name) %>% 
  mutate(chemical = "NH4")

long_record_p <- duration_dat %>% 
  filter(chemical == "P") %>% 
  filter(duration >= 15) %>% 
  select(Stream_Name) %>% 
  mutate(chemical = "P")

long_record_si <- duration_dat %>% 
  filter(chemical == "DSi") %>% 
  filter(duration >= 15) %>% 
  select(Stream_Name) %>% 
  mutate(chemical = "DSi")

# sites with both Si and N, need to have long Si record
si_n <- long_record_si %>% 
  select(Stream_Name) %>% 
  left_join(long_record_nox, by="Stream_Name") 

# full join to get all sites included
si_nox_p <- si_n %>% 
  select(Stream_Name) %>% 
  left_join(long_record_p, by="Stream_Name") 

# all sites with long record of at least one
si_nox_p_nh4 <- si_nox_p %>% 
  left_join(long_record_nh4,by="Stream_Name") %>% 
  #left_join(pre,by="Stream_Name") %>% 
  select(Stream_Name)

# Plot streams that have Si, NO3/NOx and P by LTER
duration_dat %>% 
  filter(Stream_Name %in% si_n_p$Stream_Name) %>% 
  filter(chemical == "NO3"|chemical == "NOx") %>% 
  #filter(chemical == "DSi"| chemical == "NO3"|chemical == "NOx"|chemical == "P") %>% 
  #filter(chemical == "Si:P"|chemical == "Si:DIN") %>% 
  ggplot()+
  geom_segment(aes(x=min_year,xend=max_year,
                   y=Stream_Name,xend=Stream_Name, col=chemical), alpha=0.5)+
  facet_wrap(~LTER, scales="free_y")+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  ylab("Streams")+xlab("Year")

# plot chemicals by LTER
duration_dat %>% 
  filter(Stream_Name %in% long_record_n$Stream_Name) %>% 
  filter(chemical == "NO3"|chemical == "NOx") %>% 
  #filter(chemical == "Si:P"|chemical == "Si:DIN") %>% 
  ggplot()+
  geom_segment(aes(x=min_year,xend=max_year,
                   y=Stream_Name,xend=Stream_Name, col=chemical), alpha=0.5)+
  facet_wrap(~LTER, scales="free_y")+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  ylab("Streams")+xlab("Year")


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


# Analyze trends ----------------------------------------------------------

conc_slope <- annual_dat_v0 %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNConc_uM)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_uM)) %>% 
  mutate(variable = rep("FNConc"))

sig_conc_slope <- conc_slope %>% 
  #filter(chemical == "DSi") %>% 
  filter(p.value <= 0.05)


conc_mk <- annual_dat_v0 %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNConc_uM)) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNConc_uM)) %>% 
  mutate(variable = rep("FNConc"))

sig_conc_mk <- conc_mk %>% 
  filter(p.value <= 0.05)

yield_slope <- annual_dat_v0 %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>%
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNYield_10_6kmol_yr_km2)) %>% 
  mutate(variable = rep("FNYield"))

sig_yield_slope <- yield_slope %>% 
  filter(p.value <= 0.05)

yield_mk <- annual_dat_v0 %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical != "Si:DIN" & chemical != "DIN") %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNYield_10_6kmol_yr_km2)) %>% 
  mutate(variable = rep("FNYield"))

sig_yield_mk <- yield_mk %>% 
  filter(p.value <= 0.05)

# DISCHARGE
dis_slope <- annual_dat_v0 %>% 
  filter(chemical == "DSi") %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(Discharge_cms)) %>% 
  group_by(Stream_Name) %>%
  group_modify(~ modified_sens.slope(.x$Discharge_cms)) %>% 
  mutate(variable = rep("Discharge"))

sig_dis <- dis_slope %>% 
  filter(p.value <= 0.05)


## RATIOS -- CHECK THIS!
## need to adjust Si:N to be Si:NO3 because of DIN issue

# Ratios of CONCENTRATION
ratio_trend_dat <- annual_dat_v0 %>% 
  select(LTER,Stream_Name,Year,chemical,FNConc_uM) %>% 
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


# NOW ADD Si:NO3 to full dataframe of slopes
slopes <- conc_slope %>% 
  bind_rows(conc_ratio_slope)

write.csv(x = slopes, row.names = F, na = '',
          file = "sen_slopes_conc.csv")


# Ratios of YIELD
ratio_trend_dat_yield <- annual_dat %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(Stream_Name != "BARWON RIVER AT MUNGINDI") %>% 
  select(LTER,Stream_Name,Year,chemical,FNYield_10_6kmol_yr_km2) %>% 
  filter(chemical == "NO3"|chemical == "DSi"|chemical == "NOx") %>% 
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
  pivot_longer(cols = DSi:Si_NO3x,names_to = "chemical", values_to = "FNYield_10_6kmol_yr_km2")

yield_ratio_slope <- ratio_trend_dat_yield %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNYield_10_6kmol_yr_km2))

sig_yield_slope <- yield_ratio_slope %>% 
  filter(p.value <= 0.05)

yield_ratio_mk <- ratio_trend_dat_yield %>% 
  filter(!is.na(FNYield_10_6kmol_yr_km2)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNYield_10_6kmol_yr_km2))

# NOW ADD Si:NO3 to full dataframe of slopes
yield_slopes <- yield_slope %>% 
  bind_rows(yield_ratio_slope)

# plot trends for review --------------------------------------------------

sig_conc_slope %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  filter(chemical == "DSi") %>% 
  #filter(Stream_Name %in% sig_conc_slope$Stream_Name) %>% 
  filter(estimates<1000) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates, fill=major_land))+
  geom_col(position="identity")+
  geom_hline(yintercept=0)+
  theme_classic()+
  theme(legend.position = "bottom")+
  coord_flip()+
  #facet_wrap(~major_land,scales="free_x")+
  xlab("Stream")+ylab("Sen Slope value")

sig_yield_slope %>% 
  left_join(spatial_dat, by="Stream_Name") %>% 
  filter(chemical == "DSi") %>% 
  filter(estimates<1000) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates, fill=major_land))+
  geom_col(position="identity")+
  geom_hline(yintercept=0)+
  theme_classic()+
  theme(legend.position = "bottom")+
  coord_flip()+
  #facet_wrap(~major_land,scales="free_x")+
  xlab("Stream")+ylab("Sen Slope value")

# summarize direction of trends -------------------------------------------------

# concentration
si_conc_change <- slopes %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical))

si_conc_change_lu <- si_conc_change %>% 
  left_join(spatial_dat,by="Stream_Name")


# yield
yield_change <- yield_slopes %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical))

yield_change_lu <- yield_change %>% 
  left_join(spatial_dat,by="Stream_Name")


# Discharge
dis_change <- dis_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change"))

dis_cluster <- dis_change %>% 
  left_join(cluster_dat,by="Stream_Name")


## Concentration 
conc.agr = aggregate(Stream_Name~chemical*major_land*change, FUN=length, data=si_conc_change_lu)
conc.agr1 = aggregate(Stream_Name~chemical*major_land, FUN=length, data=si_conc_change_lu)
conc.agr2 <- conc.agr %>% 
  left_join(conc.agr1, by=c("chemical","major_land"))

conc.agr2$prop.streams <-  conc.agr2$Stream_Name.x/conc.agr2$Stream_Name.y

conc.agr2$chemical <- as.factor(conc.agr2$chemical)

conc.agr2 %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical == "DSi") %>% 
  #mutate(chemical = fct_relevel(chemical, "DSi","NO3","NOx","P","DIN","NH4","Si:DIN","Si:P")) %>% 
  ggplot(aes(x = major_land, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Concentration")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  ggtitle("DSi")

## Yield 
yield.agr = aggregate(Stream_Name~chemical*cluster_name*change, FUN=length, data=yield_change_cluster)
yield.agr1 = aggregate(Stream_Name~chemical*cluster_name, FUN=length, data=yield_change_cluster)
yield.agr2 <- yield.agr %>% 
  left_join(yield.agr1, by=c("chemical","cluster_name"))

yield.agr2$prop.streams <-  yield.agr2$Stream_Name.x/yield.agr2$Stream_Name.y

yield.agr2$chemical <- as.factor(yield.agr2$chemical)

yield.agr2 %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical == "Si_NO3x") %>% 
  ggplot(aes(x = cluster_name, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  ggtitle("Yield Si:NO3")

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



