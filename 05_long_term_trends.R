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

# Annual data - Kalman 
kalman_dat <- read.csv(file.path(path_2, "Full_Results_WRTDS_kalman_annual.csv"))

# spatial data
spatial_link<-"https://drive.google.com/file/d/1RlM0kn-TWeOkStIz8BMGqF6NorA2zHQw/view?usp=sharing"

spatial_folder = drive_get(as_id(spatial_link))

spatial<-drive_download(spatial_folder$drive_resource, overwrite = T)

spatial <- drive_download(file = spatial_folder$id, path = file.path(path,"Spatial_data_average_v3.csv"), overwrite = T)

spatial_dat<-read_csv(spatial$local_path)


# cluster data
# spatial data
cluster_link<-"https://drive.google.com/file/d/13_cncYkrEw4ZgSY15TGMjMVivs_sr9d5/view?usp=sharing"

cluster_folder = drive_get(as_id(cluster_link))

cluster<-drive_download(cluster_folder$drive_resource, overwrite = T)

cluster <- drive_download(file = cluster_folder$id, path = file.path(path,"Si_sites_cluster_eight.csv"), overwrite = T)

cluster_dat<-read_csv(cluster$local_path)

cluster_dat$cluster_name <- as.factor(cluster_dat$cluster_name)

glimpse(cluster_dat)


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
duration_dat <- annual_dat %>% 
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
long_record_n <- duration_dat %>% 
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

# sites with both Si and N
si_n <- long_record_si %>% 
  select(Stream_Name) %>% 
  inner_join(long_record_n, by="Stream_Name") 

# sites with all three
si_n_p <- si_n %>% 
  inner_join(long_record_p,by="Stream_Name") %>% 
  #left_join(pre,by="Stream_Name") %>% 
  select(Stream_Name)

# full join to get all sites included
si_n <- long_record_si %>% 
  select(Stream_Name) %>% 
  full_join(long_record_n, by="Stream_Name") 

si_n_p <- si_n %>% 
  full_join(long_record_p,by="Stream_Name") %>% 
  #left_join(pre,by="Stream_Name") %>% 
  select(Stream_Name)


# Plot streams that have Si, NO3/NOx and P by LTER
duration_dat %>% 
  filter(Stream_Name %in% si_n_p$Stream_Name) %>% 
  filter(chemical == "DSi"| chemical == "NO3"|chemical == "NOx"|chemical == "P") %>% 
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


# join with cluster data
long_record_cluster <- long_record_si %>% 
  left_join(cluster_dat,by="Stream_Name")

totals <- long_record_cluster %>% 
  group_by(cluster) %>% 
  summarise(n=n())


# Plot averages by cluster ------------------------------------------------

plot_dat <- kalman_dat %>%
  filter(Stream_Name %in% long_records$Stream_Name) %>% 
  #filter(chemical == "DSi") %>% 
  left_join(cluster_dat, by="Stream_Name")

glimpse(plot_dat)

plot_dat %>% 
  filter(chemical != "DIN") %>% 
  filter(!is.na(cluster)) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  ggplot(aes(as.factor(cluster_name),GenYield, fill=cluster_name)) +
  geom_boxplot()+
  facet_wrap(~chemical, scales="free", nrow=2)+
  coord_flip()+
  xlab("")+ylab("log(Yield (10^6kg/yr/km2))")+
  scale_y_log10()+
  theme(legend.position="none")+
  scale_fill_brewer(palette="Set2")



# Plot time series --------------------------------------------------------
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
kalman_dat %>% 
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
kalman_dat_v1 <- kalman_dat %>% 
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
decades <- kalman_dat_v1 %>% 
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


# Analyze trends ----------------------------------------------------------
conc_slope <- annual_dat %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNConc_mgL)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_mgL)) %>% 
  mutate(variable = rep("FNConc"))

sig_conc_slope <- conc_slope %>% 
  filter(p.value <= 0.05)

conc_mk <- annual_dat %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNConc_mgL)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNConc_mgL)) %>% 
  mutate(variable = rep("FNConc"))

sig_conc_mk <- conc_mk %>% 
  filter(p.value <= 0.05)

yield_slope <- annual_dat %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNYield)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNYield)) %>% 
  mutate(variable = rep("FNYield"))

sig_yield_slope <- yield_slope %>% 
  filter(p.value <= 0.05)

yield_mk <- annual_dat %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNYield)) %>% 
  group_by(Stream_Name,chemical) %>%
  group_modify(~ modified_mk_test(.x$FNYield)) %>% 
  mutate(variable = rep("FNYield"))

sig_yield_mk <- yield_mk %>% 
  filter(p.value <= 0.05)

dis_slope <- annual_dat %>% 
  filter(chemical == "DSi") %>% 
  filter(Stream_Name %in% long_record_si$Stream_Name) %>% 
  filter(!is.na(FNYield)) %>% 
  group_by(Stream_Name) %>%
  group_modify(~ modified_sens.slope(.x$Discharge_cms)) %>% 
  mutate(variable = rep("Discharge"))

sig_dis <- dis_slope %>% 
  filter(p.value <= 0.05)

# analyze trends by decade ------------------------------------------------

# RUN TREND TEST
conc_slope_decade <- kalman_dat_v1 %>% 
  filter(chemical != "DIN") %>% 
  filter(stream_decade_chemical %in% all_chem_decades$stream_decade_chemical) %>% 
  filter(!is.na(FNConc_mgL)) %>% 
  group_by(Stream_Name,decade,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_mgL))

yield_slope_decade <- kalman_dat_v1 %>% 
  filter(chemical != "DIN") %>% 
  filter(stream_decade_chemical %in% all_chem_decades$stream_decade_chemical) %>% 
  filter(!is.na(FNYield)) %>% 
  group_by(Stream_Name,decade,chemical) %>%
  group_modify(~ modified_sens.slope(.x$FNConc_mgL))



conc_slope_decade_v0 <- conc_slope_decade %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  left_join(cluster_dat,by="Stream_Name") %>% 
  select(Stream_Name,decade,chemical,cluster_name,change,p.value,statistic,estimates,low.conf,high.conf)


# summarize change by decade and chemical
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
  filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(decade, n, fill=change))+
  geom_col(position="stack")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  facet_grid(chemical~cluster_name)

# create dataset to show proportions by decade
conc.change = aggregate(Stream_Name~chemical*decade*cluster_name*change, FUN=length, data=conc_slope_decade_v0)
conc.total = aggregate(Stream_Name~chemical*decade*cluster_name, FUN=length, data=conc_slope_decade_v0)

# not sure why getting an error
conc.prop <- conc.change %>% 
  filter(chemical == "NO3"|chemical == "P"|chemical == "DSi") %>% 
  left_join(conc.total, by=c("chemical","decade","cluster_name")) %>% 
  mutate(prop.streams = Stream_Name.x/Stream_Name.y)

conc.prop %>% 
  filter(chemical == "DSi"|chemical == "NO3"|chemical == "P") %>% 
  ggplot(aes(x = decade, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Concentration")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  facet_grid(chemical~cluster_name)




# Plot results of trencluster_name# Plot results of trend tests ---------------------------------------------
si_conc_change <- conc_slope %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change")) %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical))

si_conc_change_lu <- si_conc_change %>% 
  left_join(cluster_dat,by="Stream_Name")

# plot in order of magnitude change# plot in order of magnitude change
si_conc_plot <- si_conc_change_lu %>% 
  filter(!is.na(cluster)) %>% 
  filter(Stream_Name != "BILLABONG CREEK AT DARLOT") %>% 
  filter(chemical == "DSi") %>% 
  filter(p.value <0.05) %>% 
  #filter(estimates <=2) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates,col=as.factor(cluster_name),fill=as.factor(cluster_name)))+
  geom_col()+
  geom_hline(yintercept=0) +
  facet_wrap(~cluster_name, scales="free",nrow=6)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        legend.position="none")+
  ylab("Sen Slope")+xlab("")+
  ggtitle("DSi")

si_conc_plot

NOx_conc_plot <- si_conc_change_lu %>% 
  filter(Stream_Name != "BILLABONG CREEK AT DARLOT") %>% 
  filter(chemical == "NO3"|chemical == "NOx") %>% 
  filter(p.value <=0.05) %>% 
  filter(estimates <=2) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates,col=as.factor(cluster),fill=as.factor(cluster)))+
  geom_col()+
  geom_hline(yintercept=0) +
  facet_wrap(~cluster, scales="free", nrow=6)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        legend.position="none")+
  ylab("Sen Slope")+xlab("")+
  ggtitle("NOx or NO3")

p_conc_plot <- si_conc_change_lu %>% 
  filter(Stream_Name != "BILLABONG CREEK AT DARLOT") %>% 
  filter(chemical == "P") %>% 
  filter(p.value <=0.05) %>% 
  filter(estimates <=2) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates,col=as.factor(cluster),fill=as.factor(cluster)))+
  geom_col()+
  geom_hline(yintercept=0) +
  facet_wrap(~cluster, scales="free", nrow=6)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        legend.position="none")+
  ylab("Sen Slope")+xlab("")+
  ggtitle("P")

SiN_conc_plot <- si_conc_change_lu %>% 
  filter(Stream_Name != "BILLABONG CREEK AT DARLOT") %>% 
  filter(chemical == "Si:DIN") %>% 
  filter(p.value <=0.05) %>% 
  filter(estimates <=2) %>% 
  filter(!is.na(cluster)) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates,col=as.factor(cluster),fill=as.factor(cluster)))+
  geom_col()+
  geom_hline(yintercept=0) +
  facet_wrap(~cluster, scales="free", nrow=6)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        legend.position="none")+
  ylab("Sen Slope")+xlab("")+
  ggtitle("Si:N ratio")

SiP_conc_plot <- si_conc_change_lu %>% 
  filter(Stream_Name != "BILLABONG CREEK AT DARLOT") %>% 
  filter(chemical == "Si:P") %>% 
  filter(p.value <=0.05) %>% 
  filter(estimates <=2) %>%
  filter(!is.na(cluster)) %>% 
  ggplot(aes(reorder(Stream_Name,estimates),estimates,col=as.factor(cluster),fill=as.factor(cluster)))+
  geom_col()+
  geom_hline(yintercept=0) +
  facet_wrap(~cluster, scales="free", nrow=6)+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        legend.position="none")+
  ylab("Sen Slope")+xlab("")+
  ggtitle("Si:P ratio")

nutrient_plot <- plot_grid(si_conc_plot,NOx_conc_plot,p_conc_plot, ncol=3, align="vh")
ratio_plot <- plot_grid(SiN_conc_plot,SiP_conc_plot, ncol=2,align="vh")

nutrient_plot
ratio_plot

## Aggregate by cluster

conc.agr = aggregate(Stream_Name~chemical*cluster_name*change, FUN=length, data=si_conc_change_lu)
conc.agr1 = aggregate(Stream_Name~chemical*cluster_name, FUN=length, data=si_conc_change_lu)
conc.agr2 <- conc.agr %>% 
  left_join(conc.agr1, by=c("chemical","cluster_name"))

conc.agr2$prop.streams <-  conc.agr2$Stream_Name.x/conc.agr2$Stream_Name.y

conc.agr2$chemical <- as.factor(conc.agr2$chemical)

conc.agr2 %>% 
  mutate(chemical = case_when(chemical == "NO3"| chemical == "NOx" ~ "NO3", 
                              .default = chemical)) %>% 
  filter(chemical == "Si:P") %>% 
  #mutate(chemical = fct_relevel(chemical, "DSi","NO3","NOx","P","DIN","NH4","Si:DIN","Si:P")) %>% 
  ggplot(aes(x = cluster_name, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=0), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Concentration")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")+
  ggtitle("Si:P")



si_conc_change %>% 
  filter(p.value<0.05) %>% 
  filter(estimates>-1) %>% 
  filter(estimates<1) %>% 
  ggplot(aes(estimates))+
  geom_histogram()+
  facet_wrap(~chemical)

# x axis - land cover/biome, panels 

## Yield 
si_yield_change <- si_yield %>% 
  #filter(p.value <0.05) %>% 
  mutate(change = case_when(p.value < 0.05 & estimates > 0 ~ "increase",
                            p.value < 0.05 & estimates < 0 ~ "decrease",
                            p.value >= 0.05 ~ "no change"))

si_yield_change_lu <- si_yield_change %>% 
  left_join(spatial_dat_avg,by="Stream_Name") %>% 
  filter(major_land != "Bare") %>% 
  filter(!is.na(major_land))

yield.agr = aggregate(Stream_Name~chemical*major_land*change, FUN=length, data=si_yield_change_lu)
yield.agr1 = aggregate(Stream_Name~chemical*major_land, FUN=length, data=si_yield_change_lu)
yield.agr2 <- yield.agr %>% 
  left_join(yield.agr1, by=c("major_land","chemical"))

yield.agr2$prop.streams <-  yield.agr2$Stream_Name.x/yield.agr2$Stream_Name.y


yield = 
  yield.agr2 %>% 
  filter(chemical == "DIN"|chemical == "P"|chemical == "DSi") %>% 
  ggplot(aes(x = chemical, y=prop.streams,fill=change))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90), legend.title=element_blank())+
  xlab("")+ylab("Proportion of Sites")+ggtitle("Yield")+
  scale_fill_viridis(discrete=TRUE, option="magma")+
  theme(legend.position="right")

yield

library(cowplot)
plot_grid(conc,yield,nrow=2,align="vh",nrow=1,rel_heights=c(1.25,1))

# ridge plot
ridges_all <- si_conc_change_lu |> 
  filter(p.value<=0.05) %>% 
  filter(chemical == "DSi") %>% 
  filter(variable == "FNConc") %>% 
  # ordered so that productivity variables in column 1 and others in column 2
  #mutate(chemical = fct_relevel(chemical, "DSi","NOx","P","Si:P","Si:DIN")) |> 
  #mutate(cluster=fct_relevel(cluster,"2","5","6","3","1","4")) |> 
  ggplot(aes(estimates,y=cluster,fill=cluster))+
  geom_density_ridges(
    aes(
      point_color=cluster,
      point_fill=cluster),
    #point_shape=21,
    #vline_color=cluster),
    alpha=0.2,
    point_alpha=0.5, 
    point_size=0.5,
    jittered_points=TRUE)+
  geom_vline(xintercept=0)+
  #scale_fill_manual(values=yel_green, guide="none")+
  #scale_color_manual(values=yel_green)+
  #scale_discrete_manual(aesthetics = c("point_fill","point_color"), values=yel_green)+
  facet_wrap(~chemical,scales="free", ncol=2)+
  theme_minimal()+
  theme(legend.position="none")+
  xlab("Sen Slope")+ylab("")

ridges_all

