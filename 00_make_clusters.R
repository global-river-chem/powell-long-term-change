## spatial clustering of sites

library(tidyverse)
library(vegan)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(ggridges)
library(cowplot)
library(ggforce)
library(googledrive)

source(file.path("tools", "biostats.R"))

# spatial data
spatial_link<-"https://drive.google.com/file/d/1RlM0kn-TWeOkStIz8BMGqF6NorA2zHQw/view?usp=sharing"
spatial_folder = drive_get(as_id(spatial_link))
spatial<-drive_download(spatial_folder$drive_resource, overwrite = T)
spatial <- drive_download(file = spatial_folder$id, path = file.path(path,"Spatial_data_average_v3.csv"), overwrite = T)
spatial_dat<-read_csv(spatial$local_path)

# Model for all sites -----------------------------------------------------


# Prep data 
all <- spatial_dat |> 
  select(Stream_Name,Research_Network, climate_zone_name,major_land,watershed_slope,watershed_elevation_median,watershed_elevation_max,
         drainage_area,temperature,precipitation,evapotranspiration,npp,greenup,permafrost,
         snow_prop_area,Max_Daylength,soil_inceptisols,soil_andisols,soil_aridisols,soil_entisols,
         soil_histosol,soil_mollisols,soil_oxisol,soil_spodosols,soil_ultisols,soil_vertisols,
         rocks_volcanic,rocks_sedimentary,rocks_carbonate_evaporite,rocks_metamorphic,rocks_plutonic,
         land_Cropland,land_Forest,land_Grassland_Shrubland,land_Impervious,land_Water,
         land_Ice_Snow,land_Wetland_Marsh)


all_complete <- as.data.frame(all[complete.cases(all),])


all_v1 <- all_complete |> 
  column_to_rownames(var="Stream_Name") 

# log the values - PCA assumes normal distributions so safe to log to make sure that assumption is met

log_plus <- function(x) {
  # Example: subtract the mean of the column
  log(x+1)
}

all_log <- all_v1 |> 
  select(-temperature,-Research_Network,-climate_zone_name,-major_land) |> 
  mutate(across(everything(),log_plus))

all_log_v1 <- bind_cols(all_log, all_v1$temperature) |> 
  rename(temperature = ...34)



# Run PCA to reduce variables for clustering ------------------------------

# PCA analysis
spatial_pca  <- prcomp(all_log_v1, scale=TRUE)
summary(spatial_pca) # proportion variance explained by axes
pca.structure(spatial_pca, all_log_v1, dim=5, cutoff=0.5) # correlations of variables with axes

pca_results_v1 <- as.data.frame(pca.structure(spatial_pca, all_log_v1, dim=5, cutoff=0.5)) # correlations of variables with axes

all_log_v2 <- all_log_v1 |> 
  select(watershed_elevation_median,watershed_slope,
         evapotranspiration,npp,greenup,snow_prop_area,
         soil_oxisol,
         rocks_plutonic,rocks_volcanic,
         land_Grassland_Shrubland,land_Wetland_Marsh,
         land_Cropland,land_Forest,
         temperature)


spatial_pca_v2  <- prcomp(all_log_v2, scale=TRUE)
summary(spatial_pca_v2) # proportion variance explained by axes
pca.structure(spatial_pca_v2, all_log_v2, dim=5, cutoff=0.5) # correlations of variables with axes

pca_results_v2 <- as.data.frame(pca.structure(spatial_pca_v2, all_log_v2, dim=5, cutoff=0.5)) # correlations of variables with axes

# final variables for cluster analysis - combo of Sidney's paper and variables with r>0.6
final_vars <- all_log_v2 |> 
  select(watershed_elevation_median,watershed_slope,
         evapotranspiration,npp,greenup,snow_prop_area,
         soil_oxisol,
         rocks_plutonic,rocks_volcanic,
         land_Grassland_Shrubland,
         land_Cropland,land_Forest,
         temperature)


# save PCA tables
write.csv(pca_results_v1, "PCA_all_variables.csv")
write.csv(pca_results_v2, "PCA_final_variables.csv")


## Plot PCA
st=as.data.frame(spatial_pca_v2$x[,1:3])
st$LTER<-as.factor(all_complete$Research_Network)
st$climate <- as.factor(all_complete$climate_zone_name)
st$major_land <- as.factor(all_complete$major_land)
st$cluster <- as.factor(cluster_info$cluster)

# pull out vectors - for PCA from pca$rotation object
sp=as.data.frame(spatial_pca_v2$rotation[,1:3])
#yz=as.data.frame(rda_sum$biplot[,1:2])

#plot
#setwd("C:/Users/kjankowski/OneDrive - DOI/Documents/Manuscripts/SiSynthesis/Figures")
#pdf("Figure 3b RDA.pdf", height=6,width=8, family="Times")

climate_plot <- ggplot() +
  #geom_mark_ellipse(st, mapping = aes(PC1,PC2, fill=climate, col=climate))+
  geom_point(data = st,aes(PC1,PC2,color=climate),size=4)+
  theme_bw()+
  theme_classic()+
  theme(legend.position = "right")+
  ylim(-4,5)+xlim(-8,6)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_segment(data = sp,aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "open"),linetype=1, linewidth=0.2,colour = "black")+
  geom_text_repel(data = sp,aes(PC1*8,PC2*8,label=row.names(sp)),colour="black")+
  #geom_segment(data = yz,aes(x = 0, y = 0, xend = PC1, yend = PC3, lty=sig_col), 
  #            arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
  #                         type = "closed"), size=0.6)+
  #geom_text_repel(data = yz,aes(PC1,PC3,label=rowname), size=5)+
  xlab("PC1, 44% variation")+ylab("PC2, 23% variation")+
  scale_color_viridis(discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)

climate_plot

land_plot <- ggplot() +
  #geom_mark_ellipse(st, mapping = aes(PC1,PC2, fill=climate, col=climate))+
  geom_point(data = st,aes(PC1,PC2,color=major_land),size=4)+
  theme_bw()+
  theme_classic()+
  theme(legend.position = "right")+
  ylim(-4,5)+xlim(-8,6)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_segment(data = sp,aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "open"),linetype=1, linewidth=0.2,colour = "black")+
  geom_text_repel(data = sp,aes(PC1*8,PC2*8,label=row.names(sp)),colour="black")+
  #geom_segment(data = yz,aes(x = 0, y = 0, xend = PC1, yend = PC3, lty=sig_col), 
  #            arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
  #                         type = "closed"), size=0.6)+
  #geom_text_repel(data = yz,aes(PC1,PC3,label=rowname), size=5)+
  xlab("PC1, 44% variation")+ylab("PC2, 23% variation")+
  scale_color_viridis(discrete=TRUE, option="magma")+
  scale_fill_viridis(discrete=TRUE,option="magma")

land_plot

cluster_plot <- ggplot() +
  #geom_mark_ellipse(st, mapping = aes(PC1,PC2, fill=climate, col=climate))+
  geom_point(data = st,aes(PC1,PC2,color=cluster),size=4)+
  theme_bw()+
  theme_classic()+
  theme(legend.position = "right")+
  ylim(-4,5)+xlim(-8,6)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_segment(data = sp,aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "open"),linetype=1, linewidth=0.2,colour = "black")+
  geom_text_repel(data = sp,aes(PC1*8,PC2*8,label=row.names(sp)),colour="black")+
  #geom_segment(data = yz,aes(x = 0, y = 0, xend = PC1, yend = PC3, lty=sig_col), 
  #            arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
  #                         type = "closed"), size=0.6)+
  #geom_text_repel(data = yz,aes(PC1,PC3,label=rowname), size=5)+
  xlab("PC1, 44% variation")+ylab("PC2, 23% variation")+
  scale_color_viridis(discrete=TRUE, option="magma")+
  scale_fill_viridis(discrete=TRUE,option="magma")

cluster_plot

sites <- plot_grid(climate_plot, land_plot, nrow=1)

ggsave(plot=sites,filename="site_pca.png", height=6,width=14,units="in")

# Evaluate number of clusters ------------------------------------------------------

# have to assign a unique ID ("uniq_id") to each observation to ensure it joins correctly with 
# cluster assignment after the clustering analysis is performed

data_trans_all <-data.stand(final_vars, method="standardize",
                            margin="column",plot=F)

######### evaluate number of clusters
# calculate euclidean distance
data_eucd_all<-vegdist(data_trans_all,method="euclidean",
                       na.rm=TRUE)

# Dissimilarity and silhouette width
nhclus.scree(data_eucd_all,max.k=20)

# Calinski criterion
clust_kmeans_all_compare <- cascadeKM(data_trans_all,inf.gr=3,sup.gr=13,iter=1000,criterion = "calinski")
plot(clust_kmeans_all_compare)

# randomization loop 
# number of clusters to test
no.clus <- 12
# number of random datasets to test
no.ran <- 999
# number of variables in the dataset
no.var <- ncol(data_trans_all)

# matrix to store within SSE values frome the kmean clustering on the original data (first row)
# and randomized data (rest of the rows) where each column is the number of clusters
sse <- matrix(,nrow=no.ran+1,ncol=no.clus-1)

# empty matrix to store randomized dataset
randata <- matrix(,nrow=nrow(data_trans_all),ncol=ncol(data_trans_all))


for(j in 2:no.clus){
  # performing k-means using original data and storing in the first row of SSE
  all_kmeans <- kmeans(data_trans_all, centers=j,iter.max=1000,nstart=5)
  sse[1,j-1] <- sum(all_kmeans$withinss) 
  
  print(j)
  
  # for loop that repeats the clustering a total of k times (i.e. no.ran)
  for(k in 1:no.ran){
    
    # randomly select variables from each column of the dataset
    for (c in 1:no.var){
      randata[,c] <- data_trans_all[,c][sample(nrow(data_trans_all))]
    }
    
    # perform k-means on random data and store in the SSE matrix
    ran.kmeans <- kmeans(randata,centers=j,iter.max=1000,nstart=5)
    sse[k+1,j-1] <- sum(ran.kmeans$withinss)
    
    print(k)
  }
}

diff <- colMeans(sse)-sse[1,]
plot(seq(2,no.clus,by=1),diff,xlab="No. clusters",ylab="mean(ranSSE)-actualSSE",pch=16)



# Run clustering tests ----------------------------------------------------
set.seed(13)
clust_kmeans_all <-kmeans(data_trans_all,centers=8,iter.max=10000,nstart=25)

# get cluster assignment for each observation
clusters <- as.data.frame(clust_kmeans_all$cluster) |> 
  rownames_to_column(var="Stream_Name") |> 
  rename(cluster=`clust_kmeans_all$cluster`)

clusters |> group_by(cluster) |> summarise(n=n())

# join with original dataset by line number "uniq_id"
cluster_info <- clusters|> 
  left_join(all_complete,by="Stream_Name")

cluster_info_v2 <- cluster_info |> 
  mutate(cluster_name = case_when(cluster == 1 ~ "cropland",
                                  cluster == 2 ~ "cropland grassland",
                                  cluster == 3 ~ "forested steep volcanic",
                                  cluster == 4 ~ "snowy high latitude",
                                  cluster == 5 ~ "snowy high elevation",
                                  cluster == 6 ~ "grassland Australia",
                                  cluster == 7 ~ "forested warm flat",
                                  cluster == 8 ~ "tropical forested"))

## read out results
write.csv(cluster_info_v2, "Si_sites_clusters_eight.csv", row.names=FALSE)

cluster_filtering <- cluster_info |> 
  filter(cluster == 5)


# Plotting ----------------------------------------------------------------


# cluster ridge and panel plots ------------------------------------------------------

# plot of actual values
plot_means <- cluster_info |> 
  select(Stream_Name,cluster,watershed_elevation_median,watershed_slope,
         evapotranspiration,npp,greenup,snow_prop_area,
         soil_oxisol,rocks_plutonic,rocks_volcanic,
         land_Cropland,land_Forest,land_Grassland_Shrubland,temperature) |>
  pivot_longer(watershed_elevation_median:temperature,names_to="variable",values_to="value") 

plot_means$variable <- as.factor(plot_means$variable)


plot_means |> 
  ggplot(aes(as.factor(cluster),value))+
  geom_boxplot()+
  facet_wrap(~variable, scales="free")

# Ridge plots
ridges_crop <- plot_means |> 
  #mutate(variable = fct_relevel(variable, "CHL","Turb","BGA","NO3","DO","Cond",
  #                             "pH","Temp","fDOM")) |> 
  #mutate(cluster=fct_relevel(cluster,"2","5","6","3","1","4")) |> 
  filter(variable == "land_Cropland") |> 
  ggplot(aes(value,y=cluster,fill=cluster))+
  geom_density_ridges()+
  scale_fill_manual(values=yel_green)+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.y=element_blank())+
  xlim(0,100)+xlab("Chl (ug/L)")

ridges_crop


ridges_all <- plot_means |> 
  # ordered so that productivity variables in column 1 and others in column 2
  #mutate(variable = fct_relevel(variable, "CHL","fDOM","BGA","Turb",
  #                             "DO","NO3","pH","Temp")) |> 
  #mutate(cluster=fct_relevel(cluster,"2","5","6","3","1","4")) |> 
  #filter(variable != "Cond") |> 
  ggplot(aes(value,y=cluster,fill=cluster))+
  geom_density_ridges2(
    aes(
      point_color=cluster,
      point_fill=cluster),
    #point_shape=21,
    #vline_color=cluster),
    alpha=0.2,
    point_alpha=0.5, 
    point_size=0.5,
    jittered_points=TRUE,
    inherit.aes=TRUE)+
  scale_fill_manual(values=yel_green, guide="none")+
  scale_color_manual(values=yel_green)+
  scale_discrete_manual(aesthetics = c("point_fill","point_color"), values=yel_green)+
  facet_wrap(~variable,scales="free", ncol=2)+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.y=element_blank())+
  xlab("")

ridges_all


ggsave(plot = cluster_values, filename = "cluster_values_boxplots.png", width=8,height=12,units="in")
ggsave(plot = ridges_all, filename = "cluster_values_ridges.png", width=10,height=12,units="in")


# Model for UMR -----------------------------------------------------------

# Prep data 
umr <- spatial_dat |> 
  filter(Research_Network == "UMR") |> 
  select(Stream_Name,Research_Network, climate_zone_name,major_land,watershed_slope,watershed_elevation_median,watershed_elevation_max,
         drainage_area,temperature,precipitation,evapotranspiration,npp,greenup,
         snow_prop_area,soil_inceptisols,soil_aridisols,soil_entisols,
         soil_histosol,soil_mollisols,soil_oxisol,soil_spodosols,soil_ultisols,soil_vertisols,
         rocks_volcanic,rocks_sedimentary,rocks_carbonate_evaporite,rocks_metamorphic,rocks_plutonic,
         land_Cropland,land_Forest,land_Grassland_Shrubland,land_Impervious,land_Water,
         land_Wetland_Marsh)


umr_complete <- umr[complete.cases(umr),]

rownames(umr_complete) <- umr_complete$Stream_Name

umr_v1 <- umr_complete |> 
  select(-Stream_Name)

# log the values - PCA assumes normal distributions so safe to log to make sure that assumption is met

log_plus <- function(x) {
  # Example: subtract the mean of the column
  log(x+1)
}

umr_log <- umr_v1 |> 
  select(-temperature,-Research_Network,-climate_zone_name,-major_land) |> 
  mutate(across(everything(),log_plus))

umr_log_v1 <- bind_cols(umr_log, umr_v1$temperature) |> 
  rename(temperature = ...30)


# PCA analysis
spatial_pca  <- prcomp(umr_log_v1, scale=TRUE)
summary(spatial_pca) # proportion variance explained by axes
pca.structure(spatial_pca, umr_log_v1, dim=5, cutoff=0.7) # correlations of variables with axes

umr_log_v2 <- umr_log_v1 |> 
  select(watershed_slope, watershed_elevation_max,drainage_area,precipitation,
         npp,greenup,snow_prop_area, soil_aridisols,soil_entisols,soil_histosol,
         soil_mollisols,soil_spodosols,soil_vertisols,rocks_volcanic,rocks_plutonic,
         rocks_metamorphic,land_Water,
         land_Cropland,land_Forest,temperature)


spatial_pca_v2  <- prcomp(umr_log_v2, scale=TRUE)
summary(spatial_pca_v2) # proportion variance explained by axes
pca.structure(spatial_pca_v2, umr_log_v2, dim=5, cutoff=0.7) # correlations of variables with axes

biplot(spatial_pca_v2,pch=16)

## Plot PCA - do this later
st=as.data.frame(spatial_pca_v2$x[,1:3])
st$LTER<-as.factor(all_complete$Research_Network)
st$climate <- as.factor(all_complete$climate_zone_name)
st$major_land <- as.factor(all_complete$major_land)

# pull out vectors - for PCA from pca$rotation object
sp=as.data.frame(spatial_pca_v2$rotation[,1:3])
#yz=as.data.frame(rda_sum$biplot[,1:2])


climate_plot <- ggplot() +
  #geom_mark_ellipse(st, mapping = aes(PC1,PC2, fill=climate, col=climate))+
  geom_point(data = st,aes(PC1,PC2,color=climate),size=4)+
  theme_bw()+
  theme_classic()+
  theme(legend.position = "right")+
  ylim(-4,5)+xlim(-8,6)+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=12))+
  geom_segment(data = sp,aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
               arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
                             type = "open"),linetype=1, linewidth=0.2,colour = "black")+
  geom_text_repel(data = sp,aes(PC1*8,PC2*8,label=row.names(sp)),colour="black")+
  #geom_segment(data = yz,aes(x = 0, y = 0, xend = PC1, yend = PC3, lty=sig_col), 
  #            arrow = arrow(angle=22.5,length = unit(0.35,"cm"),
  #                         type = "closed"), size=0.6)+
  #geom_text_repel(data = yz,aes(PC1,PC3,label=rowname), size=5)+
  xlab("PC1, 32% variation")+ylab("PC2, 18% variation")+
  scale_color_viridis(discrete=TRUE)+
  scale_fill_viridis(discrete=TRUE)

climate_plot

ggsave(plot=sites,filename="site_pca.png", height=6,width=14,units="in")

# clustering all------------------------------------------------------

# have to assign a unique ID ("uniq_id") to each observation to ensure it joins correctly with 
# cluster assignment after the clustering analysis is performed

data_trans_umr <-data.stand(umr_log_v2, method="standardize",
                            margin="column",plot=F)

######### evaluate number of clusters
# calculate euclidean distance
data_eucd_umr<-vegdist(data_trans_umr,method="euclidean",
                       na.rm=TRUE)

# Dissimilarity and silhouette width
nhclus.scree(data_eucd_umr,max.k=15)

# Calinski criterion
clust_kmeans_umr_compare <- cascadeKM(data_trans_umr,inf.gr=2,sup.gr=10,iter=1000,criterion = "calinski")
plot(clust_kmeans_umr_compare)

# randomization loop 
# number of clusters to test
no.clus <- 10
# number of random datasets to test
no.ran <- 999
# number of variables in the dataset
no.var <- 20

# matrix to store within SSE values frome the kmean clustering on the original data (first row)
# and randomized data (rest of the rows) where each column is the number of clusters
sse <- matrix(,nrow=no.ran+1,ncol=no.clus-1)

# empty matrix to store randomized dataset
randata <- matrix(,nrow=nrow(data_trans_umr),ncol=ncol(data_trans_umr))

for(j in 2:no.clus){
  # performing k-means using original data and storing in the first row of SSE
  umr_kmeans <- kmeans(data_trans_umr, centers=j,iter.max=1000,nstart=5)
  sse[1,j-1] <- sum(umr_kmeans$withinss) 
  
  print(j)
  
  # for loop that repeats the clustering a total of k times (i.e. no.ran)
  for(k in 1:no.ran){
    
    # randomly select variables from each column of the dataset
    for (c in 1:no.var){
      randata[,c] <- data_trans_umr[,c][sample(nrow(data_trans_umr))]
    }
    
    # perform k-means on random data and store in the SSE matrix
    ran.kmeans <- kmeans(randata,centers=j,iter.max=1000,nstart=5)
    sse[k+1,j-1] <- sum(ran.kmeans$withinss)
    
    print(k)
  }
}

diff <- colMeans(sse)-sse[1,]
plot(seq(2,no.clus,by=1),diff,xlab="No. clusters",ylab="mean(ranSSE)-actualSSE",pch=16)


######## Run clustering for umr pools
clust_kmeans_umr <-kmeans(data_trans_umr,centers=5,iter.max=10000,nstart=25)

# get cluster assignment for each observation
clusters <- as.data.frame(clust_kmeans_umr$cluster) |> 
  rownames_to_column(var="Stream_Name") |> 
  rename(cluster=`clust_kmeans_umr$cluster`)

clusters |> group_by(cluster) |> summarise(n=n())

write.csv(clusters,"umr_cluster_assignments.csv", row.names=FALSE)

# join with original dataset by line number "uniq_id"
cluster_info <- clusters|> 
  left_join(umr_complete,by="Stream_Name")

vars <- colnames(data_trans_umr)

cluster_info |> 
  pivot_longer(watershed_slope:land_Wetland_Marsh, names_to="variable",values_to="value") |> 
  filter(variable %in% vars) |> 
  ggplot(aes(as.factor(cluster),value))+
  geom_boxplot()+
  facet_wrap(~variable, scales="free")

