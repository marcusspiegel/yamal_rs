library(raster)
library(rgdal)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(parallel)
library(sf)
library(exactextractr)
library(ggsignif)

setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing")

################################################################################
# Preprocessing rasters into data frames
# # Raster layers
# vsi <- raster('MODIS/withTaz/VSI/Monthly_1km/SensTotalW.tif')
# geomorphons <- raster('ArcticDEM/withTaz/2020-06-11/Geomorphons_withTaz_32m_50cell_0.5deg_flatSplitThresholdExtendVKernelMajCRCRCR_UTM.tif')
# landcover <- raster('Pangaea/Landcover/ZAM_LCP_LANDC_SEN12_V02_20150815_20180830_T01.tif')
# 
# # Vector polygons
# yp <- readOGR(dsn = 'Brigades/yp.shp', layer = 'yp')
# yp$Id <- 'yp'
# taz <- readOGR(dsn = 'Brigades/taz.shp', layer = 'taz')
# taz$Id <- 'taz'
# yam <- readOGR(dsn = 'Brigades/sovkhoz_yamalskii.shp', layer = 'sovkhoz_yamalskii')
# yam$Id <- 'yam'
# 
# # Iterate over polygons
# region_extracted <- lapply(c(taz,yam,yp), function(polygon){
# 
#   # Crop/mask to polygons and reproject to VSI
#   landforms <- mask(crop(geomorphons, polygon), polygon)
#   #veg_type <- mask(crop(landcover, polygon), polygon)
# 
#   beginCluster()
#   landforms <- projectRaster(landforms, crs=crs(vsi), method='ngb')
#   #veg_type <- projectRaster(veg_type, crs=crs(vsi), method='ngb')
#   endCluster()
#   gc()
# 
#   polygon_vsiproj <- spTransform(polygon, crs(vsi))
#   sensitivity <- mask(crop(vsi, polygon_vsiproj), polygon_vsiproj)
# 
#   # Turn VSI into polygons
#   sensitivity <- rasterToPolygons(sensitivity, na.rm=FALSE) #Consider in the future doing na.rm = FALSE
#   sensitivity <- st_as_sf(sensitivity)
# 
#   # Extract landforms and landcover
#   landform_VSI <- exact_extract(landforms, sensitivity)
#   #veg_type_VSI <- exact_extract(veg_type, sensitivity)
# 
#   # Write them out (later moved to 'Integrated' directory)
#   saveRDS(landform_VSI,paste0('landform_VSI_',polygon$Id,'.rds'))
#   #saveRDS(veg_type_VSI,paste0('veg_type_VSI_',polygon$Id,'.rds'))
# #  saveRDS(sensitivity$SensTotalW,paste0('VSI_as_polygons_',polygon$Id,'.rds'))
#   return(length(sensitivity$SensTotalW) == length(veg_type_VSI)) #sanity check
# })

################################################################################
# Figure 4A  
# VSI in uplands vs. lowlands vs. slopes vs. "rugged"
landforms_VSI <- lapply(c('taz','yp'), function(x){
  landform_list <- readRDS(paste0('Integrated/multivariateInRegions/landform_VSI_',x,'.rds'))
  reduce_landform <- lapply(landform_list, function(y){
    if (sum(y$coverage_fraction,na.rm=TRUE)==0){
      return(NA)
    }
    lowlandFrac <- sum(y$coverage_fraction[y$value>=8],na.rm=TRUE)/sum(y$coverage_fraction,na.rm=TRUE)
    uplandFrac <- sum(y$coverage_fraction[y$value<=4],na.rm=TRUE)/sum(y$coverage_fraction,na.rm=TRUE)
    slopeFrac <- sum(y$coverage_fraction[between(y$value,5,7)],na.rm=TRUE)/sum(y$coverage_fraction,na.rm=TRUE)
    thresh_single <- 2/3
    thresh_all <- 0.25
    if (lowlandFrac > thresh_single){
      'Lowlands'
    }
    else if (uplandFrac > thresh_single){
      'Uplands'
    }
    else if (slopeFrac > thresh_single){
      'Slopes'
    }
    else if (lowlandFrac > thresh_all & uplandFrac > thresh_all & slopeFrac > thresh_all){
      'Rugged'
    }
    else
      NA
  })
  reduce_landform_column <- do.call(rbind, reduce_landform)
  sensitivities <- readRDS(paste0('Integrated/multivariateInRegions/VSI_as_polygons_',x,'.rds'))
  sensitivity_by_landform <- data.frame(reduce_landform_column, sensitivities)
  names(sensitivity_by_landform) <- c('Landform','VSI')
  sensitivity_by_landform$Region <- x
  sensitivity_by_landform
})
landforms_VSI <- do.call(rbind, landforms_VSI)
landforms_VSI_full <- landforms_VSI
#landforms_VSI <- landforms_VSI[complete.cases(landforms_VSI),]

landforms_VSI$Region[landforms_VSI$Region == 'taz'] <- 'Taz'
landforms_VSI$Region[landforms_VSI$Region == 'yp'] <- 'LMA'
landforms_VSI <- filter(landforms_VSI, Region != 'yam')
landforms_VSI$Landform[is.na(landforms_VSI$Landform)] <- 'Other'

#Move 'rugged' into the 'other' category
#landforms_VSI$Landform[landforms_VSI$Landform == 'Rugged'] <- 'Other'
landforms_VSI <- filter(landforms_VSI, Landform %in% c('Uplands','Slopes', 'Lowlands'))

landforms_VSI$Landform <- factor(landforms_VSI$Landform, levels = c('Lowlands','Slopes', 'Uplands'))

give.n <- function(x){
  out <- c(y = quantile(x,0.75)-2, label = length(x))
  names(out) <- c('y','label')
  out
}

plot_a <- landforms_VSI %>% 
  ggplot(aes(x=Landform, y=VSI, fill=Region)) +
  geom_boxplot(outlier.shape=NA) +
  geom_signif(
    y_position = c(54,53,50), xmin = c(0.8,1.8,2.8), xmax = c(1.2,2.2,3.2),
    annotation = c("***","***","***"), tip_length = 0
  ) +
  xlab('Landform Class') +
  ylim(0,60) +
  ggtitle('Topography') +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge2(width = 0.75)) +
  scale_fill_discrete(breaks=c("LMA", "Taz"),
                      labels=c("Yar-Sale/Panayevsk  ", "Taz Peninsula")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),
        legend.text=element_text(size=14),axis.title = element_text(size = 14),
        plot.title=element_text(size=16)) +
  theme(legend.title=element_blank()) +
  theme(legend.position="top")

wilcox.test(VSI ~ Region, data = landforms_VSI[landforms_VSI$Landform == 'Upland',])

landform_aov <- anova_test(VSI~Region*Landform, data=landforms_VSI,type=2)
landform_aov <- aov(VSI~Region*Landform, data=landforms_VSI)
TukeyHSD(landform_aov,which="Landform") #All significant except lowland-slope (p=0.6088)
################################################################################
# Figure 4B
#VSI in different veg_types
sparse <- 1:2
tundra <- 3:5
shrub <- 6:10
meadow <- 11
floodplain <- 12:16
barren <- 17:18
veg_type_nums <- list(sparse, tundra, shrub, meadow, floodplain, barren)
veg_type_names <- c('Sparse','Prostrate Shrub','Tall Shrub','Meadow','Inundated','Rare Veg')
thresh_single <- 0.5
veg_type_VSI <- lapply(c('taz','yp'), function(x){
  veg_type_list <- readRDS(paste0('Integrated/multivariateInRegions/veg_type_VSI_',x,'.rds'))
  reduce_veg_type <- lapply(veg_type_list, function(y){
    if (sum(y$coverage_fraction,na.rm=TRUE)==0){
      return(NA)
    }
    veg_fractions <- unlist(lapply(1:length(veg_type_nums), function(veg){
      fraction <- sum(y$coverage_fraction[y$value %in% veg_type_nums[[veg]]],na.rm=TRUE)/sum(y$coverage_fraction,na.rm=TRUE)
    }))
    if (length(which(veg_fractions > thresh_single)) > 0){
      veg_type_names[which(veg_fractions > thresh_single)]
    }
    else {
      NA
    }
  })
  reduce_veg_type_column <- do.call(rbind, reduce_veg_type)
  sensitivities <- readRDS(paste0('Integrated/multivariateInRegions/VSI_as_polygons_',x,'.rds'))
  sensitivity_by_landform <- data.frame(reduce_veg_type_column, sensitivities)
  names(sensitivity_by_landform) <- c('Veg_Type','VSI')
  sensitivity_by_landform$Region <- x
  sensitivity_by_landform
})
veg_type_VSI <- do.call(rbind, veg_type_VSI)
veg_type_VSI_full <- veg_type_VSI
#veg_type_VSI <- veg_type_VSI[complete.cases(veg_type_VSI),]

# veg_type_VSI %>%
#   filter(Region == 'yam') %>%
#   group_by(Veg_Type) %>%
#   summarise(count = n())
veg_type_VSI$Region[veg_type_VSI$Region == 'taz'] <- 'Taz'
veg_type_VSI$Region[veg_type_VSI$Region == 'yp'] <- 'LMA'

#veg_type_VSI$Veg_Type[is.na(veg_type_VSI$Veg_Type)] <- 'Other'
veg_type_VSI$Veg_Type <- factor(veg_type_VSI$Veg_Type, levels = c('Inundated','Prostrate Shrub', 'Sparse', 'Tall Shrub'))
veg_type_VSI$Veg_Type <- recode(veg_type_VSI$Veg_Type,"Inundated"="Inund/Landslide",
                                       "Prostrate Shrub"="Low Stat Shrub", 'Tall Shrub' = 'Tall Woody')
veg_type_VSI <- veg_type_VSI[!is.na(veg_type_VSI$Veg_Type),]

plot_b <- veg_type_VSI %>% 
  #filter(Veg_Type != 'Meadow' & Veg_Type != 'Rare Veg') %>% 
  ggplot(aes(x=Veg_Type, y=VSI, fill=Region)) +
  geom_boxplot(outlier.shape=NA) +
  geom_signif(
    y_position = c(58,53,44,56), xmin = c(0.8,1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
    annotation = c("*","***","**","NS"), tip_length = 0
  ) +
  xlab('Landcover Class') +
  ggtitle('Landcover') +
  ylim(0,60) +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge2(width = 0.75)) +
  scale_fill_discrete(breaks=c("LMA", "Taz"),
                      labels=c("Yar-Sale/Panayevsk  ", "Taz Peninsula")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),
        legend.text=element_text(size=14),axis.title = element_text(size = 14),
        plot.title=element_text(size=16)) +
  theme(legend.position="top")


# veg_type_aov <- veg_type_VSI %>% 
#   filter(Region %in% c('taz','yp')) %>% 
#   filter(Veg_Type %in% c('shrub','sparse','tundra')) %>% 
#   drop_na() %>% 
#   anova_test(VSI~Region*Veg_Type, type=3, white.adjust=T )

wilcox.test(VSI ~ Region, data = veg_type_VSI[veg_type_VSI$Veg_Type == 'Tall Shrub',])

landcover_aov <- anova_test(VSI~Region*Veg_Type, data=veg_type_VSI, type=2)
landcover_aov <- aov(VSI~Region*Veg_Type, data=veg_type_VSI)
TukeyHSD(landcover_aov,which="Veg_Type") #All significant except tall shrub-inundated (p=0.0866)

################################################################################
# Figure 4C-E
# VSI in landforms and veg types
landforms_veg_type_VSI <- landforms_VSI_full
landforms_veg_type_VSI$Veg_Type <- veg_type_VSI_full$Veg_Type
landforms_veg_type_VSI <- drop_na(landforms_veg_type_VSI)
landforms_veg_type_VSI$Region[landforms_veg_type_VSI$Region == 'taz'] <- 'Taz'
landforms_veg_type_VSI$Region[landforms_veg_type_VSI$Region == 'yp'] <- 'LMA'

landforms_veg_type_VSI <- landforms_veg_type_VSI %>% 
  filter(Veg_Type %in% c('Inundated','Prostrate Shrub','Sparse','Tall Shrub')) %>%
  filter(Landform %in% c('Uplands','Lowlands','Slopes'))
landforms_veg_type_VSI$Veg_Type <- recode_factor(landforms_veg_type_VSI$Veg_Type,"Inundated"="I/L",
                                                 "Prostrate Shrub"="Low Stat Shr","Tall Shrub" = "Tall Woody")
landforms_veg_type_VSI$Veg_Type <- factor(landforms_veg_type_VSI$Veg_Type, 
                                          levels = c('I/L','Low Stat Shr', 'Sparse', 'Tall Woody'))

plot_landforms <- lapply(unique(landforms_veg_type_VSI$Landform), function(x){
  landforms_veg_type_VSI %>% 
    filter(Landform %in% x) %>%
    group_by(Veg_Type,Region) %>% 
    filter(n() > 2) %>% 
    ggplot(aes(x=Veg_Type, y=VSI, fill=Region)) +
    geom_boxplot(outlier.shape=NA) +
    xlab('Landcover Class') +
    ggtitle(x) +
    ylim(0,60) +
    stat_summary(fun.data = give.n, geom = "text", fun = median, 
                 position = position_dodge2(width = 0.75)) +
    scale_fill_discrete(breaks=c("LMA", "Taz"),
                        labels=c("Yar-Sale/Panayevsk  ", "Taz Peninsula")) +
    theme_bw()+
    theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),
          legend.text=element_text(size=14),axis.title = element_text(size = 14),
          plot.title=element_text(size=16)) +
    theme(legend.position="none")
})
#Use Wilcox test for pairwise comparisons
test_landforms <- lapply(unique(landforms_veg_type_VSI$Landform), function(x){
  curr_plot <- landforms_veg_type_VSI %>% 
    filter(Landform %in% x) %>%
    group_by(Veg_Type,Region) %>% 
    filter(n() > 2)
  curr_veg <- lapply(unique(curr_plot$Veg_Type), function(y){
    test_result <- wilcox.test(VSI ~ Region, data = curr_plot[curr_plot$Veg_Type == y,])
    c(x,as.character(y), test_result$p.value)
  })
})
#Add significances
plot_landforms[[1]] <- plot_landforms[[1]] +
  geom_signif(
    y_position = c(58,53,45,50), xmin = c(0.8,1.8,2.8,3.8), xmax = c(1.2,2.2,3.2,4.2),
    annotation = c("NS","***","*","NS"), tip_length = 0
  )
plot_landforms[[2]] <- plot_landforms[[2]] +
  geom_signif(
    y_position = c(52,39,49), xmin = c(0.8,1.8,2.8), xmax = c(1.2,2.2,3.2),
    annotation = c("***","NS","NS"), tip_length = 0
  )
plot_landforms[[3]] <- plot_landforms[[3]] +
  geom_signif(
    y_position = c(51,38), xmin = c(0.8,1.8), xmax = c(1.2,2.2),
    annotation = c("***","NS"), tip_length = 0
  )

plot_c <- ggarrange(plot_landforms[[1]],plot_landforms[[2]],plot_landforms[[3]],widths=c(6,5,4),
                    ncol=3,nrow=1,align="h",labels=c("C","D","E"),font.label=list(size=16))
plot_ab <- ggarrange(plot_a, plot_b, ncol=2,nrow=1,align="h",widths=c(5,6),common.legend = T,labels=c("A","B"),font.label=list(size=16))
plot_abc <- ggarrange(plot_ab,plot_c,ncol=1,nrow=2)

#ggsave('VSIcompareRegions_2022-11-18.pdf',width=12,height=8)


# data.frame(landforms_veg_type_filtered %>%
#              group_by(Region,Landform,Veg_Type) %>% 
#              get_summary_stats(VSI, type="common"))

# landforms_veg_type_filtered %>% 
#   group_by(Region,Landform,Veg_Type) %>% 
#   shapiro_test(VSI)

# VSI.aov <- landforms_veg_type_filtered %>% 
#   anova_test(VSI ~Region*Landform, type=3, white.adjust=T )
threeway_aov <- anova_test(VSI~Region*Veg_Type*Landform, data=landforms_veg_type_VSI)


################################################################################
#Incorporate reindeer for landcover comparison

yp <- readOGR(dsn = 'Brigades/yp.shp', layer = 'yp')
yp$Id <- 'yp'

vsi <- raster('MODIS/withTaz/VSI/Monthly_1km/SensTotalW.tif')

reindeer_june <- raster('Reindeer/Monthly/Merged/June_AllReindeerDensity.tif')
reindeer_july <- raster('Reindeer/Monthly/Merged/July_AllReindeerDensity.tif')
reindeer_august <- raster('Reindeer/Monthly/Merged/August_AllReindeerDensity.tif')
reindeer_jja <- raster('Reindeer/Monthly/Merged/JJA_AllReindeerDensity.tif')
reindeer <- list(reindeer_june, reindeer_july, reindeer_august, reindeer_jja)

reindeer_june_area <- raster('Reindeer/Monthly/Merged/June_AllReindeerDensity_LandArea.tif')
reindeer_july_area <- raster('Reindeer/Monthly/Merged/July_AllReindeerDensity_LandArea.tif')
reindeer_august_area <- raster('Reindeer/Monthly/Merged/August_AllReindeerDensity_LandArea.tif')
reindeer_jja_area <- raster('Reindeer/Monthly/Merged/JJA_AllReindeerDensity_LandArea.tif')
reindeer_area <- list(reindeer_june_area, reindeer_july_area, reindeer_august_area, reindeer_jja_area)

reindeer_yp <- lapply(reindeer, crop, spTransform(yp, crs(vsi)))
reindeer_area_yp <- lapply(reindeer_area, crop, spTransform(yp, crs(vsi)))
vsi <- crop(vsi, spTransform(yp,crs(vsi)))

reindeer_none <- is.na(reindeer_yp[[1]])&is.na(reindeer_yp[[2]])&is.na(reindeer_yp[[3]])
reindeer_none[reindeer_none == 0] <- NA

################################################################################
# Figure 5A
#Compare VSI within landcovers for different pasture areas (~treatments)

veg_type_list <- readRDS(paste0('Integrated/multivariateInRegions/veg_type_VSI_withNA_yp.rds'))

reduce_veg_type <- lapply(veg_type_list, function(y){
  if (sum(y$coverage_fraction,na.rm=TRUE)==0){
    return(NA)
  }
  veg_fractions <- unlist(lapply(1:length(veg_type_nums), function(veg){
    fraction <- sum(y$coverage_fraction[y$value %in% veg_type_nums[[veg]]],na.rm=TRUE)/sum(y$coverage_fraction,na.rm=TRUE)
  }))
  if (length(which(veg_fractions > thresh_single)) > 0){
    veg_type_names[which(veg_fractions > thresh_single)]
  }
  else {
    NA
  }
})
reduce_veg_type_column <- do.call(rbind, reduce_veg_type)

#Create a data frame for each treatment with VSI, veg, and reindeer density
veg_type_VSI_reindeer <- lapply(c(reindeer_yp,reindeer_none), function(x){ #Removed reindeer_may
  df <- data.frame(as.vector(vsi),reduce_veg_type_column,as.vector(x))
  names(df) <- c('VSI','Veg_Type','Reindeer')
  df
})

#Add the months and reindeer pasture areas (note the drop_na within the order)
months <- c('June','July','August','JJA')
for (i in 1:length(reindeer_area_yp)){
  veg_type_VSI_reindeer[[i]]$Area <- as.vector(reindeer_area_yp[[i]])
  veg_type_VSI_reindeer[[i]] <- drop_na(veg_type_VSI_reindeer[[i]])
  veg_type_VSI_reindeer[[i]]$Month <- months[i]
}
veg_type_VSI_reindeer[[5]] <- drop_na(veg_type_VSI_reindeer[[5]])
veg_type_VSI_reindeer[[5]]$Area <- NA
veg_type_VSI_reindeer[[5]]$Month <- 'None'

#Concatenate a single data frame
veg_type_VSI_reindeer <- do.call(rbind, veg_type_VSI_reindeer)

#Add the Taz landcover
add_Taz <- veg_type_VSI[veg_type_VSI$Region == 'Taz',]
names(add_Taz)[3] <- 'Month'
add_Taz$Reindeer <- NA
add_Taz$Area <- NA
add_Taz$Veg_Type <- recode_factor(add_Taz$Veg_Type,"Inund/Landslide"="Inundated",
                                  "Low Stat Shrub"="Prostrate Shrub","Tall Woody"="Tall Shrub")
veg_type_VSI_reindeer <- rbind(veg_type_VSI_reindeer, add_Taz)

veg_type_VSI_reindeer$Month <- factor(veg_type_VSI_reindeer$Month, 
                                      levels = c(months,'None','Taz'))
veg_type_VSI_reindeer$Veg_Type <- factor(veg_type_VSI_reindeer$Veg_Type,
                                         levels = c('Inundated','Prostrate Shrub','Sparse','Tall Shrub'))
veg_type_VSI_reindeer$Veg_Type <- recode_factor(veg_type_VSI_reindeer$Veg_Type,
                                                "Inundated"="Inundated/Landslide",
                                                "Prostrate Shrub"="Low Stature Shrub Tundra",
                                                "Sparse" = "Sparse Vegetation",
                                                "Tall Shrub" = "Tall Woody Vegetation")

plot_lc <- veg_type_VSI_reindeer %>% 
  filter(Veg_Type != 'Meadow' & Veg_Type != 'Rare Veg' & Veg_Type != 'Other') %>% 
  filter(Month != 'JJA') %>% 
  ggplot(aes(x=Veg_Type, y=VSI, fill=Month)) +
  geom_boxplot(outlier.shape=NA) +
  geom_signif(
    y_position = c(52,52,54,57.5,57.5,60,60), 
    xmin = c(2.0,2.15,1.7,1.85,1.85,1.7,1.7), 
    xmax = c(2.3,2.3,2.0,2.15,2.3,2.15,2.3),
    annotation = c("***","","**","","***","","***"), 
    tip_length = 0.01, vjust = 0.5
  )+
  geom_signif(
    y_position = c(60), 
    xmin = c(3.85), 
    xmax = c(4.15),
    annotation = c("**"), 
    tip_length = c(0.14,0.01), vjust = 0.5
  )+
  geom_signif(
    y_position = c(44,47,50,53), 
    xmin = c(3.15,3,2.85,2.7), 
    xmax = c(3.3,3.3,3.3,3.3),
    annotation = c("*","*","***","*"), 
    tip_length = c(0.03,0.01,0.05,0.03,0.08,0.05,0.14,0.07), vjust = 0.5
  )+
  xlab('Landcover Class') +
  ylim(0,60) +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge2(width = 0.75)) +
  scale_fill_discrete(breaks=c("June", "July","August","None","Taz"),
                      labels=c("June Pastures","July Pastures","August Pastures","No Pastures","Taz Peninsula")) +
  theme_bw()+
  theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),
        legend.text=element_text(size=13),axis.title = element_text(size = 14),
        plot.title=element_text(size=16)) +
  theme(legend.title=element_blank()) +
  theme(legend.position="top")

veg_type_VSI_reindeer_aov <- filter(veg_type_VSI_reindeer,Month != 'JJA')

reindeer_aov <- aov(VSI ~ Month, data = veg_type_VSI_reindeer_aov[veg_type_VSI_reindeer_aov$Veg_Type == 'Sparse Vegetation',])
summary(reindeer_aov)
TukeyHSD(reindeer_aov)

################################################################################
# Figure 5B,C
#Plot VSI against reindeer density in the tall shrub areas
plot_titles <- c('June Pastures: Tall Woody Vegetation','July Pastures','August Pastures','June-August Pastures: Tall Woody Vegetation')
plot_VSI_density <- lapply(1:length(months), function(x){
  veg_type_VSI_reindeer %>% 
    filter(Veg_Type == 'Tall Woody Vegetation') %>% 
    filter(Month == months[x]) %>% 
    filter(Area >= 10) %>% 
    ggplot(aes(x=Reindeer, y=VSI)) +
    geom_point() +
    geom_smooth(method = "lm") +
    #scale_size(range = c(1, 4),breaks = seq(100,500,100), limits=c(0,600), name=expression(paste("Area ", (km^2)))) +
    ylim(10,70) +
    {if (months[x] == 'JJA'){
      xlab(expression(paste((Reindeer%.%Months/km^2))))
    } else{
      xlab(expression(paste((Reindeer/km^2))))
    }} +
    ylab('VSI')+
    ggtitle(plot_titles[x])+
    theme_bw()+
    theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),
          legend.text=element_text(size=13),axis.title = element_text(size = 14),
          plot.title=element_text(size=14)) 
})

plot_points <- ggarrange(plot_VSI_density[[1]],plot_VSI_density[[4]], 
                         ncol=2, nrow=1, common.legend = T,legend='top',align='v',
                         labels=c('B','C'),font.label = list(size = 15))

ggarrange(plot_lc, plot_points, ncol=1, nrow=2, common.legend=F,
          labels=c('A',''),font.label = list(size = 15))
#ggsave('VSIcompareReindeer_2022-04-15.pdf',height=8,width=14)

reindeer_regression <- veg_type_VSI_reindeer %>% 
  filter(Veg_Type == 'Tall Woody Vegetation') %>% 
  filter(Month == 'JJA') %>% 
  filter(Area >= 10) %>% 
  lm(VSI ~ Reindeer, data = .)
summary(reindeer_regression)

################################################################################
# Figure S3
#Output the "pure pixels" as polygons
vsi <- raster('MODIS/withTaz/VSI/Monthly_1km/SensTotalW.tif')

# Vector polygons
yp <- readOGR(dsn = 'Brigades/yp.shp', layer = 'yp')
yp$Id <- 'yp'
taz <- readOGR(dsn = 'Brigades/taz.shp', layer = 'taz')
taz$Id <- 'taz'
yam <- readOGR(dsn = 'Brigades/sovkhoz_yamalskii.shp', layer = 'sovkhoz_yamalskii')
yam$Id <- 'yam'

# Iterate over polygons
vsi_polygons <- lapply(c(taz,yam,yp), function(polygon){
  
  polygon_vsiproj <- spTransform(polygon, crs(vsi))
  sensitivity <- mask(crop(vsi, polygon_vsiproj), polygon_vsiproj)
  
  # Turn VSI into polygons
  sensitivity <- rasterToPolygons(sensitivity)
})

taz_polygons <- cbind(vsi_polygons[[1]], landforms_VSI$Landform[landforms_VSI$Region=='taz'], 
                      veg_type_VSI$Veg_Type[landforms_VSI$Region=='taz'])
names(taz_polygons) <- c('VSI','Landform','Veg_Type')

yp_polygons <- cbind(vsi_polygons[[3]], landforms_VSI$Landform[landforms_VSI$Region=='yp'], 
                     veg_type_VSI$Veg_Type[landforms_VSI$Region=='yp'])
names(yp_polygons) <- c('VSI','Landform','Veg_Type')

yam_polygons <- cbind(vsi_polygons[[2]], landforms_VSI$Landform[landforms_VSI$Region=='yam'], 
                      veg_type_VSI$Veg_Type[landforms_VSI$Region=='yam'])
names(yam_polygons) <- c('VSI','Landform','Veg_Type')

shapefile(taz_polygons,'PurePixels_taz.shp')
shapefile(yam_polygons,'PurePixels_yam.shp')
shapefile(yp_polygons,'PurePixels_yp.shp')


################################################################################
# Figure S4B
#Try Drozdov
drozdov <- shapefile('Geology/RussianArcticSurfaceGeology.shp')
spTransform(drozdov, crs(vsi))

drozdov_regions <- lapply(c(taz,yp), function(polygon){
  polygon_drozdovproj <- spTransform(polygon, crs(drozdov))
  geology_drozdovproj <- crop(drozdov, polygon_drozdovproj)
  geology <- spTransform(geology_drozdovproj, crs(vsi))
  polygon_vsiproj <- spTransform(polygon, crs(vsi))
  sensitivity <- mask(crop(vsi,polygon_vsiproj),polygon_vsiproj)
  geology_sf <- st_as_sf(geology)
  drozdov_sensitivity <- exact_extract(sensitivity, geology_sf)
  list(geology, drozdov_sensitivity,polygon$Id)
})

drozdov_VSI <- lapply(drozdov_regions, function(polygon){
  lithology <- polygon[[1]]$Lithology_
  landscape <- polygon[[1]]$Landscape_
  VSI_values <- polygon[[2]]
  VSI_df <- lapply(1:length(VSI_values), function(x){
    VSI_values[[x]]$Lithology <- lithology[x]
    VSI_values[[x]]$Landscape <- landscape[x]
    VSI_values[[x]]$Region <- polygon[[3]]
    VSI_values[[x]]
  })
  do.call(rbind, VSI_df)
})
drozdov_VSI <- do.call(rbind, drozdov_VSI)
drozdov_VSI <- drop_na(drozdov_VSI)

drozdov_VSI$Region[drozdov_VSI$Region == 'taz'] <- 'Taz'
drozdov_VSI$Region[drozdov_VSI$Region == 'yp'] <- 'LMA'
drozdov_VSI %>% 
  filter(!Lithology %in% c("Peat and Sand")) %>% 
  mutate(Lithology = case_when(
    Lithology == 'Clay and Silt and Sand' ~ 'Clay/Silt\n Sand',
    Lithology == 'Clay and Silt' ~ 'Clay/Silt',
    Lithology == 'Peat and Sand' ~ 'PSa',
    Lithology == 'Peat, Sand, Clay and Silt' ~ 'Peat\n Sand\n Clay/Silt',
    Lithology == 'Peat, Clay, Silt, Sand' ~ 'Peat\n Clay/Silt\n Sand',
    Lithology == 'Sand' ~ 'Sand',
    Lithology == 'Sand, Clay and Silt' ~ 'Sand\n Clay/Silt'
  )) %>%
  ggplot(aes(x=Lithology,y=value,weight=coverage_fraction, fill=Region)) +
  geom_boxplot(outlier.shape=NA) +
  geom_signif(
    y_position = c(50,52,54,51,53), xmin = c(0.8,1.8,2.8,3.8,5.8), xmax = c(1.2,2.2,3.2,4.2,6.2),
    annotation = c("***","***","NS","NS","***"), tip_length = 0
  ) +
  scale_fill_discrete(breaks=c("LMA", "Taz"),
                      labels=c("Yar-Sale/Panayevsk", "Taz Peninsula")) +
  ylim(0,60) +
  ylab('VSI') +
  stat_summary(fun.data = give.n, geom = "text", fun = median, 
               position = position_dodge2(width = 0.75))+
  theme_bw()+
  theme(legend.title=element_blank()) +
  theme(legend.position="top")
ggsave('drozdovVSI_2022-03-12.pdf',width=10,height=6)

wilcox.test(value ~ Region, data = drozdov_VSI[drozdov_VSI$Lithology == 'Sand, Clay and Silt',])
drozdov_aov <- anova_test(value~Region*Lithology, 
                            data=drozdov_VSI[!drozdov_VSI$Lithology %in% c('Peat and Sand'),], type=2)

drozdov_aov <- aov(value~Region*Lithology, data=drozdov_VSI[!drozdov_VSI$Lithology %in% c("Peat and Sand"),])
TukeyHSD(drozdov_aov,which="Lithology") #All significant except lowland-slope (p=0.6088)


