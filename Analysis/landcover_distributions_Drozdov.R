library(raster)
library(rgdal)
library(tidyverse)
library(ggpubr)
library(parallel)
library(sf)
library(exactextractr)

setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing")

drozdov <- shapefile('Geology/RussianArcticSurfaceGeology.shp')
landcover <- raster('Pangaea/Landcover/ZAM_LCP_LANDC_SEN12_V02_20150815_20180830_T01.tif')

yp <- readOGR(dsn = 'Brigades/yp.shp', layer = 'yp')
yp$Id <- 'Yar-Sale/Panayevsk'
taz <- readOGR(dsn = 'Brigades/taz.shp', layer = 'taz')
taz$Id <- 'Taz Peninsula'

sparse <- 1:2
tundra <- 3:5
shrub <- 6:10
meadow <- 11
floodplain <- 12:16
barren <- 17:18
cover_types <- list(sparse,tundra,shrub,meadow,floodplain,barren)
cover_names <- c('Sparse Vegetation','Low Stature Shrub Tundra','Tall Woody Vegetation','Meadow','Inundated/Landslide','Rare Veg')

lithologies <- c('Clay and Silt and Sand', 'Clay and Silt','Peat and Sand','Peat, Sand, Clay and Silt',
                 'Peat, Clay, Silt, Sand','Sand','Sand, Clay and Silt')
lithologies_short <- c('CS','C','PS','PSC','PCS','S','SC')

landcover_by_lithology_polygon <- lapply(c(yp,taz), function(polygon){
  poly <- spTransform(polygon, crs(drozdov))
  drozdov_poly <- crop(drozdov, poly)
  drozdov_proj <- spTransform(drozdov_poly, crs(landcover))
  
  indices <- lapply(lithologies, function(x){
    which(drozdov_proj$Lithology_ == x)
  })
  
  landcover_by_lithology <- lapply(1:length(indices), function(x){
    landcover_extract <- exact_extract(landcover,st_as_sf(drozdov_proj[indices[[x]],]))
    landcover_extract <- do.call(rbind, landcover_extract)
    
    landcover_counts <- lapply(1:length(cover_types), function(y){
      count_pixels <- sum(landcover_extract$coverage_fraction[landcover_extract$value %in% cover_types[[y]]])
      data.frame(count_pixels, cover_names[y], polygon$Id)
    })
    out <- do.call(rbind, landcover_counts)
    names(out) <- c('Pixels','Landcover','Region')
    out$Lithology <- lithologies_short[[x]]
    out
  })
  do.call(rbind,landcover_by_lithology)
})
landcover_totals <- do.call(rbind, landcover_by_lithology_polygon)

landcover_totals <- filter(landcover_totals, Lithology!='PS') 

plot_regions <- lapply(c('Yar-Sale/Panayevsk','Taz Peninsula'), function(region){
  landcover_totals %>%
    filter(Region == region) %>% 
    filter(!Landcover %in% c('Rare Veg', 'Meadow')) %>%
    mutate(Percent_All = Pixels/sum(Pixels)) %>% 
    group_by(Lithology) %>% 
    mutate(Percentage = Pixels/sum(Pixels)) %>%
    mutate(n_region = sum(Pixels)) %>% 
    drop_na() %>% 
    #group_by(Position) %>% 
    arrange(Lithology, desc(Landcover)) %>% 
    mutate(lab_ypos = cumsum(Percent_All) - 0.5 * Percent_All) %>%
    mutate(lab_ypos = #To deal with overlapping numbers
             case_when(
               lab_ypos < 0.003 ~ -0.003,
               TRUE ~ lab_ypos
             )) %>% 
    mutate(lab_ypos2 = sum(Percent_All) +0.01) %>% 
    ggplot(aes(x = Lithology, y = Percent_All)) +
    geom_col(aes(fill = Landcover)) +
    scale_fill_manual(values=c("#999999", "#228B22", "#90EE90","#CF9FFF")) +
    scale_x_discrete(labels=c("C" = "Clay/Silt", "CS" = "Clay/Silt\n Sand", "PCS" = "Peat\n Clay/Silt\n Sand",
                              "PSC" = "Peat\n Sand\n Clay/Silt", "S" = "Sand", "SC" = "Sand\n Clay/Silt")) +
    geom_text(aes(y = lab_ypos, label = sprintf('%0.1f',Percentage*100), group = Landcover)) +
    geom_text(aes(y = lab_ypos2, label = paste0('n=',sprintf("%0.1e",n_region)), group = Landcover), size = 3) +
    ylim(-0.01,0.45) +
    ylab('Percentage of Pixels') +
    theme_bw()+
    ggtitle(region)
})
drozdov_landcover <- ggarrange(plot_regions[[1]], plot_regions[[2]], ncol=2, nrow=1, common.legend = T, align="h", widths=c(6,5))
#ggsave('drozdovLandcoverDistributions.pdf',width=10,height=7)
