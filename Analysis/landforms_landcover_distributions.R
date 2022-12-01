library(raster)
library(rgdal)
library(parallel)
library(tidyverse)
library(ggpubr)
library(sf)
library(exactextractr)
library(rgeos)
library(RColorBrewer)
library(viridis)

setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing")

landforms <- raster('ArcticDEM/withTaz/2020-06-11/Geomorphons_withTaz_32m_50cell_0.5deg_flatSplitThresholdExtendVKernelMajCRCRCR_UTM20.tif')
landcover <- raster('Pangaea/Landcover/ZAM_LCP_LANDC_SEN12_V02_20150815_20180830_T01.tif')

#----
# beginCluster()
# landforms <- projectRaster(geomorphons, landcover, method = 'ngb') #Align the rasters
# endCluster()
# gc()
# writeRaster(landforms, 'ArcticDEM/withTaz/2020-06-11/Geomorphons_withTaz_32m_50cell_0.5deg_flatSplitThresholdExtendVKernelMajCRCRCR_UTM20.tif',datatype='INT1U')
#----

################################################################################
# Figure 2A,C,E
# Landcover breakdown within landforms in YP vs Taz
sparse <- 1:2
tundra <- 3:5
shrub <- 6:10
meadow <- 11
floodplain <- 12:16
barren <- 17:18
cover_types <- list(sparse,tundra,shrub,meadow,floodplain,barren)
cover_names <- c('Sparse Vegetation','Low Stature Shrub Tundra','Tall Woody Vegetation','Meadow','Inundated/Landslide','Rare Veg')

uplands <- 0:4
slopes <- 5:7
lowlands <- 8:11
landscape_positions <- list(uplands, slopes, lowlands)
landscape_names <- list('upland','slope','lowland')


yp <- readOGR(dsn = 'Brigades/yp.shp', layer = 'yp')
yp$Id <- 'LMA'
taz <- readOGR(dsn = 'Brigades/taz.shp', layer = 'taz')
taz$Id <- 'Taz'

regions <- c(yp,taz)

veg_distribution_landforms <- lapply(regions, function(polygon){
  polygon <- st_as_sf(polygon)
  landforms_in_region <- exact_extract(landforms, polygon)[[1]]
  landcover_in_region <- exact_extract(landcover, polygon)[[1]]
  landforms_counts <- lapply(1:length(landscape_positions), function(x){
    landcover_counts <- lapply(1:length(cover_types), function(y){
      count_pixels <- sum(landforms_in_region$coverage_fraction[landforms_in_region$value %in% landscape_positions[[x]] & 
                                                                  landcover_in_region$value%in% cover_types[[y]]])
      out <- data.frame(count_pixels, landscape_names[x], cover_names[y], polygon$Id)
    })
    out <- do.call(rbind, landcover_counts)
    names(out) <- c('Pixels','Landforms','Landcover','Region')
    out
  })
  landforms_landcover_counts <- do.call(rbind, landforms_counts)
})
landforms_landcover_regions <- do.call(rbind, veg_distribution_landforms)

landforms_landcover_regions$Region <- factor(landforms_landcover_regions$Region)
plot_regions <- lapply(unique(landforms_landcover_regions$Landforms), function(lf){
  landforms_landcover_regions %>%
    filter(Landforms == lf) %>% 
    filter(!Landcover %in% c('Rare Veg', 'Meadow')) %>% #Rare veg < 8%, Meadow < 1%
    group_by(Region) %>% 
    mutate(Percentage = Pixels/sum(Pixels)*100) %>% 
    arrange(Landforms, desc(Landcover)) %>% 
    mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage) %>%
    mutate(n_region = sum(Pixels)) %>% 
    ggplot(aes(x = Region, y = Percentage)) +
    geom_col(aes(fill = Landcover),width=0.4) +
    scale_fill_manual(values=c("#999999", "#228B22", "#90EE90","#CF9FFF")) +
    ylab('Percentage of Pixels') +
    scale_x_discrete(labels=c("LMA" = "Yar-Sale/Panayevsk", "Taz" = "Taz Peninsula")) +
    xlab('Region') +
    geom_text(aes(y = lab_ypos, label = paste0(sprintf('%0.1f',Percentage)), group = Landcover),size=3.5) +
    geom_text(aes(y = 110, label = paste0('n=',sprintf("%0.1e",n_region)), group = Landcover), size = 3) +
    theme_bw()+
    theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),legend.text=element_text(size=10)) 
  
})

plotTitles <- c('Uplands','Slopes','Lowlands')
for (i in 1:length(plot_regions)){
  plot_regions[[i]] <- plot_regions[[i]] + ggtitle(plotTitles[i])
}

pr <- ggarrange(plot_regions[[1]],plot_regions[[2]],plot_regions[[3]], 
          ncol=1, nrow=3, common.legend = T, align="h")


################################################################################
# Do the statistics for region vs landforms vs landcover
landforms_landcover_regions_filter <- filter(landforms_landcover_regions, Landcover %in% c('Tall Shrub', 'Sparse'))
cont_table <- xtabs(Pixels ~ Landcover + Region  + Landforms, data = landforms_landcover_regions_filter)

### Test the Mutual Independence, step by step 
### compute the expected frequences
E=array(NA,dim(cont_table))
for (i in 1:dim(cont_table)[1]) {
  for (j in 1:dim(cont_table)[2]) {
    for (k in 1:dim(cont_table)[3]) {
      E[i,j,k]=(margin.table(cont_table,3)[k]*margin.table(cont_table,2)[j]*margin.table(cont_table,1))[i]/(sum(cont_table))^2
    }}}
E

### compute the X^2, and G^2
df=(length(cont_table)-1)-(sum(dim(cont_table))-3)
X2=sum((cont_table-E)^2/E)
X2
1-pchisq(X2,df)
G2=2*sum(cont_table*log(cont_table/E))
G2
1-pchisq(G2,df)

### Test for Mutual indpendence (and other models) by considering anlaysis of all two-way tables
landcover_region=margin.table(cont_table,c(1,2))
result=chisq.test(landcover_region, correct=F)
result
result$expected

region_landform=margin.table(cont_table,c(2,3))
result=chisq.test(region_landform, correct=F)
result
result$expected

landcover_landform=margin.table(cont_table,c(1,3))
result=chisq.test(landcover_landform, correct=F)
result
result$expected

library(vcd)
## compute the log(oddsraio), oddsratio and its 95% CI using {vcd} package
or_region_landform = oddsratio(region_landform, log=FALSE)
ci_region_landform = confint(or_region_landform)

or_landcover_region = oddsratio(landcover_region, log=FALSE)
ci_landcover_region = confint(or_landcover_region)

or_landcover_landform = oddsratio(landcover_landform, log=FALSE)
ci_landcover_landform = confint(or_landcover_landform)

##################
#Log-linear models
##################
landforms_landcover_regions_filter$Pixels <- round(landforms_landcover_regions_filter$Pixels)

#Saturated Model
#modelSat=glm(Pixels~Landcover+Landforms+Region+Landcover*Landforms+Landcover*Region+Landforms*Region
#            +Landcover*Landforms*Region,family=poisson(link=log),data = landforms_landcover_regions_filter)
modelSat = glm(Pixels~Landcover*Landforms*Region,family=poisson(link=log),data = landforms_landcover_regions_filter)
summary(modelSat)

#Null Model
modelNull=glm(Pixels~Landcover+Landforms+Region+Landforms*Region
              ,family=poisson(link=log),data = landforms_landcover_regions_filter)
summary(modelNull)

#Conditional, Landforms have no effect on landcover after region is included
modelCon1=glm(Pixels~Landcover+Landforms+Region+Landcover*Region+Landforms*Region
              ,family=poisson(link=log),data = landforms_landcover_regions_filter)
summary(modelCon1)

#Conditional, Region has no effect on landcover after landforms are included
modelCon2=glm(Pixels~Landcover+Landforms+Region+Landcover*Landforms+Landforms*Region
              ,family=poisson(link=log),data = landforms_landcover_regions_filter)
summary(modelCon2)

#Homogenous Association
modelHom=glm(Pixels~Landcover+Landforms+Region+Landcover*Landforms+Landcover*Region+Landforms*Region
             ,family=poisson(link=log),data = landforms_landcover_regions_filter)
summary(modelHom)

library(MuMIn)
Weights(c(AIC(modelSat),AIC(modelNull),AIC(modelCon1),AIC(modelCon2),AIC(modelHom)))

#Compare AIC and BIC
models <- list(modelNull, modelCon1, modelCon2, modelHom)
n_params <- c(5,6,6,7)
aic_vals <- sapply(1:length(models), function(x){
  G2 <- 2*(logLik(modelSat)-logLik(models[[x]]))
  G2[1]-2*n_params[x]
})
bic_vals <- sapply(models, function(x){
  AIC(x,k=log(sum(landforms_landcover_regions_filter$Pixels)))-
    AIC(modelSat,k=log(sum(landforms_landcover_regions_filter$Pixels)))
})
model_names <- c('Null','No_Landforms','No_Region','HomogenousAssoc')
data.frame(model_names, aic_vals, bic_vals)


################################################################################
#Figure 2B,D,F  
#Repeat but with reindeer polygons

reindeerJune <- shapefile('Reindeer/Monthly/Merged/June_AllReindeer.shp')
reindeerJuly <- shapefile('Reindeer/Monthly/Merged/July_AllReindeer.shp')
reindeerAugust <- shapefile('Reindeer/Monthly/Merged/August_AllReindeer.shp')
reindeerNone <- bind(reindeerJune, reindeerJuly, reindeerAugust)

reindeerMonth <- list(reindeerJune, reindeerJuly, reindeerAugust,reindeerNone)
reindeerMonth <- lapply(reindeerMonth, spTransform, crs(landforms))
reindeerMonthUnion <- lapply(reindeerMonth, gUnaryUnion)
reindeerMonthUnion[[5]] <- gDifference(yp,reindeerMonthUnion[[4]])

month_names <- c('June','July','August','All','None')

veg_distribution_landforms_month <- lapply(1:length(reindeerMonthUnion), function(i){
  polygon <- st_as_sf(reindeerMonthUnion[[i]])
  landforms_in_region <- exact_extract(landforms, polygon)[[1]]#There should be only one polygon so [[1]] unlists
  landcover_in_region <- exact_extract(landcover, polygon)[[1]]
  landforms_counts <- lapply(1:length(landscape_positions), function(x){
    landcover_counts <- lapply(1:length(cover_types), function(y){
      count_pixels <- sum(landforms_in_region$coverage_fraction[landforms_in_region$value %in% landscape_positions[[x]] & 
                                                                  landcover_in_region$value%in% cover_types[[y]]])
      out <- data.frame(count_pixels, landscape_names[x], cover_names[y], month_names[i])
    })
    out <- do.call(rbind, landcover_counts)
    names(out) <- c('Pixels','Landforms','Landcover','Month')
    out
  })
  landforms_landcover_counts <- do.call(rbind, landforms_counts)
})
landforms_landcover_month <- do.call(rbind, veg_distribution_landforms_month)


landforms_landcover_month$Month <- factor(landforms_landcover_month$Month, levels=c('June','July','August','None','All'))
plot_months <- lapply(unique(landforms_landcover_month$Landforms), function(lf){
  landforms_landcover_month %>%
    filter(Month != 'All') %>%
    filter(Landforms == lf) %>% 
    filter(!Landcover %in% c('Rare Veg', 'Meadow')) %>% #Rare veg < 8%, Meadow < 1%
    group_by(Month) %>% 
    mutate(Percentage = Pixels/sum(Pixels)*100) %>% 
    arrange(Landforms, desc(Landcover)) %>% 
    mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage) %>%
    mutate(n_month = sum(Pixels)) %>% 
    ggplot(aes(x = Month, y = Percentage)) +
    geom_col(aes(fill = Landcover),width=0.5) +
    scale_fill_manual(values=c("#999999", "#228B22", "#90EE90","#CF9FFF")) +
    ylab('Percentage of Pixels') +
    scale_x_discrete(labels=c("June" = "June Pastures", "July" = "July Pastures",
                                "August" = "August Pastures", "None" = "No Pastures")) +
    xlab('Yar-Sale/Panayevsk Pasture Areas') +
    geom_text(aes(y = lab_ypos, label = paste0(sprintf('%0.1f',Percentage)), group = Landcover),size=3.5) +
    geom_text(aes(y = 110, label = paste0('n=',sprintf("%0.1e",n_month)), group = Landcover), size = 3) +
    theme_bw()+
    theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),legend.text=element_text(size=10)) 
    
})

plotTitles <- c('Uplands','Slopes','Lowlands')
for (i in 1:length(plot_months)){
  plot_months[[i]] <- plot_months[[i]] + ggtitle(plotTitles[i])
}

pm <- ggarrange(plot_months[[1]],plot_months[[2]],plot_months[[3]], 
          ncol=1, nrow=3, common.legend = T, align="h")

ggarrange(plot_regions[[1]],plot_months[[1]],plot_regions[[2]],plot_months[[2]],plot_regions[[3]],plot_months[[3]],
          ncol=2,nrow=3,common.legend=T,align='hv',labels="AUTO",widths = c(2, 3))

#ggsave('landcoverDistributions_2022-03-13.pdf',width=10,height=7)

################################################################################
# Do the statistics for month vs landforms vs landcover
landforms_landcover_month_filter <- landforms_landcover_month %>% 
  filter(Landcover %in% c('Tall Shrub', 'Sparse')) %>% 
  filter(Month != 'All')

landforms_landcover_month_filter$Pixels <- round(landforms_landcover_month_filter$Pixels)

#Saturated Model
#modelSat=glm(Pixels~Landcover+Landforms+Region+Landcover*Landforms+Landcover*Region+Landforms*Region
#            +Landcover*Landforms*Region,family=poisson(link=log),data = landforms_landcover_regions_filter)
modelSat = glm(Pixels~Landcover*Landforms*Month,family=poisson(link=log),data = landforms_landcover_month_filter)
summary(modelSat)

#Null Model
modelNull=glm(Pixels~Landcover+Landforms+Month+Landforms*Month
              ,family=poisson(link=log),data = landforms_landcover_month_filter)
summary(modelNull)

#Conditional, Landforms have no effect on landcover after Month is included
modelCon1=glm(Pixels~Landcover+Landforms+Month+Landcover*Month+Landforms*Month
              ,family=poisson(link=log),data = landforms_landcover_month_filter)
summary(modelCon1)

#Conditional, Month has no effect on landcover after landforms are included
modelCon2=glm(Pixels~Landcover+Landforms+Month+Landcover*Landforms+Landforms*Month
              ,family=poisson(link=log),data = landforms_landcover_month_filter)
summary(modelCon2)

#Homogenous Association
modelHom=glm(Pixels~Landcover+Landforms+Month+Landcover*Landforms+Landcover*Month+Landforms*Month
             ,family=poisson(link=log),data = landforms_landcover_month_filter)
summary(modelHom)


################################################################################
#Figure 3A,B
#Now incorporate the reindeer densities

#1. Calculate anticipated ratios based on landforms
sparse_frac_june <- landforms_landcover_month %>%
  filter(Month == 'June') %>% 
  group_by(Landforms) %>%
  summarise(ratio = Pixels[Landcover=='Sparse Vegetation']/(Pixels[Landcover=='Tall Woody Vegetation']+Pixels[Landcover=='Sparse Vegetation'])) %>% 
  pivot_wider(ratio,names_from=Landforms,values_from=ratio)

#2. Change reindeer numbers into densities
orig_poly_veg <- exact_extract(landcover, st_as_sf(reindeerMonth[[1]])) #Find land pixels
land_area <- sapply(orig_poly_veg, function(x){
  sum(x$coverage_fraction[x$value %in% unlist(cover_types)])*res(landcover)[1]*res(landcover)[2]/10^6
})
reindeer_density <- reindeerMonth[[1]]$Reindeer/land_area

#3. Intersect overlapping polygons
reindeerJuneUnion <- raster::union(reindeerMonth[[1]])
reindeerJune_sf <- st_as_sf(reindeerJuneUnion)
june_veg <- exact_extract(landcover, reindeerJune_sf)
june_lf <- exact_extract(landforms, reindeerJune_sf)

june_veg_density <- lapply(1:length(june_veg),function(x){
  june_extract <- data.frame(june_lf[[x]]$value, june_veg[[x]]$value, june_veg[[x]]$coverage_fraction)
  names(june_extract) <- c('landforms','veg','coverage')
  
  june_extract <- june_extract[complete.cases(june_extract),] #Filter out NA
  june_extract <- june_extract[june_extract$veg %in% unlist(cover_types),] #Filter out water
  n <- sum(june_extract$coverage[june_extract$veg %in% unlist(cover_types)])
  
  #4. Get the proportion of landforms
  upland_count <- sum(june_extract$coverage[june_extract$landforms %in% uplands])
  slope_count <- sum(june_extract$coverage[june_extract$landforms %in% slopes])
  lowland_count <- sum(june_extract$coverage[june_extract$landforms %in% lowlands])
  
  #5. Get the numbers of sparse and tall shrub pixels
  sparse_count <- sum(june_extract$coverage[june_extract$veg %in% sparse])
  shrub_count <- sum(june_extract$coverage[june_extract$veg %in% shrub])
  
  upland_sparse_count <- sum(june_extract$coverage[(june_extract$veg %in% sparse)
                                                   &(june_extract$landforms %in% uplands)])
  upland_shrub_count <- sum(june_extract$coverage[(june_extract$veg %in% shrub)
                                                  &(june_extract$landforms %in% uplands)])
  
  
  #6. Get the reindeer numbers and density for the polygon
  #The union has one row for each of the new polygons and one column for each of the old polygons.
  #A "1" entry in a column indicates that old polygon comprises the new polygon for that row
  polys_included <- which(data.frame(reindeerJuneUnion[x,1:(ncol(reindeerJuneUnion)-1)]) == 1)
  density <- sum(reindeer_density[polys_included])
  reindeer <- sum(reindeerMonth[[1]]$Reindeer[polys_included])
  
  c(n,lowland_count,slope_count,upland_count,sparse_count, shrub_count, 
    upland_sparse_count, upland_shrub_count, density, reindeer)
  
})
june_veg_density <- data.frame(do.call(rbind, june_veg_density))
names(june_veg_density) <- c('n','lowland','slope','upland','sparse','shrub',
                             'upland_sparse', 'upland_shrub','reindeer_density','reindeer')

#7. Find the expected sparse frac based on the landform dist and subtract from the actual sparse frac
june_reindeer_effect <- june_veg_density %>% 
  filter(n > 25000) %>% #If resolution is 20x20, then this implies > 10km^2
  mutate(expected_sparse = (lowland*sparse_frac_june$lowland+slope*sparse_frac_june$slope+upland*sparse_frac_june$upland)/
           (lowland+slope+upland)) %>% 
  mutate(sparse_diff = (sparse/(sparse+shrub))-expected_sparse)

#8. Plot the difference in the sparse frac against the reindeer density with points sized by area
getPalette=colorRampPalette(brewer.pal(9,"Greens"))

june_reindeer_effect$area <- june_reindeer_effect$n*res(landcover)[1]*res(landcover)[2]/10^6
p1 <- ggplot(june_reindeer_effect,aes(x=reindeer_density, y=sparse_diff, size=area)) +
  geom_hline(yintercept=0,color='darkgray') +
  geom_point() +
  scale_size(range = c(1, 4), name=expression(paste("Area ", (km^2)))) +
  geom_density_2d(show.legend=F,color='black') +
  geom_density_2d_filled(alpha = 0.4,show.legend=F) +
  scale_fill_manual(values=getPalette(11)) +
  xlab(expression(paste((Reindeer/km^2)))) +
  ylab('')+
  ylim(-0.8,0.4) +
  ggtitle('June: All Landforms')+
  theme_bw()+
  theme(text = element_text(size=11))

#7. Find the expected sparse frac based on the landform dist and subtract from the actual sparse frac
june_reindeer_effect_uplands <- june_veg_density %>% 
  filter(n > 25000) %>% #If resolution is 20x20, then this implies > 10km^2
  mutate(sparse_diff = (upland_sparse/(upland_sparse+upland_shrub))-sparse_frac_june$upland)

#8. Plot the difference in the sparse frac against the reindeer density with points sized by area
june_reindeer_effect_uplands$area <- june_reindeer_effect_uplands$n*res(landcover)[1]*res(landcover)[2]/10^6
p2 <- ggplot(june_reindeer_effect_uplands,aes(x=reindeer_density, y=sparse_diff, size=area)) +
  geom_hline(yintercept=0,color='darkgray') +
  geom_point() +
  scale_size(range = c(1, 4), name=expression(paste("Area ", (km^2)))) +
  geom_density_2d(show.legend=F,color='black') +
  geom_density_2d_filled(alpha = 0.4,show.legend=F) +
  scale_fill_manual(values=getPalette(14)) +
  xlab(expression(paste((Reindeer/km^2)))) +
  ylab('')+
  ylim(-0.8,0.4) +
  ggtitle('June: Uplands')+
  theme_bw()+
  theme(text = element_text(size=11))


################################################################################
#Figure 3C,D
#Now incorporate the reindeer densities
#And repeat but instead of June, use all the reindeer

#1. Calculate anticipated ratios based on landforms
sparse_frac <- landforms_landcover_month %>%
  filter(Month == 'All') %>% 
  group_by(Landforms) %>%
  summarise(ratio = Pixels[Landcover=='Sparse Vegetation']/(Pixels[Landcover=='Tall Woody Vegetation']+Pixels[Landcover=='Sparse Vegetation'])) %>% 
  pivot_wider(ratio,names_from=Landforms,values_from=ratio)

#2. Calculate reindeer numbers and densities in intersected polygons
#Doing something funky because union not working for reindeerAll
unions <- lapply(1:3, function(x){
  
  orig_poly_veg <- exact_extract(landcover, st_as_sf(reindeerMonth[[x]])) #Find land pixels
  land_area <- sapply(orig_poly_veg, function(y){
    sum(y$coverage_fraction[y$value %in% unlist(cover_types)])*res(landcover)[1]*res(landcover)[2]/10^6
  })
  month_density <- reindeerMonth[[x]]$Reindeer/land_area
  
  month_union <- raster::union(reindeerMonth[[x]])
  month_union$reindeer_density <- 0
  month_union$reindeer_total <- 0
  for (y in 1:nrow(month_union)){
    polys_included <- which(data.frame(month_union[y,1:nrow(reindeerMonth[[x]])]) == 1)
    month_union$reindeer_density[y] <- sum(month_density[polys_included])
    month_union$reindeer_total[y] <- sum(reindeerMonth[[x]]$Reindeer[polys_included])
  }
  month_union[,(nrow(reindeerMonth[[x]])+1):ncol(month_union)]
})
unionJJ <- raster::union(unions[[1]], unions[[2]])
unionJJA <- raster::union(unionJJ, unions[[3]])
unionJJA$count <- rowSums(data.frame(unionJJA[,c(1,4,7)]),na.rm=TRUE)
unionJJA$reindeer_density <- rowSums(data.frame(unionJJA[,c(2,5,8)]),na.rm=TRUE)
unionJJA$reindeer_total <- rowSums(data.frame(unionJJA[,c(3,6,9)]),na.rm=TRUE)
unionJJA <- unionJJA[,7:9]

reindeer_sf <- st_as_sf(unionJJA)
reindeer_veg <- exact_extract(landcover, reindeer_sf)
reindeer_lf <- exact_extract(landforms, reindeer_sf)

reindeer_veg_density <- lapply(1:length(reindeer_veg),function(x){
  reindeer_extract <- data.frame(reindeer_lf[[x]]$value, reindeer_veg[[x]]$value, reindeer_veg[[x]]$coverage_fraction)
  names(reindeer_extract) <- c('landforms','veg','coverage')
  
  reindeer_extract <- reindeer_extract[complete.cases(reindeer_extract),] #Filter out NA
  reindeer_extract <- reindeer_extract[reindeer_extract$veg %in% unlist(cover_types),] #Filter out water
  n <- sum(reindeer_extract$coverage[reindeer_extract$veg %in% unlist(cover_types)])
  
  #4. Get the proportion of landforms
  upland_count <- sum(reindeer_extract$coverage[reindeer_extract$landforms %in% uplands])
  slope_count <- sum(reindeer_extract$coverage[reindeer_extract$landforms %in% slopes])
  lowland_count <- sum(reindeer_extract$coverage[reindeer_extract$landforms %in% lowlands])
  
  #5. Get the numbers of sparse and tall shrub pixels
  sparse_count <- sum(reindeer_extract$coverage[reindeer_extract$veg %in% sparse])
  shrub_count <- sum(reindeer_extract$coverage[reindeer_extract$veg %in% shrub])
  
  upland_sparse_count <- sum(reindeer_extract$coverage[(reindeer_extract$veg %in% sparse)
                                                       &(reindeer_extract$landforms %in% uplands)])
  upland_shrub_count <- sum(reindeer_extract$coverage[(reindeer_extract$veg %in% shrub)
                                                      &(reindeer_extract$landforms %in% uplands)])
  
  #6. Get the reindeer numbers and density for the polygon
  density <- unionJJA$reindeer_density[[x]]
  reindeer <- unionJJA$reindeer_total[[x]]
  
  c(n,lowland_count,slope_count,upland_count,sparse_count, shrub_count, 
    upland_sparse_count, upland_shrub_count, density, reindeer)
  
})
reindeer_veg_density <- data.frame(do.call(rbind, reindeer_veg_density))
names(reindeer_veg_density) <- c('n','lowland','slope','upland','sparse','shrub',
                                 'upland_sparse', 'upland_shrub','reindeer_density','reindeer')

#7. Find the expected sparse frac based on the landform dist and subtract from the actual sparse frac
reindeer_effect <- reindeer_veg_density %>% 
  filter(n > 25000) %>% #If resolution is 20x20, then this implies > 10km^2
  mutate(expected_sparse = (lowland*sparse_frac$lowland+slope*sparse_frac$slope+upland*sparse_frac$upland)/
           (lowland+slope+upland)) %>% 
  mutate(sparse_diff = (sparse/(sparse+shrub))-expected_sparse)

#8. Plot the difference in the sparse frac against the reindeer density with points sized by area
reindeer_effect$area <- reindeer_effect$n*res(landcover)[1]*res(landcover)[2]/10^6
p3 <- ggplot(reindeer_effect,aes(x=reindeer_density, y=sparse_diff, size=area)) +
  geom_hline(yintercept=0,color='darkgray') +
  geom_point() +
  scale_size(range = c(1, 4), name=expression(paste("Area ", (km^2)))) +
  geom_density_2d(show.legend=F,color='black') +
  geom_density_2d_filled(alpha = 0.4,show.legend=F) +
  scale_fill_manual(values=getPalette(11)) +
  xlab(expression(paste((Reindeer%.%Months/km^2)))) +
  ylab('')+
  ylim(-0.8,0.4) +
  ggtitle('June-August: All Landforms')+
  theme_bw()+
  theme(text = element_text(size=11))

#7. Find the expected sparse frac based on the landform dist and subtract from the actual sparse frac
reindeer_effect_uplands <- reindeer_veg_density %>% 
  filter(n > 25000) %>% #If resolution is 20x20, then this implies > 10km^2
  mutate(sparse_diff = (upland_sparse/(upland_sparse+upland_shrub))-sparse_frac$upland)

#8. Plot the difference in the sparse frac against the reindeer density with points sized by area
reindeer_effect_uplands$area <- reindeer_effect_uplands$n*res(landcover)[1]*res(landcover)[2]/10^6
p4 <- ggplot(reindeer_effect_uplands,aes(x=reindeer_density, y=sparse_diff, size=area)) +
  geom_hline(yintercept=0,color='darkgray') +
  geom_point() +
  scale_size(range = c(1, 4), name=expression(paste("Area ", (km^2)))) +
  geom_density_2d(show.legend=F,color='black') +
  geom_density_2d_filled(alpha = 0.4,show.legend=F) +
  scale_fill_manual(values=getPalette(8)) +
  xlab(expression(paste((Reindeer%.%Months/km^2)))) +
  ylab('')+
  ylim(-0.8,0.4) +
  ggtitle('June-August: Uplands')+
  theme_bw()+
  theme(text = element_text(size=11))

p <- ggarrange(p1,p2,p3,p4, 
               ncol=2, nrow=2, common.legend = T, labels="AUTO",hjust=-2,legend="top")
annotate_figure(p,
                bottom = text_grob("Animal Density", size = 12,vjust=-0.5),
                left = text_grob("Difference from Expectation in Sparse:Tall Woody Fraction", size = 12,rot = 90,vjust=2.5)
)
#ggsave('landcoverVsReindeerDensity_2022-11-18.pdf',height = 7, width = 7)
