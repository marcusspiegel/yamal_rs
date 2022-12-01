library(raster)
library(rgdal)
library(spatialEco)

# setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing")
setwd("E:/Data/Biogeosciences/Yamal/")

################################################################################
#Geomorphons distributions
vsi <- raster('MODIS/withTaz/VSI/rasters_Monthly_NoT1/SensTotalW.tif')
geomorphons <- raster('ArcticDEM/landforms/withTaz/2020-06-11/Geomorphons_withTaz_32m_50cell_0.5deg_flatSplitThresholdExtendVKernelMajCRCRCR_UTM.tif')

################################################################################

panayevsk <- readOGR(dsn = 'Panayevsk enterprise brigades.shp', layer = 'Panayevsk enterprise brigades')
yarsale <- readOGR(dsn = 'Yar-Sale enterprise brigades 2020_correct.shp', layer = 'Yar-Sale enterprise brigades 2020_correct')
ypCorridors <- readOGR(dsn = 'Yar-Sale & Panayevsk brigades areas.shp', layer = 'Yar-Sale & Panayevsk brigades areas')
private <- readOGR(dsn = 'privateHerders_polygon.shp', layer = 'privateHerders_polygon')

################################################################################
#Preprocess brigades

yarsaleDF <- data.frame(yarsale)
yarsaleDF$Number <- gsub(pattern='[a-z].*$','',yarsaleDF$Name)
panayevskDF <- data.frame(panayevsk)
panayevskDF$Number <- gsub(pattern='[a-z].*$','',panayevskDF$Name)

ypCorridorsDF <- data.frame(ypCorridors)
ypCorridorsDF$Enterprise <- gsub(pattern=(' [[:digit:]]*$'),'',ypCorridorsDF$Brigade)
ypCorridorsDF$Number <- as.numeric(gsub(pattern=('([[:alpha:]]|-)* '),'',ypCorridorsDF$Brigade))

for (i in 1:nrow(ypCorridorsDF)){
  if (ypCorridorsDF$Enterprise[i] == 'Yar-Sale'){
    ypCorridorsDF$Total[i] <- yarsaleDF$R..total[yarsaleDF$Number==ypCorridorsDF$Number[i]]
  }
  else {
    ypCorridorsDF$Total[i] <- panayevskDF$R..total[panayevskDF$Number==ypCorridorsDF$Number[i]]
  }
}
ypCorridors$Total <- ypCorridorsDF$Total

ypCorridors$Date[ypCorridors$Date=='VII-VIIII'] <- 'VII-VIII' #Sasha Typo

ypCorridors_May <- ypCorridors[ypCorridors$Date=='V',]
ypCorridors_Jun <- ypCorridors[ypCorridors$Date=='VI',]
ypCorridors_Jul <- ypCorridors[(ypCorridors$Date=='VII'|ypCorridors$Date=='VII-VIII'),]
ypCorridors_Aug <- ypCorridors[(ypCorridors$Date=='VIII'|ypCorridors$Date=='VII-VIII'),]
shapefile(ypCorridors_May,'ypBrigadeAreasMay.shp',overwrite=TRUE)
shapefile(ypCorridors_Jun,'ypBrigadeAreasJune.shp',overwrite=TRUE)
shapefile(ypCorridors_Jul,'ypBrigadeAreasJuly.shp',overwrite=TRUE)
shapefile(ypCorridors_Aug,'ypBrigadeAreasAugust.shp',overwrite=TRUE)

###############################################################################
#Preprocess private
private$Name <- as.character(private$Name)
private$Name[41] <- '20c'
JuneOpts <- c('June-August','June-mid July','June','June-July','All the year')
JulyOpts <- c('June-August','June-mid July','June-July','All the year','July-August','July','mid July-August','July_August')
AugustOpts <- c('June-August','All the year','July-August','mid July-August','August','July_August')

privateJun <- private[private$Time%in%JuneOpts,]
privateJul <- private[private$Time%in%JulyOpts,]
privateAug <- private[private$Time%in%AugustOpts,]
shapefile(privateJun, 'privateHerders_June.shp',overwrite=TRUE)
shapefile(privateJul, 'privateHerders_July.shp',overwrite=TRUE)
shapefile(privateAug, 'privateHerders_August.shp',overwrite=TRUE)

################################################################################

lapply(c('June','July','August'),function(x){
  brigade <- readOGR(dsn = paste0('ypBrigadeAreas',x,'.shp'), layer = paste0('ypBrigadeAreas',x))
  private <- readOGR(dsn = paste0('privateHerders_',x,'.shp'), layer = paste0('privateHerders_',x))
  
  names(brigade)[4] <- "Reindeer" #Change Total to Reindeer like in private
  all <- bind(brigade, private)
  shapefile(all, paste0('Monthly/',x,'_AllReindeer.shp'))
  # all_utm <- spTransform(all, crs(vsi))
  # all_rast <- rasterize(all_utm,vsi,field='Reindeer',fun='sum')
  # writeRaster(all_rast, paste0('Monthly/',x,'_AllReindeer.tif'))
  # 
})
# june_rast2 <- june_rast
# june_rast2[is.na(june_rast2)] <- 0
# kernel <- gaussian.kernel(sigma = 2, n = 25)
# june_smooth <- focal(june_rast2, w=kernel, na.rm=TRUE)
################################################################################

setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing/Reindeer/Monthly")
vsi <- raster('../../MODIS/withTaz/VSI/Monthly_1km/SensTotalW.tif')
lapply(c('June','July','August'),function(x){
  shp <- shapefile(paste0(x,'_AllReindeer.shp'))
  
  all_utm <- spTransform(shp, crs(vsi))
  all_rast <- rasterize(all_utm,vsi,field='Reindeer',fun='sum')
  writeRaster(all_rast, paste0(x,'_AllReindeer.tif'))
})

################################################################################
library(exactextractr)

setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing/Reindeer/Monthly/Merged")
vsi <- raster('../../../MODIS/withTaz/VSI/Monthly_1km/SensTotalW.tif')
landcover <- raster('../../../Pangaea/Landcover/ZAM_LCP_LANDC_SEN12_V02_20150815_20180830_T01.tif')

lapply(c('June','July','August'),function(x){
  shp <- shapefile(paste0(x,'_AllReindeer.shp'))
  
  orig_poly_veg <- exact_extract(landcover, st_as_sf(shp)) #Find land pixels
  land_area <- sapply(orig_poly_veg, function(y){
    sum(y$coverage_fraction[y$value %in% 1:18])*res(landcover)[1]*res(landcover)[2]/10^6 #1:18 are land (not water) cover types
  })
  density <- shp$Reindeer/land_area
  
  poly_union <- raster::union(shp)
  
  veg_again <- exact_extract(landcover, st_as_sf(poly_union))
  union_land_area <- sapply(veg_again, function(y){
    sum(y$coverage_fraction[y$value %in% 1:18])*res(landcover)[1]*res(landcover)[2]/10^6 #1:18 are land (not water) cover types
  })
  
  union_density <- lapply(1:length(poly_union),function(y){
    #The union has one row for each of the new polygons and one column for each of the old polygons.
    #A "1" entry in a column indicates that old polygon comprises the new polygon for that row
    polys_included <- which(data.frame(poly_union[y,1:(ncol(poly_union)-1)]) == 1)
    density <- sum(density[polys_included])
    reindeer <- sum(shp$Reindeer[polys_included])
    
    c(density, reindeer) #Keep reindeer totals for debugging purposes
  })
  union_density <- do.call(rbind, union_density)
  
  new_shp <- poly_union
  new_shp$area <- union_land_area
  new_shp <- new_shp[,ncol(new_shp)] #Just keep the area column
  new_shp$density <- union_density[,1]
  new_shp$reindeer <- union_density[,2]
  
  new_proj <- spTransform(new_shp, crs(vsi))
  new_rast <- rasterize(new_proj,vsi,field='density',fun='sum') #Specifying the function shouldn't matter b/c of union
  new_rast_area <- rasterize(new_proj,vsi,field='area')
  
  writeRaster(new_rast, paste0(x,'_AllReindeerDensity.tif'))
  writeRaster(new_rast_area,paste0(x,'_AllReindeerDensity_LandArea.tif'))
})

shp <- shapefile('JJA_AllReindeer.shp')
orig_poly_veg <- exact_extract(landcover, st_as_sf(shp)) #Find land pixels
land_area <- sapply(orig_poly_veg, function(y){
  sum(y$coverage_fraction[y$value %in% 1:18])*res(landcover)[1]*res(landcover)[2]/10^6 #1:18 are land (not water) cover types
})
shp$area <- land_area
new_proj <- spTransform(shp, crs(vsi))
new_rast <- rasterize(new_proj,vsi,field='rndr_dn') #Specifying the function shouldn't matter b/c of union
new_rast_area <- rasterize(new_proj,vsi,field='area')

writeRaster(new_rast, paste0('JJA','_AllReindeerDensity.tif'))
writeRaster(new_rast_area,paste0('JJA','_AllReindeerDensity_LandArea.tif'))




