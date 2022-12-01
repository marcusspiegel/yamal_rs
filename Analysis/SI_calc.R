library(raster)
library(rgdal)

setwd("/lustre/soge1/projects/BiogeosciencesLab/Yamal_RemoteSensing")

yp <- readOGR(dsn = 'Brigades/yp.shp', layer = 'yp')
taz <- readOGR(dsn = 'Brigades/taz.shp', layer = 'taz')

chelsa <- list.files('MODIS/withTaz/Projected/Monthly_1km/tif/Temperature',full.names = T)

yp <- spTransform(yp, crs(raster(chelsa[[1]])))
taz <- spTransform(taz, crs(raster(chelsa[[1]])))

out <- lapply(chelsa, function(x){
  r <- raster(x)
  r_yp <- mask(r,yp)
  r_taz <- mask(r,taz)
  c(cellStats(r_yp,mean,na.rm=T),cellStats(r_taz,mean,na.rm=T))
})

temps <- do.call(rbind, out)
temps2 <- temps[-seq(4,80,4),]
june_means <- mean(temps[seq(1,80,4),2])/10-273.15