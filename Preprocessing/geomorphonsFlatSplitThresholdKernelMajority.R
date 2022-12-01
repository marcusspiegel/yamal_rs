library(raster)
library(parallel)

setwd("E:/Data/Biogeosciences/Yamal/ArcticDEM/landforms/")

geomorphons <- raster("Geomorphons_withTaz_32m_50cell_0.5deg.tif")
dem <- raster("../mosaic/ArcticDEM_32m_withTaz.tif")

fullMatrix <- matrix(NA, nrow=nrow(geomorphons), ncol=ncol(geomorphons))
gMatrix <- as.matrix(geomorphons)
demMatrix <- as.matrix(dem)

lowlands <- 11
uplands <- 0

#First, let's threshold by dem values
gMatrix[(demMatrix < 0) & (gMatrix == 1)] <- lowlands
rm(demMatrix)
outRaster <- raster(gMatrix, template = geomorphons)
writeRaster(outRaster,'2019-07-26/Geomorphons_yamal_32m_50cell_0.5deg_flatSplitThreshold.tif')

#row by row
progBar <- txtProgressBar(min = 0, max = nrow(geomorphons), style = 3)
for (x in 1:nrow(geomorphons)){
  currRow <- gMatrix[x,]
  
  #First find all the stretches of flat pixels
  flats <- which(currRow == 1)
  setTxtProgressBar(progBar, x)
  if (length(flats)==0){
    fullMatrix[x,] <- currRow
    next
  }
  flatStarts <- c()
  flatStartsInds <- c()
  flatEnds <- c()
  flatEndsInds <- c()
  inFlat <- FALSE
  j <- 1
  
  for (i in 1:length(flats)) {
    if (i == length(flats)){ #Last flat pixel corner case
      if (inFlat == FALSE){
        flatStarts[j] <- flats[i]
        flatStartsInds[j] <- i
      }
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      break
    } else if (inFlat == FALSE) {
      flatStarts[j] <- flats[i]
      flatStartsInds[j] <- i
      inFlat <- TRUE
    }
    if ((flats[i+1]-flats[i])>5){ #This can be changed to filter lone pixels
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      j <- j+1
      inFlat <- FALSE
    }
  }
  
  #Then for each stretch, check the pixels on either end
  for (i in 1:length(flatStarts)){
    #Corner cases to ignore flat pixels
    if ((flatStarts[i] == 1)||(flatEnds[i] == length(currRow))){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    
    end1 <- 0
    end2 <- 0
    
    #left end of flat stretch
    k <- flatStarts[i]-1
    if (is.na(currRow[k])){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    #Already adjacent to green pixels from thresholding
    if (currRow[k] == lowlands){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    while (currRow[k] < 8){
      k <- k-1
      end1 <- end1 + 1
      if ((k == 0) || (is.na(currRow[k])))
        break
    }
    
    #right end of flat stretch
    k <- flatEnds[i]+1
    if (is.na(currRow[k])){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    #Already adjacent to green pixels from thresholding
    if (currRow[k] == lowlands){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    while (currRow[k] < 8){
      k <- k+1
      end2 <- end2 + 1
      if (is.na(currRow[k]))
        break
    }
    
    #end1 <- currRow[flatStarts[i]-1]
    #end2 <- currRow[flatEnds[i]+1]
    #if (is.na(end1)||is.na(end2))
    #  next
    if ((end1>=5)&&(end2>=5))
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- uplands
    else
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
    
  }
  fullMatrix[x,] <- currRow
}

#col by col
fullMatrixV <- matrix(NA, nrow=nrow(geomorphons), ncol=ncol(geomorphons))
progBar <- txtProgressBar(min = 0, max = nrow(geomorphons), style = 3)
for (x in 1:ncol(geomorphons)){
  currRow <- gMatrix[,x]
  
  #First find all the stretches of flat pixels
  flats <- which(currRow == 1)
  setTxtProgressBar(progBar, x)
  if (length(flats)==0){
    fullMatrixV[,x] <- currRow
    next
  }
  flatStarts <- c()
  flatStartsInds <- c()
  flatEnds <- c()
  flatEndsInds <- c()
  inFlat <- FALSE
  j <- 1
  
  for (i in 1:length(flats)) {
    if (i == length(flats)){ #Last flat pixel corner case
      if (inFlat == FALSE){
        flatStarts[j] <- flats[i]
        flatStartsInds[j] <- i
      }
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      break
    } else if (inFlat == FALSE) {
      flatStarts[j] <- flats[i]
      flatStartsInds[j] <- i
      inFlat <- TRUE
    }
    if ((flats[i+1]-flats[i])>5){ #This can be changed to filter lone pixels
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      j <- j+1
      inFlat <- FALSE
    }
  }
  
  #Then for each stretch, check the pixels on either end
  for (i in 1:length(flatStarts)){
    #Corner cases to ignore flat pixels
    if ((flatStarts[i] == 1)||(flatEnds[i] == length(currRow))){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    
    end1 <- 0
    end2 <- 0
    
    #left end of flat stretch
    k <- flatStarts[i]-1
    if (is.na(currRow[k])){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    #Already adjacent to green pixels from thresholding
    if (currRow[k] == lowlands){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    while (currRow[k] < 8){
      k <- k-1
      end1 <- end1 + 1
      if ((k == 0) || (is.na(currRow[k])))
        break
    }
    
    #right end of flat stretch
    k <- flatEnds[i]+1
    if (is.na(currRow[k])){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    #Already adjacent to green pixels from thresholding
    if (currRow[k] == lowlands){
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
      next
    }
    while (currRow[k] < 8){
      k <- k+1
      end2 <- end2 + 1
      if (is.na(currRow[k]))
        break
    }
    
    #end1 <- currRow[flatStarts[i]-1]
    #end2 <- currRow[flatEnds[i]+1]
    #if (is.na(end1)||is.na(end2))
    #  next
    if ((end1>=5)&&(end2>=5))
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- uplands
    else
      currRow[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
    
  }
  fullMatrixV[,x] <- currRow
}

#Combine the two full matrices to prioritise highlands (uplands) 
fullMatrix[fullMatrixV == uplands] <- uplands


outRaster <- raster(fullMatrix, template = geomorphons)
writeRaster(outRaster,'2019-07-26/Geomorphons_yamal_32m_50cell_0.5deg_flatSplitThresholdExtendV.tif')


###############################################################################
#Let's try a majority kernel

Mode <- function(x) {
  x <- x[!is.na(x)]
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[(tab == max(tab))]
}

ksize <- 5 #11x11 kernel
#fullMatrix <- as.matrix(raster('2019-07-16/Geomorphons_yamal_32m_50cell_flatSplitThresholdExtendV.tif'))
flatMatrix <- fullMatrix
flatMatrix[(fullMatrix > uplands) & (fullMatrix < lowlands)] <- NA

ntiles <- 125
rowsPerTile <- nrow(flatMatrix)/ntiles

dir.create("2019-07-26/kernelTiles")

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
parLapply(cl, 1:ntiles, function(i){
  rowStart <- max(((i-1)*rowsPerTile+1-ksize),1)
  rowEnd <- min((i*rowsPerTile+ksize),nrow(fullMatrix))
  rows <- rowStart:rowEnd
  #progBar <- txtProgressBar(min = 0, max = nrow(fullMatrix), style = 3)
  outMatrix <- matrix(NA, nrow=length(rows), ncol=ncol(fullMatrix))
  for (x in 1:length(rows)){
    row <- rows[x]
    for (y in 1:ncol(fullMatrix)) {
      if (!is.na(flatMatrix[row,y])){
        lowX <- max(1,(row-ksize))
        upX <- min(nrow(flatMatrix),(row+ksize))
        lowY <- max(1,(y-ksize))
        upY <- min(ncol(flatMatrix),(y+ksize))
        kernel <- flatMatrix[lowX:upX,lowY:upY]
        kmode <- Mode(kernel)
        if (length(kmode)>1)
          kmode <- uplands
        outMatrix[x,y]<-kmode
      }
    }
    #setTxtProgressBar(progBar, x)
  }
  saveRDS(outMatrix, file = paste0("2019-07-26/kernelTiles/tile",i,".rds"))
})
stopCluster(cl)
gc()

tileMatrix <- matrix(NA, nrow=nrow(fullMatrix), ncol=ncol(fullMatrix))

progBar <- txtProgressBar(min = 0, max = ntiles, style = 3)
tile <- readRDS("2019-07-26/kernelTiles/tile1.rds")
tileMatrix[1:rowsPerTile,] <- tile[1:rowsPerTile,]
for (i in 2:(ntiles)){
  tile <- readRDS(paste0("2019-07-26/kernelTiles/tile", i, ".rds"))
  rowStart <- (i-1)*rowsPerTile+1
  rowEnd <- i*rowsPerTile
  tileMatrix[rowStart:rowEnd,] <- tile[(ksize+1):(ksize+rowsPerTile),]
  setTxtProgressBar(progBar, i)
}


fullMatrix[!is.na(tileMatrix)]<-tileMatrix[!is.na(tileMatrix)]

fullRast <- raster(fullMatrix, template=geomorphons)
writeRaster(fullRast, '2019-07-26/Geomorphons_yamal_32m_50cell_0.5deg_flatSplitThresholdExtendVKernel.tif')


###############################################################################

#Next, reassign split values as the majority from each cluster of flat pixels
#in each column, then each row

#fullMatrix <- as.matrix(raster('2019-07-19/Geomorphons_yamal_32m_50cell_0.5deg_flatSplitThresholdExtendVKernelMajCRCR.tif'))


#columns
progBar <- txtProgressBar(min = 0, max = ncol(fullMatrix), style = 3)
for (x in 1:ncol(fullMatrix)){
  curr <- fullMatrix[,x]
  
  #First find all the stretches of flat pixels
  flats <- which((curr == uplands) | (curr == lowlands))
  if (length(flats)==0)
    next
  flatStarts <- c()
  flatStartsInds <- c()
  flatEnds <- c()
  flatEndsInds <- c()
  inFlat <- FALSE
  j <- 1
  
  for (i in 1:length(flats)) {
    if (i == length(flats)){
      if (inFlat == FALSE){
        flatStarts[j] <- flats[i]
        flatStartsInds[j] <- i
      }
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      break
    } else if (inFlat == FALSE) {
      flatStarts[j] <- flats[i]
      flatStartsInds[j] <- i
      inFlat <- TRUE
    } 
    if ((flats[i+1]-flats[i])>5){ #This can be changed to filter lone pixels
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      j <- j+1
      inFlat <- FALSE
    }
  }
  
  #Then for each stretch, make the stretch the majority of the pixels
  for (i in 1:length(flatStarts)){
    meanFlats <- mean(curr[flats[flatStartsInds[i]:flatEndsInds[i]]])
    if (meanFlats > mean(c(uplands, lowlands)))
      curr[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
    else
      curr[flats[flatStartsInds[i]:flatEndsInds[i]]] <- uplands
  }
  fullMatrix[,x] <- curr
  setTxtProgressBar(progBar, x)
}

#rows
progBar <- txtProgressBar(min = 0, max = nrow(fullMatrix), style = 3)
for (x in 1:nrow(fullMatrix)){
  curr <- fullMatrix[x,]
  
  #First find all the stretches of flat pixels
  flats <- which((curr == uplands) | (curr == lowlands))
  if (length(flats)==0)
    next
  flatStarts <- c()
  flatStartsInds <- c()
  flatEnds <- c()
  flatEndsInds <- c()
  inFlat <- FALSE
  j <- 1
  
  for (i in 1:length(flats)) {
    if (i == length(flats)){
      if (inFlat == FALSE){
        flatStarts[j] <- flats[i]
        flatStartsInds[j] <- i
      }
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      break
    } else if (inFlat == FALSE) {
      flatStarts[j] <- flats[i]
      flatStartsInds[j] <- i
      inFlat <- TRUE
    } 
    if ((flats[i+1]-flats[i])>5){ #This can be changed to filter lone pixels
      flatEnds[j] <- flats[i]
      flatEndsInds[j] <- i
      j <- j+1
      inFlat <- FALSE
    }
  }
  
  #Then for each stretch, make the stretch the majority of the pixels
  for (i in 1:length(flatStarts)){
    meanFlats <- mean(curr[flats[flatStartsInds[i]:flatEndsInds[i]]])
    if (meanFlats > mean(c(uplands, lowlands)))
      curr[flats[flatStartsInds[i]:flatEndsInds[i]]] <- lowlands
    else
      curr[flats[flatStartsInds[i]:flatEndsInds[i]]] <- uplands
  }
  fullMatrix[x,] <- curr
  setTxtProgressBar(progBar, x)
}

fullRast <- raster(fullMatrix, template=geomorphons)
writeRaster(fullRast, '2019-07-26/Geomorphons_yamal_32m_50cell_0.5deg_flatSplitThresholdExtendVKernelMajCRCRCR.tif')
