##########################################################
# Range function to standardise between 0 and 100
##########################################################

range01 <- function(y){	
	x<-y[which(is.finite(y)==T)]
	((y-min(x))/(max(x)-min(x)))*100
	}



##########################################################
### MAKE A CSV RASTER FUNCTION FOR A VARIABLE
##########################################################

makeRaster <- function(x, landsea, resultPath, var){

  #MPS: I have NA's outside of 66,66,74,74 in my mask (regardless if they are land or water)
  #Let's just make these all sea pixels (landsea=1) to avoid comparing NA's
  landsea[is.na(landsea)] <- 1
	result <- landsea
	result[result == 1] <- -9999  # Make the sea pixels -9999 like in Seddon's mask 
	result[result == 0] <- x  # 0 (not 1) is land in my landsea mask
	resultMatrix <- matrix(result, nrow = 1080, byrow = FALSE)
	resultMatrix[is.na(resultMatrix)] <- -8888  # converts the NA values in the results to -8888 so that they can be distinguished
	write.table(resultMatrix, file = paste(resultPath, var, ".csv", sep="" ), sep= ",", col.names=FALSE, row.names= FALSE)
}

##########################################################
# TRIM THE ANOMALIES FOR A VARIABLE
##########################################################

cutAnomalies <- function(x){
	# Cut the anomolies			
			toPlot <- c(x)
			toPlot1<- toPlot[!is.na(toPlot)]
			toPlot <-toPlot1
			sortedLow<-sort(toPlot)
	 		lowCut<-sortedLow[length(sortedLow)*0.005]
	 		sortedHigh<-sort(toPlot, decreasing=T)
	 		highCut<-sortedHigh[length(sortedHigh)*0.005]
	 		x[which(x>highCut)]<-highCut
	 		x[which(x<lowCut)]<-lowCut
	 		x			
}

##########################################################
# MAKE A CSV MATRIX
# Wrapper function to do the full landSea calculation and standardisation with the range between 0-100, then make the appropriate csv matrices
##########################################################

landSeaProj <- function(resultMatrix, landsea, fullStand = TRUE, resultPath, raw = FALSE){
	
	# resultMatrix = results that need to be projected on the landsea mask 
	# fullStand = Do you want to standardise between 0 and 100 against all the variables (TRUE), or for each variable individually (FALSE)
	# resultPath = path to save the results raster (as .csv file)
	# raw = Do you want to output the raw values when making the .csv matrix?
	
	print("Making projection to landsea matrix")
	progBar <- txtProgressBar(min = 0, max = length(colnames(resultMatrix)), style = 3)
 	
	if (raw == TRUE){
		varName <- colnames(resultMatrix)
		for(i in 1: length(varName)){
			makeRaster(resultMatrix[,i], landsea, resultPath, var= varName[i]) # Make the raster of each variable
			setTxtProgressBar(progBar, i)	
	}
	} else {
		anomCut <- apply(resultMatrix, 2, function(x)cutAnomalies(x)) # trim the anomalies
		if(fullStand == FALSE) stand0_100 <- apply(anomCut, 2, function(x)range01(x)) else stand0_100 <- range01(anomCut) # standardise between 0 and 100
		varName <- colnames(stand0_100)
		for(i in 1: length(varName)) {
			makeRaster(stand0_100[,i], landsea, resultPath, var= varName[i]) # Make the raster of each variable
			setTxtProgressBar(progBar, i)
		}
	}
}



##########################################################
# Function for compiling variance anomalies
##########################################################

compileVarTilesAnon <- function(tilePath, var = "meanVecSig"){

resultMatrix <- matrix(NA, nrow=1, ncol=4)
print(paste("Compiling tiles of the variable", var, sep=" ")) 
progBar <- txtProgressBar(min =0, max = 90, style = 3)

for(tile in 1:90) {
	load(paste(tilePath, "/tile_", tile, ".RData", sep=""))
	namesResults <- names(vsiResult[[1]])
	dimension <- 4
	
	vsi <- vsiResult
	rm(vsiResult)
	xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 3.var
	varNum<-which(namesResults== var) # Identify correct item in list for each pixel
	
	for(i in 1:length(vsi)) {
				toPlot<-vsi[[i]] # extract a given pixel
				varCoef <- toPlot[[varNum]] # extract the coefficient matrix 
				xToPlot[i,] <- varCoef
								
	} # End of pixel forloop
	resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
	setTxtProgressBar(progBar, tile)
}
	resultMatrix <- resultMatrix[-1,] 
	colnames(resultMatrix) <- c("evi", "temp", "cld", "aet") 
	 # All you need to do now is line them up with the landsea mask (see next section)
	close(progBar)
	resultMatrix

}

##########################################################
# Function for estimating variance anomalies
##########################################################

anomCalcNew <- function(
		var = "temp", 
		filePath = paste(mainDir, "/results/anom", sep= "")) {
	
	colSel <- match(var, colnames(meanSig))

	x = meanSig[,colSel]
	y = sdSig[,colSel]
	
	# Compile the data to make the the plot
	require("nlme")
	
	# Residual result is the variable that you want map
	residualResult<-rep(NA, length(x))
	I <- which(is.na(x) == FALSE)
	
	nls2 <- nls(y~b0+b1*(x)+ b2*x^2,start=c(b0=0,b1=0, b2=0) )
			
	colList <- list(evi = rgb(100,100 ,100, 50, maxColorValue = 255), temp =  rgb(100, 0 ,0, 50, maxColorValue = 255),
			    aet = rgb(0, 0 , 100, 50, maxColorValue = 255), cld = rgb(0, 100 ,0, 50, maxColorValue = 255))
	
	smp <- sample(I, 1000)

	plot(x[smp], y[smp], col= eval(parse(text = paste("colList$", var, sep=""))), pch=16, xlab = "Mean", ylab = "Variance", main = var)
	
	xPred <- seq(floor(min(x, na.rm= TRUE)), ceiling(max(x, na.rm= TRUE)))
	yPred <- predict(nls2, newdata = data.frame(x= xPred))

	lines(xPred, yPred)
		
	residualResult[I] <-residuals(nls2)
     	anomCut <-cutAnomalies(residualResult) # trim the anomalies
	stand0_100 <- range01(anomCut)
	
	makeRaster(x= stand0_100, landsea=landsea, resultPath= filePath, var)
	
	result <- list(xPred = xPred, yPred = yPred, x = x[smp], y = y[smp], varName = var) 
	result	
}





##########################################################
######### Compiles coefficient variables 
##########################################################

compileCoefTiles <- function(tilePath, absVal = TRUE){

resultMatrix <- matrix(NA, nrow=1, ncol=3)
print(paste("Compiling coefficient tiles"))
progBar <- txtProgressBar(min = 0, max = 90, style = 3)

for(tile in 1:90) {
	load(paste(tilePath, "/tile_", tile, ".RData", sep=""))
	
	vsi<-vsiResult
	rm(vsiResult)
	namesResults <- names(vsi[[1]])
	var <- "sig.var.inf.sum"
	dimension <- 3
	
	xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 4.var
	varNum<-which(namesResults== var) # Identify correct item in list for each pixel
	
	for(i in 1:length(vsi)) {
				toPlot<-vsi[[i]] # extract a given pixel
				varCoef <- toPlot[[varNum]] # extract the coefficient matrix 
				
				if(absVal== TRUE){
				xToPlot[i,] <- apply(abs(varCoef), 1, function(x)mean(x, na.rm=TRUE)) 
				} else{
				xToPlot[i,] <- apply((varCoef), 1, function(x)mean(x, na.rm=TRUE)) # if don't want to use absolute coefficient values
					}				
	} # End of pixel forloop
	resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
setTxtProgressBar(progBar, tile)
}
resultMatrix <- resultMatrix[-1,] 
colnames(resultMatrix) <- rownames(varCoef) 
resultMatrix # All you need to do now is line them up with the landsea mask (see next section)
}


##########################################################
####### Compiles coefficient CI 
##########################################################

compileCITiles <- function(tilePath){

resultMatrix <- matrix(NA, nrow=1, ncol=3)
print(paste("Compiling coefficient CI tiles"))
progBar <- txtProgressBar(min = 0, max = 90, style = 3)

for(tile in 1:90) {
	load(paste(tilePath, "/tile_", tile, ".RData", sep=""))
	
	vsi<-vsiResult
	rm(vsiResult)
	namesResults <- names(vsi[[1]])
	up <- "upCIsum"
	down <- "lowCI.sum"
	dimension <- 3
	
	xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 3.var
	upNum<-which(namesResults== up) # Identify correct item in list for each pixel
	downNum<-which(namesResults== down)
	
	for(i in 1:length(vsi)) {
				toPlot<-vsi[[i]] # extract a given pixel
				upCI <- toPlot[[upNum]] # extract the coefficient matrix 
				downCI <- toPlot[[downNum]]				
				
				range <- abs(upCI- downCI)
				rangeAv <- apply(range, 1, function(x)mean(x, na.rm=TRUE))
				xToPlot[i,] <- rangeAv				
	} # End of pixel forloop
	resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
	setTxtProgressBar(progBar, tile)
}
resultMatrix <- resultMatrix[-1,] 
colnames(resultMatrix) <- rownames(upCI) 
resultMatrix # All you need to do now is line them up with the landsea mask (see next section)
}






##########################################################
####### Compiles coefficient variables without using t-1
##########################################################

compileCoefTilesNONT1 <- function(tilePath, absVal = TRUE){

resultMatrix <- matrix(NA, nrow=1, ncol=3)
for(tile in 1:90) {
	load(paste(tilePath, "/tile_", tile, ".RData", sep=""))
	
	vsi<-vsiResult
	rm(vsiResult)
	namesResults <- names(vsi[[1]])
	var <- "sig.var.inf.sum"
	dimension <- 3
	
	xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 3.var
	varNum<-which(namesResults== var) # Identify correct item in list for each pixel
	
	for(i in 1:length(vsi)) {
				toPlot<-vsi[[i]] # extract a given pixel
				varCoef <- toPlot[[varNum]][-4,] # extract the coefficient matrix 
				
				if(absVal== TRUE){
				xToPlot[i,] <- apply(abs(varCoef), 1, function(x)mean(x, na.rm=TRUE)) # NOTE USING THE ABSOLUTE COEFFICIENT VALUES HERE
				} else{
				xToPlot[i,] <- apply((varCoef), 1, function(x)mean(x, na.rm=TRUE)) # if don't want to use absolute coefficient values
					}				
	} # End of pixel forloop
	resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
}
resultMatrix <- resultMatrix[-1,] 
colnames(resultMatrix) <- rownames(varCoef) 
resultMatrix # All you need to do now is line them up with the landsea mask (see next section)
}



##########################################################
######## Compiles coefficient variables as a percentage
##########################################################

compileCoefTilesPERCENTAGE <- function(tilePath, absVal = TRUE){

resultMatrix <- matrix(NA, nrow=1, ncol=3)

for(tile in 1:90) {
	load(paste(tilePath, "/tile_", tile, ".RData", sep=""))
	
	vsi<-vsiResult
	rm(vsiResult)
	namesResults <- names(vsi[[1]])
	var <- "sig.var.inf.sum"
	dimension <- 3
	
	xToPlot <- matrix(NA, nrow=length(vsi), ncol= dimension) # Create your results matrix, in this instance will be n.land-pixels x 3.var
	varNum<-which(namesResults== var) # Identify correct item in list for each pixel
	
	for(i in 1:length(vsi)) {
				toPlot<-vsi[[i]] # extract a given pixel
				varCoef <- toPlot[[varNum]] # extract the coefficient matrix 
				
				
				result <- apply(abs(varCoef), 1, function(x)mean(x, na.rm=TRUE)) # NOTE USING THE ABSOLUTE COEFFICIENT VALUES HERE
				resultSum <- sum(result)
				xToPlot[i,]	<- 	(result/resultSum)*100		

						
	} # End of pixel forloop
	resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
}
resultMatrix <- resultMatrix[-1,] 
colnames(resultMatrix) <- rownames(varCoef) 
resultMatrix # All you need to do now is line them up with the landsea mask (see next section)
}








##########################################################
######## FUNCTIONS FOR ESTIMATING THE NUMBER OF MONTHS WITH SIG. CLIMATE COEF.
##########################################################

compileSigMonthTiles <- function(tilePath, var = "sig.var.inf.sum"){

months <- c("129","145","161","177","193","209","225","241","257")

resultMatrix <- matrix(NA, nrow=1, ncol=length(months))

for(tile in 1:90) {
	load(paste(tilePath, "/tile_", tile, ".RData", sep=""))
	
	namesResults <- names(vsiResult[[1]])
	
	xToPlot <- matrix(NA, nrow=length(vsiResult), ncol=length(months)) # Create your results matrix
	varNum<-which(namesResults== var) # Identify correct item in list for each pixel
	
	for(i in 1:length(vsiResult)) {
				toPlot<-vsiResult[[i]] # extract a given pixel
				varInfSum <- toPlot[[varNum]][1,]
				varInfSum[!is.na(varInfSum)] <- 1 				
				xToPlot[i,] <- varInfSum
	
	} # End of pixel forloop
	resultMatrix <- rbind(resultMatrix, xToPlot) # Combine all the results from the different variables
}
resultMatrix <- as.matrix(resultMatrix[-1,] , ncol=1)
#colnames(resultMatrix) <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec") 
colnames(resultMatrix) <- months
rowSums(resultMatrix, na.rm= TRUE) # All you need to do now is line them up with the landsea mask (see next section)
}

##### TURN THE MONTH SIG INTO RASTER

landSeaMonthProj <- function(resultMatrix, landsea, resultPath){
	
	varName <- colnames(resultMatrix)
	for(i in 1: length(varName)) makeRaster(resultMatrix[,i], landsea, resultPath, var= varName[i]) # Make the raster of each variable
}




