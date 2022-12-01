
###################################################################################
# Small functions vsi calc
###################################################################################


# colSds function
colSds<-function(x) apply(x, 2, function(x) sd(x, na.rm=T))

# detrend function
det.func<-function(x) {
	x<-t(x)
	x<- x-rowMeans(x, na.rm=T)	
	t(x)
}

# z score calc function
z.func<-function(x){
	x<-t(x)
	x<-x-rowMeans(x, na.rm=T)	
	x<-x/ apply(x, 1, function(x) sd(x, na.rm=T))
	t(x)
}

###################################################################################
# Main functions for vsi calc
###################################################################################


vsi <-function(evi.ts, temp.ts, cld.ts, aet.ts, eviThresh=1000){

# Define variable names required in analysis
	years <- c("2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014","2015","2016","2017","2018","2019")
	months <- c("jun", "jul", "aug","sep")
 	vars <- c("temp", "cld", "aet")
	
	# Convert temperature to Kelvin	
	#temp.ts[which(temp.ts==-999)] <- NA
	#temp.ts <- (temp.ts + 15000)*0.01
 	temp.ts <- temp.ts/10
		
	# Combine dataset	
	full <- list(evi.ts = evi.ts, temp.ts = temp.ts, cld.ts = cld.ts, aet.ts = aet.ts)
	
	# Replace no data values (-999) with NAs
	#full <- lapply(full,function(x) replace(x,  which(x==-999), NA))
		
	# Build the year-month matrices	
	full <- lapply(full, function(x)matrix(x, nrow = length(years), byrow = TRUE, dimnames = list(years, months)))
	
	# Calculate the climatological means +sds
	means <- lapply(full, function(x)colMeans(x, na.rm=T))
	sds <- lapply(full, function(x)colSds(x))
	
	# Calculate detrended matrices
	fullDet <- lapply(full, function(x)det.func(x))
	fullZ <- lapply(full, function(x)z.func(x))
	
	# # Make t+1 variable	
	# t <- as.vector(t(fullZ$evi.ts))
	# t1 <- c(NA, t[-length(t)])
	# fullZ$t1 <- matrix(t1, nrow=length(years), byrow = TRUE, dimnames=list(years, months))

      ##################
	#  masking step  #
	##################

	# Find the months which should be masked (due to mean temperature being less than 0 and mean evi < 1000
	
	maskTemp <- which(means$temp.ts <= 273)
	maskEvi <- which(means$evi.ts <= eviThresh) 
	
	maskMonths<-c(maskEvi, maskTemp)
	maskMonths<- sort(unique(maskMonths)) # months to be excluded
	
	#	mask<-c(1,2,3,4,5,6,7,8,9,10,11,12)  
	mask<-1:length(months)  
	if(length(maskMonths >0)) mask<- mask[-maskMonths] #months to be included
	
	# Mask the time series	
	full.mask <- lapply(full, function(x)x[,mask])
	fullDet.mask <- lapply(fullDet, function(x)x[,mask])	
	fullZ.mask <- lapply(fullZ, function(x)x[,mask])
	
	# Do an NA check- there are some pixels that flatline in evi or in aet/pet. Want to identify these pixels and remove
	screeningMatrix <- matrix(unlist(fullZ.mask), ncol=4)
	NAcheck<- apply(screeningMatrix, 2, function(x) length(x[is.na(x)==T]))/nrow(screeningMatrix) # Gives you the proportion NAs in each column in the zscores.
	
	NAcheckResult <- length(which(NAcheck == 1)) # Could change the 1 to an argument- then can easily adjust this value in the settings to see the influence of these NAs. Currently removing if full of nas

	#################
	# P C A   R E G #
	#################
  
	# Make objects that are needed to sucessfully run the PCA regression and store the results	 
  pc1<-matrix(NA, nrow=length(years), ncol=length(months), dimnames=list(years, months))
	pc2<-matrix(NA, nrow=length(years), ncol=length(months), dimnames=list(years, months))
	pc3<-matrix(NA, nrow=length(years), ncol=length(months), dimnames=list(years, months))
	#pc4<-matrix(NA, nrow=length(years), ncol=length(months), dimnames=list(years, months))

	PCAloadings<-array(NA, dim=c(3,3,length(months)), dimnames=list(c("1","2", "3"), c(vars), months))
	rsquared<-rep(NA, length(months)); names(rsquared)<-months
	pVal<-rep(NA, length(months)); names(rsquared)<-months
	
	sig.var.inf.sum<-matrix(NA, nrow=3, ncol=length(months), dimnames=list(c(vars), months))
	upCI.sum <- matrix(NA, nrow=3, ncol=length(months), dimnames=list(c(vars), months))
	lowCI.sum <- matrix(NA, nrow=3, ncol=length(months), dimnames=list(c(vars), months))
	 
  if(NAcheckResult == 0) {

	# OK, run the PCA analysis for each month.
		
		PCmask<- seq(1, length(months))
		PCmask[maskMonths]<-NA
			
		for(iii in 1:length(months)){
						
			if(is.na(PCmask[iii])==FALSE){
			
				# select the relevant monthly data
				tempPCA <- fullZ$temp.ts[,iii] 
				cldPCA <- fullZ$cld.ts[,iii]
				aetPCA <- fullZ$aet.ts[,iii]
				#t1PCA <- fullZ$t1[,iii]
				#if(length(t1PCA[is.na(t1PCA)]) == length(years)) t1PCA <- rep(0, length(years))				


				# combine these into a matrix
				varPCA<-cbind(tempPCA, cldPCA, aetPCA); colnames(varPCA)<-c(vars) 
				
				# In one pixel found evi Zscores NaNs because variance was zero in November. Replaced NaNs with 0
				varPCA[is.nan(varPCA)]<-0
					
				# remove rows with NAs in them (required for PCA)
				varPCA = na.omit(varPCA)					
				
				# run the PCA. Only doing a PCA if there are more than 4 complete rows
				if(nrow(varPCA) >5){
				 					
					pca.model<-princomp(varPCA, na.action=na.omit) 
								
					# estimate the loadings of each of the variables on the 3 different PC axes
					PCAloadings[,,iii]<-rbind(pca.model$loadings[,1], pca.model$loadings[,2], pca.model$loadings[,3]) 
					PCAloadings1<-data.frame(rbind(pca.model$loadings[,1], pca.model$loadings[,2], pca.model$loadings[,3]))
				
					# identify which rows of PC tables should be filled
					pcTableMatch = match(row.names(pca.model$scores), row.names(pc1))
			
					# Find the scores for PC1, PC2, PC3, PC4
					pc1[pcTableMatch,iii]<-pca.model$scores[,1] 
					pc2[pcTableMatch,iii]<-pca.model$scores[,2]
					pc3[pcTableMatch,iii]<-pca.model$scores[,3]

					# Regress the PC scores agacldt the evi data
					varTable <- data.frame(resp = fullZ$evi.ts[,iii], PC1 = pc1[,iii], PC2 = pc2[,iii], PC3 = pc3[,iii])
					varTable <- varTable[complete.cases(varTable), ]

					#MPS: lowering numer of complete row requirement from 6 to 4... 
					#might not have an effect still because p value will be too low, 
					#but at least consistent with conditions for PCA above
					if(nrow(varTable)>=6){	
						model <- lm(resp ~ PC1 + PC2 + PC3, data = varTable)	

						# extract the cofficients of the PC regression 
						regCoeff <- model$coefficients[-1]

						# Find the rsquared and pvalue of the PC regression
					 	
						a <- tryCatch((summary(model)$r.squared), warning = function(w) w = NA)
						if(is.na(a)){
							rm(a)
						} else {						

						rsquared[iii]<-summary(model)$r.squared	
						fstat<- summary(model)$fstatistic
						pVal[iii]<-pf(fstat[1], fstat[2], fstat[3], lower.tail=FALSE)
						cInf <- confint(model, c("PC1", "PC2", "PC3"), level = 0.9)
					
						# if a pc axis is signicant at p = 0.1 it will tramsform these back to the original coefficients by mulitplying the loadings of that axis by the regression coefficients
						sigCoef <- which(summary(model)$coefficients[-1,4] <0.1)
						if(length(sigCoef)!=0){							
							sig.var.inf<-PCAloadings1[sigCoef,]*regCoeff[sigCoef]
	 						sig.var.inf.sum[,iii] <- apply(sig.var.inf, 2, sum)
														# Find the upper and lower confidence in the same way
							upCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 2]
							upCI.sum[,iii] <- apply(upCI, 2, sum)

							lowCI <- PCAloadings1[sigCoef,]*cInf[sigCoef, 1]
							lowCI.sum[,iii] <- apply(lowCI, 2, sum)
						
						} # closes the sig if statement	
					
						} # tryCatch if statement
					} # closesnrow masking statement
													
				}	# closes varPCA masking statement
					
			} # closes the NA masking if statement
		} #closes the monthly for loop for running the pcas/ regressions
		
			
	} # closes the NA check masking statement
	
			
	 #############################################
	# # V A R I A N C E   B Y  SIG COEF #
	# #############################################
	monthCoefSel <- which(is.na(sig.var.inf.sum[1,])==FALSE)
	if(length(monthCoefSel) ==0){
		
		sdDetVecSig <- rep(NA, 4)
		meanVecSig <- rep(NA, 4)
		
	} else {
		#fullDet already demeaned by month, now stddev on demeaned values
	  #(prevents inflated stddev due to spread of monthly averages)
	  #meanVecSig is mean over all selected months, not months individually
		detVecSig <- lapply(fullDet, function(x)c(t(x[,monthCoefSel])))
		fullVecSig <- lapply(full, function(x)c(t(x[,monthCoefSel])))
		
		sdDetVecSig<- unlist(lapply(detVecSig, function(x)sd(x, na.rm=T)))
		meanVecSig <- unlist(lapply(fullVecSig, function(x)mean(x, na.rm=T)))
				
	}

	# Store the results output. Processing of these takes place in the mapping script.
	
	results<-list(mask = mask, rsquared = rsquared, sig.var.inf.sum = sig.var.inf.sum, 
	pVal= pVal,  sdDetVecSig = sdDetVecSig, meanVecSig = meanVecSig, upCIsum = upCI.sum, lowCI.sum = lowCI.sum)

	results		
}