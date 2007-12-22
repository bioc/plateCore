# utilities.R
# 
# TODO: Add comment
#
# Author: straine
###############################################################################



makePhenoData <- function(plateDesc,abName="Ab.Name",sampleType="Sample.Type",negCon="Negative.Control",plateMetaData=NULL,...) {
	
	##--------------------------------------------------------------------------
	##
	## Function: makePhenoData.R
	##
	## Description: Create an AnnotatedDataFrame for a flowSet from a data.frame.    
	##
	##--------------------------------------------------------------------------
	
	
	## Add some validity check function
	
	
	## Now get the unique channels, and remove dashes make them legal column names
	chans <- unique(plateDesc$Channel)
	chans <- gsub("-",".",chans[chans != ""])
	plateDesc$Channel <- gsub("-",".",plateDesc$Channel)
	
	## Make a plate layout data.frame
	plateLayout <- data.frame(Well.Id=unique(plateDesc$Well.Id))
	
	for(wellChan in chans) {
		temp <- subset(plateDesc,Channel==wellChan,select=c("Well.Id",abName,sampleType,negCon))
		colnames(temp) <- c("Well.Id",wellChan,paste(sampleType,wellChan,sep="."),paste(negCon,wellChan,sep="."))
		plateLayout <- merge(plateLayout,temp,by="Well.Id",all.x=TRUE)
	}
	
	
	temp <- subset(plateDesc,Channel=="",select=c("Well.Id",sampleType))
	plateLayout <- merge(plateLayout,temp,by="Well.Id",all.x=TRUE)
	
	plateLayout <- cbind(name=plateLayout$Well.Id,plateLayout)
	
	if(is.null(plateMetaData)) {
		plateMetaData <- data.frame(colnames(plateLayout))
		rownames(plateMetaData) <- colnames(plateLayout)
	}
	
	tempPheno.adf <- new("AnnotatedDataFrame", data=plateLayout, varMetadata=plateMetaData, dimLabels=c("rowNames", "colNames"))
	sampleNames(tempPheno.adf) <- plateLayout$name
	return(tempPheno.adf)
}


#########################################################################################################
##
## Function: callPops
##
## Description:  Use mclust to say whether there are 1 or 2 populations in a well/channel.  
##	The function uses the unstained wells to call the minimum quantile.  Different analyses have
##	differing levels of garbage near the axes, need to cut this off since it complicates determining
##	the number of populations.
##
## Note: This should not be used when there are multiple populations in the unstained wells.
##
#########################################################################################################

callPops <- function(data,chanCols,bySize=0.01,cutOff=1.1,minDiff=0.25) {
	
	require(mclust)
	options(warn=-1)
	
	unstain <- unlist(subset(pData(phenoData(data)),as.logical((Sample.Type=="Unstained") ),select="name"))
	
	pops <-  rep(2,length=length(chanCols))
	minQuant = vector("numeric",length=length(chanCols))
	names(minQuant) = chanCols
	
	## Set minQuant at the value which gives one population in the unstained samples
	while(mean(pops)>cutOff) {
		temp.df <- fsApply(data[unstain], function(x) {
					sapply(chanCols, function(y) {
								quantVals <- quantile(log10(exprs(x)[,y]),probs=seq(minQuant[y],0.9999,by=bySize))
								quantVals <- quantVals[quantVals>0]
								testClust <- try(Mclust(quantVals,G=1:2,modelNames=c("V")),silent=TRUE)	
								length(testClust$parameters$mean)
							})
				})
		pops <- apply(temp.df,2,mean)
		for(i in 1:length(pops)) {
			if(pops[i]>1) {minQuant[i] <- minQuant[i] + 0.01}
		}	
	}
	
	temp.df <- fsApply(data[unstain], function(x) {
				sapply(chanCols, function(y) {
							quantile(log10(exprs(x)[,y]),probs=minQuant[y])
						})
			})
	
	lowerBound <- apply(temp.df,2,mean)
	
	names(lowerBound) <- chanCols
	
	## Then run it again for all the data.
	temp.df <- fsApply(data, function(x) {
				sapply(chanCols, function(y) {
							quantVals <- quantile(log10(exprs(x)[,y]),probs=seq(0,0.9999,by=bySize))  
							if(sum(quantVals>lowerBound[y])<=2) {
								l <- length(quantVals)
								quantVals <- quantVals[c(l-2:l)]
							} else { 
								quantVals <- quantVals[quantVals>lowerBound[y]] 
							}
							testClust <- Mclust(quantVals,G=1:2,modelNames=c("V"))
							numPops <- length(testClust$parameters$mean)
							
							## I've tried doing t-tests for the difference in means, but it's picking up too much
							## Just using the ad-hoc meansDiff for now
							if(numPops == 1) {
								return(1) 
							} else {	
								meansDiff <- (testClust$parameters$mean[2] - testClust$parameters$mean[1])/mean(testClust$parameters$mean)
								if(meansDiff > minDiff ) return(2)
								return(1)	
							}
							
						})
			})
	
	temp.df
	
}


#########################################################################################################
##
## Function: replaceColHeadings
##
## Description:  FCS 2.0 have FL*-H column names, which FCS 3.0 have FITC.A,...  
##
#########################################################################################################

replaceColHeadings <- function(plateSet,orig=c("FL1-H","FL2-H","FL3-H","FL4-H"),
		replace=c("FITC.A","PE.A","PerCP.Cy5.5.A","APC.A")) {
	
	colHeadings <- exprs(plateSet[[1]])
	
	for(i in 1:length(orig)) {
		colHeadings[colHeadings==orig[[i]]] <- replace[[i]]
	}
	
	fsApply(plateSet,function(x) {
				colnames(exprs(x)) <- colHeadings
			})
	
}
#########################################################################################################
##
## Function: eventCheck
##
## Description:  Get the number of events in each well.  I anticipate modifying this for a more thorough
##	evaluation later.  
##
#########################################################################################################
eventCheck <- function(data,chan="FSC.A",minSignal=0) {
	
	fsApply(data,function(x) {
				sum(exprs(x)[,chan]>minSignal)
			})
	
}

##########################################################################################################
###
### Function: chanReg
###
### Description:  Fit linear model, regress chan2 on chan1.  
###
##########################################################################################################
#
#chanReg <- function(chan1,chan2,data) {
#	
#	require(rrcov)
#	
#	coeffs <- ltsReg(as.formula(paste(chan2," ~ ",chan1)), data.frame(exprs(data)))$coefficients
#
#	return(round(coeffs,digits=3))
#}

#########################################################################################################
##
## Function: fixAutoFl
##
## Description:  Apply Florian Hahne's correction for autofluorescence for the channels in chanCols.  The
##	approach fits a linear model of the form FL*-H ~ FSC.  The correction parameter is the mean of the slopes
## 	from the unstained samples in a plate.  
##
##  Sets the median to be the same after the correction.
##
#########################################################################################################

setMethod("fixAutoFl",signature("flowSet"), 
		function(fp,chanCols,unstain,fsc="FSC.A") {
			
			require(rrcov)
			
			## Identify the wells containing unstained samples
			if(missing(unstain)) {
				unstainWells <- unlist(subset(pData(phenoData(fp)),as.logical((Sample.Type=="Unstained") ),select="name"))
			} else {
				unstainWells <- unstain
			}
			
			## Fit a linear model to the unstained data, get the slope
			unstainFits <- fsApply(fp[unstainWells], function(x) {
						unlist(lapply(chanCols,function(y) {
											exprData <- exprs(x)[(exprs(x)[,y]> 1),]	
											ltsReg(as.formula(paste(y," ~ ",fsc,sep="")), log(data.frame(exprData)))$coefficients[fsc]
										}))
					})	
			
			## Get a matrix of mean slopes
			coeff.mat <- matrix(apply(unstainFits,2,mean),nrow=1)
			
			## Apply the correction to the channels of interest in chanCols\
			## Truncate values less than 0 at 0
			fp <- fsApply(fp,function(x) {
						initMFIs <- log(apply(exprs(x)[,chanCols],2,median))	
						exprs(x)[,chanCols] <- exp(log(exprs(x)[,chanCols])-(log(exprs(x)[,fsc]) %*% coeff.mat)) 
						diffMFI <- matrix(rep(initMFIs - log(apply(exprs(x)[,chanCols],2,median)),nrow(exprs(x))),ncol=length(chanCols),byrow=TRUE)
						exprs(x)[,chanCols] <- exp(log(exprs(x)[,chanCols]) + diffMFI)
						x@exprs[x@exprs<0] <- 0	
						x
					})
			
			return(fp)
		})




######################################################################################################################
##
## Function: read.flowPlate
##    
## Description: Read in plate fcs files.  The last three characters of the filename are usually the well id.
##  The annotated date frame describing the experiment and the sampleNames are then merged by well id.
##
######################################################################################################################

read.flowPlate = function(path=".",pattern=NULL,wellAnnotation) {
	
	files = list.files(path) 
	
	## Now get the well ids from the sample filenames
	if(is.null(pattern)) {
		Well.Id <- unlist(lapply(files, function(x) substr(x, start=nchar(x)-2,stop=nchar(x))))
	} else {
		Well.Id <- unlist(lapply(files, function(x) substr(x, start=nchar(x)-pattern[1],stop=nchar(x)-pattern[2])))
	}
	
	## Create a data frame with the sample file names and corresponding well id
	name = as.character(list.files(path))
	tempPlate.df <- data.frame(cbind(name,Well.Id),stringsAsFactors=FALSE)
	
	## Merge the tempPlate.df with a subset of the current phenoData (excluding name column from phenoData)
	## name gets overwritten with the sample file names
	tempPlate2.df <- merge(tempPlate.df,subset(pData(wellAnnotation),select=-name),by="Well.Id")
	
	## The name column is usually assumed to be in the first position for the other flowCore stuff
	name = as.character(tempPlate2.df$name)
	tempPlate2.df <- subset(tempPlate2.df,select=-name)
	tempPlate2.df <- data.frame(cbind(name,tempPlate2.df),stringsAsFactors=FALSE)
	
	## Now replace the phenoData in the wellAnnotation.adf with the new data frame
	row.names(tempPlate2.df) <- name
	pData(wellAnnotation) <- tempPlate2.df
	
	## It's still converting strings to factors for some reason, fixing it here for now.
	pData(wellAnnotation)$name <- as.character(pData(wellAnnotation)$name)
	
	## Call the read.flowSet function to actually read in the plates
	newFlowSet <- read.flowSet(path=path)
	
	## Copy over the annotated data frame describing the experiment
	phenoData(newFlowSet) <- wellAnnotation
	
	## Get rid of the dashes
	newFlowSet <- fsApply(newFlowSet,function(x) {
				newNames <- gsub("-",".",colnames(exprs(x)))
				colnames(exprs(x)) <- newNames
				x
			})
	
	gc()
	
	return(newFlowSet)
	
}

######################################################################################################################
##
## Function: getIsoGates
##    
## Description: For each isotype get the threshold for positive/negative.  unstainDelta
##	is the constant factor that is added to each gate after the 99th percentile is found, seems to generally
##	work well for our log10 transformed data.  
##
######################################################################################################################

getIsoGates <- function(data,cutOff=0.99,unstainDelta=0.024) {
	
	## Get the names of the columns with dyes
	dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	
	## Make a data frame containing the dye/well info for the isotype controls
	iso.df <- subset(pData(phenoData(data)),as.logical(Sample.Type=="Isotype"),select=c("Well.Number",dyeCols))
	
	## Get the list of dyes, wells, and column names corresponding to the dyes.  Note that this assumes one dye 
	## per well.
	isoDyes <- apply(iso.df,1,function(x) x[dyeCols[x[dyeCols]!=""]] )
	isoWells <- names(isoDyes)
	isoCols <- apply(iso.df,1,function(x) gsub("\\.","-",sub("\\.dye","",dyeCols[x[dyeCols]!=""]) ))
	
	## Create a hash for storing the rectangle gates.
	isoHash <- new.env(hash=TRUE,parent=emptyenv())
	
	## For each Isotype control, create a gate based on the channel of interest and foward scatter
	lapply(1:length(isoCols),function(x) {
				thresh <- quantile(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]]),probs=c(cutOff)) + unstainDelta 
				rangeF <- ((range(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]])))[1])	
				if(rangeF>0) {rangeF <- 0}
				names(thresh) <- isoCols[x]			
				names(rangeF) <- isoCols[x]
				assign(isoDyes[[x]],new("rectangleGate",filterId="rectangleGate",parameters=isoCols[x],min=rangeF,max=thresh), env=isoHash)
			})
	
	return(isoHash)
}

######################################################################################################################
##
## Function: getLinIsoGates
##    
## Description: For each isotype get the threshold for positive/negative.  It seems to work better if I first get the isotype
##	(negative control) gate values using the linear scale, MAD doesn't do a good job with log10 transformed data.
##	numMads is the number of mads above the median used to set the gate.   
##
######################################################################################################################

getLinIsoGates <- function(data,numMads=5,Samp.Type="Isotype") {
	
	
	## Get the names of the columns with dyes
	dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	
	## Make a data frame containing the dye/well info for the isotype controls
	iso.df <- subset(pData(phenoData(data)),as.logical(Sample.Type==Samp.Type),select=c("Well.Number",dyeCols))
	
	## Get the list of dyes, wells, and column names corresponding to the dyes.  Note that this assumes one dye 
	## per well.
	isoDyes <- apply(iso.df,1,function(x) x[dyeCols[x[dyeCols]!=""]] )
	isoWells <- names(isoDyes)
	isoCols <- apply(iso.df,1,function(x) sub("\\.dye","",dyeCols[x[dyeCols]!=""]))
	
	## Create a hash for storing the rectangle gates.
	isoHash <- new.env(hash=TRUE,parent=emptyenv())
	
	## For each Isotype control, create a gate based on the channel of interest and foward scatter
	lapply(1:length(isoCols),function(x) {
				mfi <- median((data[[isoWells[x]]]@exprs[,isoCols[x]]))
				mfi.mad <- mad((data[[isoWells[x]]]@exprs[,isoCols[x]]))
				isoMad <- numMads
				while(quantile(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]]),probs=0.975)>(mfi+isoMad*mfi.mad)) {
					isoMad <- isoMad + 0.1
				}
				thresh <- c(1,mfi+isoMad*mfi.mad) 
				rangeF <- (range(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]])))[1]	
				if(rangeF>0) {rangeF <- 0}
				rangeF <- c(0,rangeF)
				names(thresh) <- c("FSC.A",isoCols[x])			
				names(rangeF) <- c("FSC.A",isoCols[x])
				assign(isoDyes[[x]],new("rectangleGate",filterId="rectangleGate",parameters=c("FSC.A",isoCols[x]),min=rangeF,max=thresh), env=isoHash)
			})
	
	return(isoHash)
	
}

######################################################################################################################
##
## Function: getControlGates
##    
## Description: For each isotype get the threshold for positive/negative.  It seems to work better if I first get the isotype
##	(negative control) gate values using the linear scale, MAD doesn't do a good job with log10 transformed data.
##	numMads is the number of mads above the median used to set the gate.   
##
######################################################################################################################

getControlGates <- function(data,numMads=5,controlWells,channel) {
	
	## For each Isotype control, create a gate based on the channel of interest and foward scatter
	cutOff <- mean(sapply(controlWells,function(x) {
						mfi <- median((data[[x]]@exprs[,channel]))
						mfi.mad <- mad((data[[x]]@exprs[,channel]))
						thresh <-mfi+numMads*mfi.mad
						names(thresh) <- channel		
						thresh
					}))
	
}

######################################################################################################################
##
## Function: applyIsoGates
##    
## Description: Apply the isotype gates (actually the not of the gates) to identify the positive cells.  The *.dye
##	column in phenoData indicated which dye is in the channel, and the isotypes are specific to a dye/channel.
##
######################################################################################################################

applyIsoGates <- function(data,isoGates) {
	
	## Get the names of the columns with dyes
	dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	
	## Filter the cells 
	apply(pData(phenoData(data)),1,function(x) {
				lapply(dyeCols,function(y) {
							if((x[y]!="") & (x[y] %in% ls(isoGates))) {
								tempFilt <- get(x[y],env=isoGates)
								filter(data[[x["name"]]], !tempFilt)
							} 
							else { NA }
						})
			})
	
}

######################################################################################################################
##
## Function: checkIsoGates
##    
## Description: Make sure that the isotype gates are not much higher than the negative test wells, otherwise use the
##	unstained sample to gate.
##
######################################################################################################################

checkIsoGates <- function(data,plateSumm.df,isoGates,unstain,numMads=5) {
	
	## Identify the wells containing unstained samples
	if(missing(unstain)) {
		unstain <- unlist(subset(pData(phenoData(data)),as.logical((Sample.Type=="Unstained") ),select="name"))[[1]]
	}
	
	## Get the names of the columns with dyes
	dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	isoCols <- sapply(dyeCols,function(x) sub("\\.dye","",x))
	
	unstainGates <- lapply(1:length(isoCols),function(x) {
				mfi <- median((data[[unstain]]@exprs[,isoCols[x]]))
				mfi.mad <- mad((data[[unstain]]@exprs[,isoCols[x]]))
				isoMad <- numMads
				while(quantile(unlist(data[[unstain]]@exprs[,isoCols[x]]),probs=0.99)>(mfi+isoMad*mfi.mad)) {
					isoMad <- isoMad + 0.1
				}
				thresh <- c(1,mfi+isoMad*mfi.mad) 
				rangeF <- (range(unlist(data[[unstain]]@exprs[,isoCols[x]])))[1]	
				if(rangeF>0) {rangeF <- 0}
				rangeF <- c(0,rangeF)
				names(thresh) <- c("FSC.A",isoCols[x])			
				names(rangeF) <- c("FSC.A",isoCols[x])
				new("rectangleGate",filterId=paste(isoCols[x]),parameters=c("FSC.A",isoCols[x]),min=rangeF,max=thresh)
			})
	names(unstainGates) <- isoCols
	
	## Set up the gates by dye, since that's how the isotype gates are contstructed
	unstainGates <- sapply(ls(isoGates),function(x) {
				tempFilt <- get(x,env=isoGates)
				chan <- paste(parameters(tempFilt)[[2]])
				unstainGates[[chan]]
			})
	
	## Get the test wells for each isotype control that are <=1% positive
	testWells <- sapply(ls(isoGates),function(x) {
				tempFilt <- get(x,env=isoGates)		
				sampNames <- sampleNames(data)[pData(phenoData(data))[,paste(parameters(tempFilt)[[2]],".dye",sep="")]==x]	
				sampNames <- sampNames[pData(phenoData(data[sampNames]))[,"Sample.Type"]=="Test"]
				chanPP <- paste(parameters(tempFilt)[[2]],".pp",sep="")
				sampNames[plateSumm.df[plateSumm.df$name %in% sampNames,chanPP]<=1]
			})
	
	lapply(ls(isoGates), function(x) {
				if(length(testWells[[x]])>0) {
					unstainPP <- sapply(testWells[[x]],function(y) {	
								unPP <- 100*sum(filter(data[[y]],!unstainGates[[x]])@subSet)/(nrow(exprs(data[[y]]))+1)
								isoPP <- 100*sum(filter(data[[y]],!isoGates[[x]])@subSet)/(nrow(exprs(data[[y]]))+1)
								if(unPP <= 2 & unPP <= isoPP) FALSE
								else TRUE
							})
					if(sum(unstainPP)==0) assign(x,unstainGates[[x]],env=isoGates) 
				}
			})
	
	return(isoGates)	
}

checkIsoGates2 <- function(data,plateSumm.df,isoGates) {
	
	## Get the test wells for each isotype control that are <=1% positive
	testWells <- sapply(ls(isoGates),function(x) {
				tempFilt <- get(x,env=isoGates)		
				sampNames <- sampleNames(data)[pData(phenoData(data))[,paste(parameters(tempFilt)[[2]],".dye",sep="")]==x]	
				sampNames <- sampNames[pData(phenoData(data[sampNames]))[,"Sample.Type"]=="Test"]
				chanPP <- paste(parameters(tempFilt)[[2]],".pp",sep="")
				sampNames[plateSumm.df[plateSumm.df$name %in% sampNames,chanPP]<=1]
			})
	
	lapply(ls(isoGates), function(x) {
				if(length(testWells[[x]])>0) {
					tempGate <- get(x,env=isoGates)
					delta <- tempGate@max[2]/200
					while(delta>0) {
						pp <- sapply(testWells[[x]],function(y) {
									100*sum(filter(data[[y]],!tempGate)@subSet)/nrow(exprs(data[[y]]))
								})
						if(mean(pp)>=2) {
							assign(x,tempGate,env=isoGates)
							delta <- -1
						} else { 
							tempGate@max[2] <- tempGate@max[2] - delta
						}
						
					}
				}
			})
	
	return(isoGates)	
}

######################################################################################################################
##
## Function: getMarkers
##    
## Description: Take in an annotated data frame describing a plate and the plate MFI values, return a list of genes (markers) and their
## associated MFI values.
##
######################################################################################################################

getMarkers <- function(data,summ.df,chanCols) {
	
	## Markers (genes) are only in the "Test" samples
	genes.df <- subset(pData(phenoData(data))[,c("Well.Id",chanCols,"Sample.Type")],Sample.Type=="Test",select=-Sample.Type)
	
	## Make a list of data frames showing the well/channel for each marker 	
	genes2.df <- lapply(chanCols,function(x) {
				temp.df <- cbind(genes.df[,c(x,"Well.Id")],x)
				colnames(temp.df) <- c("Antigen","Well.Id","Channel")
				temp.df 
			})
	
	## Now bind the list of dataframes into a single data frame
	temp.df <- data.frame(matrix(nrow=0,ncol=3))
	colnames(temp.df) <- c("Antigen","Well.Id","Channel")
	for(i in 1:length(genes2.df)) {
		temp.df <- rbind(temp.df,genes2.df[[i]])
	}
	geneMap.df <- temp.df[temp.df[,"Antigen"]!="",]
	
	rownames(summ.df) <- as.character(summ.df$Well.Id)
	
	## Now get the marker MFI, etc. values from the plate data
	out.mat <- apply(geneMap.df,1,function(x) {
				matrix( c(x["Antigen"],x[["Well.Id"]],summ.df[x[["Well.Id"]],c(grep(x[["Channel"]],colnames(summ.df)))]) ,nrow=1 )
			})
	
	cols <- c("Antigen","Well.Id",sub(paste(chanCols[1],"\\.",sep=""),"",colnames(summ.df)[grep(chanCols[1],colnames(summ.df))])) 
	
	out.mat <- matrix(unlist(out.mat,recursive=TRUE),ncol=length(cols),byrow=TRUE)
	colnames(out.mat) <- cols
	out.df <- data.frame(out.mat,stringsAsFactors=FALSE)
	fjPP <- out.df$pp
	ppCol <- which(cols=="pp")	
	cbind(out.df[,1:ppCol],fjPP,out.df[,(ppCol+1):length(cols)])
}

######################################################################################################################
##
## Function: calcIsoRatio
##    
## Description:  Calculate the ratios of MFI of test wells vs isotypes and unstained for total, positive, and negative mfis
##
######################################################################################################################

calcIsoRatio <- function(data,summ.df,channels) {
	
	mfiCols <- unlist(lapply(channels,function(x) paste(x,".ratioMFI",sep="")))
	posCols <- unlist(lapply(channels,function(x) paste(x,".MFIratio.pos",sep="")))
	negCols <- unlist(lapply(channels,function(x) paste(x,".MFIratio.neg",sep="")))
	
	temp.df <- data.frame(matrix(0,nrow=nrow(summ.df),ncol=length(mfiCols)*3))
	colnames(temp.df) <- c(mfiCols,posCols,negCols)
	
	## Get the names of the columns with dyes
	dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]	
	
	
	## Make a data frame containing the dye/well info for the isotype controls
	iso.df <- subset(pData(phenoData(data)),as.logical(Sample.Type=="Isotype"),select=c("Well.Id",dyeCols))
	
	## Get the list of dyes, wells, and column names corresponding to the dyes.  Note that this assumes one dye 
	## per well.
	isoDyes <- apply(iso.df,1,function(x) x[dyeCols[x[dyeCols]!=""]] )
	isoWells <- as.character(iso.df[,1])
	isoCols <- apply(iso.df,1,function(x) sub("dye","MFI",dyeCols[x[dyeCols]!=""]) )
	
	## Create a hash for storing the rectangle gates.
	isoHash <- new.env(hash=TRUE,parent=emptyenv())
	
	lapply(1:length(isoCols),function(x) {
				assign(isoDyes[[x]],summ.df[isoWells[x],isoCols[[x]]],env=isoHash)
			})
	
	
	iso.list <- apply(pData(phenoData(data))[,dyeCols],1, function(x) {
				lapply(x,function(y) {
							if(y!="") { 
								get(y,env=isoHash) 
							}
							else{ NA } 
						})
			})
	
	iso.df <- data.frame(matrix(unlist(iso.list),ncol=length(channels),byrow=TRUE))	
	
	temp.df[,mfiCols] <- round(summ.df[,unlist(lapply(channels,function(x) gsub("-","\\.",paste(x,".MFI",sep=""))))]/iso.df,digits=3)
	temp.df[,posCols] <- round(summ.df[,grep(".*posMFI$",colnames(summ.df))]/iso.df,digits=3)
	temp.df[,negCols] <- round(summ.df[,grep(".*negMFI$",colnames(summ.df))]/iso.df,digits=3)
	
	names(temp.df) <- sub("-","\\.",names(temp.df))
	
	summ.df <- cbind(summ.df,temp.df)
	
	
	## Also calculate the unstained ratios, for now we'll just average over unstained wells
	unstainWells <- unlist(subset(pData(phenoData(data)),as.logical((Sample.Type=="Unstained") ),select="name"))
	
	unstainMFIs <- fsApply(plateSet[unstainWells], function(x) {
				apply(exprs(x)[,channels],2,median)
			})	
	
	unstainMFIs.mat <- matrix(apply(unstainMFIs,2,mean),nrow=1)
	
	mfiCols2 <- unlist(lapply(channels,function(x) paste(x,".unstainRatioMFI",sep="")))
	
	temp2.df <- data.frame(matrix(0,nrow=nrow(summ.df),ncol=length(mfiCols2)))
	
	temp2.df <- summ.df[,unlist(lapply(channels,function(x){paste(x,".MFI",sep="")}))]/unstainMFIs.mat
	colnames(temp2.df) <- c(mfiCols2)
	
	temp2.df <- round(temp2.df,digits=3)
	
	names(temp2.df) <- sub("-","\\.",names(temp2.df))
	cbind(summ.df,temp2.df)
	
	
	
}