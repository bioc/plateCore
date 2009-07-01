################################################################################
##
## Filename: flowPlate-accessors.R
##
## Author: straine
##
## Data: Oct 10, 2007
##
## Description: Methods for accessing flowPlate objects. 
##
################################################################################

################################################################################
## This method creates a virtual plate by combining 2 or more flowPlates.
## If the sample names are not unique, then they are changed using the 
## make.unique function. I need to support binding a list of flowPlates, since
## the currenty approach of listing them in the arguments is not elegant. Also,
## it would be nice to be able to do an as(flowPlate,list), similar to 
## as(flowSet,list) in flowCore.
## #############################################################################
setMethod("fpbind",signature("flowPlate","flowPlate"),
		function(p1,p2,...,plateName="VirtPlate") {
			
	## Get the number of arguments
	na <- nargs()
	argl <- list(p1,p2,...)
			
	## Get rid of anything that's not a flowPlate'
	while(na > 0 && is.null(argl[[na]])) { argl <- argl[-na]; na <- na - 1 }
	argl <- argl[sapply(argl,class)=="flowPlate"]
			
	## For now, don't bind two plates with the same name
	plateIds <- unlist(sapply(argl,function(x) unique(x@wellAnnotation$plateName))	)
	if(length(plateIds) != length(unique(plateIds))) {
		stop("Binding flowPlates with identical plateNames is not supported.")
	}
			
	## This method can only bind flowPlates where the well annotation data frames have the same column headers
	## I think I'll change this later, since well annotation changes over
	## the course of an analysis (results columns are added)
	annl <- lapply(argl,function(x) x@wellAnnotation)	
	for(i in 2:length(argl)) {
		if(!all.equal(colnames(annl[[1]]),colnames(annl[[i]]))) {
			stop("flowPlate Well Annotations must have identical column names.")
		}
	}
			
	## The sample names have to be unique in order for the phenoData on
	## the flowSet to work correctly. Can't have rows with identical names
	sampNames <- unlist(lapply(argl,sampleNames))		
	uniqNames <-make.unique(sampNames)
			
	sampLength <- lapply(argl,function(x) length(sampleNames(x)))
	
	## Vector of lists containing sample names for each plate to be merged
	uniqList <- vector("list",length=length(sampLength))		
	uniqList <- lapply(1:length(sampLength), function(i) {
		index2 <- sum(sapply(sampLength[1:i],function(x) x[[1]]))
		index1 <- index2 - sampLength[[i]] + 1
		temp <- uniqNames[index1:index2]
		names(temp) <- sampleNames(argl[[i]])
		temp
	})
			
	## Remake the annotation list with the unique names
	annl <- lapply(1:length(argl),function(i) {			
		anndf <- argl[[i]]@wellAnnotation
		anndf$name <- uniqList[[i]][unlist(anndf$name)]
		anndf$name[1]
		anndf <- data.frame(name=anndf$name,subset(anndf,select=-name),stringsAsFactors=FALSE)
	})
			
	## rbind the annotations into a single data frame
	wellAnnotation <- annl[[1]]
	for(i in 2:length(annl)) {
		wellAnnotation <- rbind(wellAnnotation,annl[[i]])
	}
			
	## Get the frames		
	frames <- unlist(lapply(1:length(argl),function(i) {			
		x <- argl[[i]]
		temp <- x@plateSet@frames
		fsList <- lapply(names(uniqList[[i]]),function(fileName) {
			temp[[fileName]]							
		})
		names(fsList) <- uniqList[[i]]
		fsList
	}))
				
	## Now make a flowPlate
	np <- flowPlate(as(frames,"flowSet"),wellAnnotation,plateName=plateName)
			
	return(np)
})

################################################################################
## Get groups, which are currently just the negative control groups. 
## These consist of the negative control well and all the associated test wells
################################################################################
setMethod("getGroups",signature("flowPlate"),function(data,type="Negative.Control",chan,...) {
			
	## Declaring Variables, otherwise R check throws an error
	Channel <- ""
	Negative.Control <- ""
	name <- ""
	Well.Id <- ""
			
	## Get the list negative control wells
	wellIds <- unique(data@wellAnnotation$Negative.Control)
	
	## If the negative control well is not in this flowPlate, then don't
	## include it.
	wellIds <- wellIds[wellIds %in% pData(phenoData(data@plateSet))$Well.Id]
			
	wells <- list()
			
	if(type=="Negative.Control") {
		wells <- lapply(wellIds,function(x) {
			wells <- unlist(subset(data@wellAnnotation,Channel==chan & Negative.Control==x,select=name))
			wells <- c(unlist(subset(data@wellAnnotation,Channel==chan & Well.Id==x,select=name)),wells)
			wells <- unique(wells)	
			if(!length(wells)) {
				return(NA)
			} else {
				return(wells)
			}
		})
	} else {
		stop("invalid option for type")
	}
			
	wells <- wells[!is.na(wells)]
	return(wells)
})


################################################################################
## This is used to transform the flowSet and  the isotype gates.  Other types
## of gates should be included later. 
################################################################################
setMethod("%on%",signature(e2="flowPlate"),function(e1,e2) {
	
	## If the data is transformed, then any gates associated with the plate
	## should also be transformed if possible.
	if("Negative.Control.Gate" %in% colnames(e2@wellAnnotation)) {
		for(y in e1@transforms) {
			wellList <- which(e2@wellAnnotation$Channel %in% y@input)	
			e2@wellAnnotation[wellList,"Negative.Control.Gate"] <- y@f(e2@wellAnnotation[wellList,"Negative.Control.Gate"])			  	
		}
	}
	
	## Apply the transformation to the flowSet
	e2@plateSet <- fsApply(e2@plateSet,"%on%",e1=e1)	
	e2
			
})

################################################################################
## Constructor for flowPlates
################################################################################
setMethod("flowPlate",signature("flowSet"),function(data,wellAnnotation,plateName="",...) {
			
	## Getting rid of R CMD Check errors
	name <- ""
	
	## Add plateName to well annotation if it doesn't exist
	if(!"plateName" %in% colnames(wellAnnotation)) {
		wellAnnot <- data.frame(wellAnnotation,plateName=plateName,stringsAsFactors=FALSE)
	} else {
		wellAnnot <- wellAnnotation
	}
	
	## Create a data.frame that corresponds to phenoData.
	## All info about a well is on a single row
	config <- makePlateLayout(wellAnnot)
			
	##Add the phenoData to the flowSet
	data <- flowPhenoMerge(data,config)
			
	## When a flowPlate is first created from a plateLayout template, names are 
	## missing and need to be added.  Names correspond to filenames
	if(colnames(wellAnnot)[1]=="Well.Id") {
		wellAnnot$name <- sapply(wellAnnot$Well.Id,function(x) {
			rownames(pData(phenoData(data))[pData(phenoData(data))$Well.Id==x,])
		}) 
	}
			
	## Only keep annotation for samples that are in  the dataset
	wellAnnot <- subset(wellAnnot,name %in% sampleNames(data) )
	
	## No create the flowPlate
	temp <- new("flowPlate")
	temp@plateName <- plateName
	temp@wellAnnotation <- wellAnnot
	temp@plateSet <- data
			
	return(temp)
})

######################################################################
## Fit linear model to FSC vs channel of interest.
## Correct for autofluoresence and then set the median back to it's orignal
## value.  The fitting is performed on log  transformed data.
######################################################################
setMethod("fixAutoFl",signature("flowPlate"), 
		function(fp,fsc,chanCols,unstain=NULL,minCut=10) {
			
	plateSet <- fp@plateSet
			
	## Identify the wells containing unstained samples
	if(is.null(unstain)) {
		unstainWells <- unique(unlist(subset(fp@wellAnnotation,Sample.Type=="Unstained",select="name")))
	} else {
		unstainWells <- unstain
	}
	
	## One or more of the unstained wells needs to have at least 10 or more 
	## events
	numEvents <- fsApply(plateSet[unstainWells],function(x) {
		nrow(exprs(x))
	})
	if (all(numEvents<10)) {
		stop("Unstained wells contain too few events!")
	} else {
		unstainWells <- unstainWells[numEvents>=10]
	}
	
	## Fit a linear model to the unstained data, get the slope
	unstainFits <- fsApply(plateSet[unstainWells], function(x) {
		unlist(lapply(chanCols,function(y) {
			## pdh: the data less than 10 seem to throw off the fit
			## eas: maybe throw a warning?  Or ignore wells with <10 if we have
			## other unstained wells?
			exprData <- exprs(x)[(exprs(x)[,y]> minCut),]	
				ltsReg(as.formula(paste(gsub("-",".",y)," ~ ",gsub("-",".",fsc),sep="")), log(data.frame(exprData)))$coefficients[gsub("-",".",fsc)]
		}))
	})	
			
	## Get a matrix of mean slopes
	coeff.mat <- matrix(apply(unstainFits,2,mean),nrow=1)
	colnames(coeff.mat) <- chanCols
			
	## Apply the correction to the channels of interest in chanCols
	## Truncate values less than 0 at 0
	plateSet <- fsApply(plateSet,function(x) {
		initMFIs <- log(apply(exprs(x)[,chanCols],2,median))
		## these are values that are 1 (minimum value before the correction, we want them to be the same after the
		## corrrection
		selectOnes <- (exprs(x) == 1)
		exprs(x)[,chanCols] <- exp(log(exprs(x)[,chanCols])-(log(exprs(x)[,fsc]) %*% coeff.mat)) 
		diffMFI <- matrix(rep(initMFIs - log(apply(exprs(x)[,chanCols],2,median)),nrow(exprs(x))),ncol=length(chanCols),byrow=TRUE)
		exprs(x)[,chanCols] <- exp(log(exprs(x)[,chanCols]) + diffMFI)
		x@exprs[x@exprs<0] <- 0	
		## no reset any values that were 1 before to be 1 again
		x@exprs[selectOnes] <- 1
			x
	})
		
	fp@plateSet <- plateSet
	return(fp)
})

################################################################################
## Compensate, but only for dyes that were added. 
## I shouldn't need to make a copy of the flowSet, but I did in case it 
## would save me trouble with environments hanging around.
################################################################################
setMethod("compensate", signature(x="flowPlate"), function(x,spillover) {
			
	temp <- x@plateSet@frames
			
	frames <- lapply(ls(temp),function(fileName) {
		well <- pData(phenoData(x@plateSet))[fileName,]		
		dyeCols <- colnames(spillover)[!is.na(well[,gsub("-",".",colnames(spillover))])]
		
		## Use the compensate method from flowCore on individual wells
		if(length(dyeCols)>=2) {
			compensate(temp[[fileName]],spillover[dyeCols,dyeCols])			
		} else {
			temp[[fileName]]
		}								
	})
			
	names(frames) <- ls(temp)	
			
	config <- makePlateLayout(x@wellAnnotation)
			
	plateSet <- flowPhenoMerge(as(frames,"flowSet"),config)
			
	x@plateSet <- plateSet
			
	return(x)
})


################################################################################
## The summaryStats method calculates some simple statistics for a flowPlate
## that has been processed using setControlGates and applyControlGates. Results
## are stored as new columns in the wellAnnotation data.frame
################################################################################
setMethod("summaryStats", signature("flowPlate"), function(data,...) {
			
	## Declare variables for R check
	Well.Id <- ""
	Channel <- ""
	plateName <- ""
	Sample.Type <- ""

	## Make a local copy of the well annotation
	wellAnnotation <- data.frame(data@wellAnnotation)
	
	## Calculate the median fluoresence intensities
	wellAnnotation$MFI <- apply(wellAnnotation,1,function(x) {
		chan <- x[["Channel"]]
		frame <- data@plateSet[[x[["name"]]]]				
		if(chan %in% colnames(exprs(frame)) && nrow(exprs(frame))>0) {
			return(median(exprs(frame)[,chan]))
		} else { NA }	
	})
	
	## Calculate the MFI ratio, which is the ratio of MFIs of test wells to 
	## their associated isotype control on a linear scale.
	wellAnnotation$MFI.Ratio <- apply(wellAnnotation,1,function(x) {
		negWell <- x[["Negative.Control"]]
		chan <- x[["Channel"]]
		plateN <- x[["plateName"]]
		frame <- data@plateSet[[x[["name"]]]]
						
		if(chan %in% colnames(exprs(frame)) && negWell!="" && nrow(exprs(frame))>0) {
			negMFI <- subset(wellAnnotation,Well.Id==negWell & Channel==chan & plateName==plateN)$MFI
			if(!is.na(negMFI) && !is.na(x[["MFI"]])) {
				return(as.numeric(x[["MFI"]])/as.numeric(negMFI))
			}
		} else { 
			return(NA)
		}	
	})	
			
			
	## Fit a robust glm to the MFI ratio vs Percent Positive to assess the
	## quality of the gating, uses glmrob from robustbase
	if("Percent.Positive" %in% colnames(wellAnnotation)) {
		df <- subset(wellAnnotation, Sample.Type=="Test")
		df$MFI.Ratio <- log(df$MFI.Ratio)
		PosCount <- round(df$Percent.Positive,0)
		NegCount <- 100-PosCount			
		robMFI <- glmrob(cbind(PosCount,NegCount) ~ MFI.Ratio, data=df,family=binomial(link="logit"))		
		x <- -robMFI$coefficients[[1]]-robMFI$coefficients[[2]]*log(wellAnnotation$MFI.Ratio)
		wellAnnotation$Predict.PP <- 100/(1 + exp(x))		
		diffR <- wellAnnotation$Percent.Positive-wellAnnotation$Predict.PP	
			
		resids <- (robMFI$residuals-mean(robMFI$residuals,na.rm=TRUE))/sd(robMFI$residuals,na.rm=TRUE)
					
		wellAnnotation$Gate.Score <- NA
		wellAnnotation$Gate.Score[wellAnnotation$Sample.Type=="Test"] <- resids
	}
			
	data@wellAnnotation <- wellAnnotation
			
	return(data)
})

################################################################################
## setControl gates the negative control wells to establish a threshold between
## positive and negative cells. Currenty, this threshold can be estimated using
## either some number of MADs or a set quantile.
################################################################################
setMethod("setControlGates", signature("flowPlate"), function(data,gateType="Negative.Control",
				threshType="MAD",numMads=5,isoquantile=.995,minCut=NULL,...) {

	if(gateType=="Negative.Control") {
		## First get the control gate for each of the isotype groups.		
		isoWells <- subset(data@wellAnnotation,(Sample.Type %in% c("Isotype","Negative.Control")) & !is.na(Channel))
			
		## Create a threshold based on median absolute deviations
		if(threshType=="MAD"){
			isoGates <- lapply(unique(isoWells$name), function(x) {
				sapply(isoWells[isoWells$name==x,"Channel"], function(i) {	
					mfi <- median(exprs(data[[x]])[,i])
					mfi.mad <- mad(exprs(data[[x]])[,i])
					thresh <- mfi+numMads*mfi.mad	
					thresh
				})		
			})
		## Quantile based threshold
		} else if(threshType=="isoQuant"){
			isoGates <- lapply(unique(isoWells$name), function(x) {
				sapply(isoWells[isoWells$name==x,"Channel"], function(i) {	
					y <- exprs(data[[x]])[,i]
					if(is.null(minCut))	{
						thresh <- quantile(y,isoquantile)
					} else {
						y = y[y>=minCut]
						thresh <- quantile(y,isoquantile)
					}
					thresh <- as.numeric(thresh)								
					})		
			})
		} else {
			stop("invalid option for threshType")
		}
	
		## Copy the gates to the corresponding test wells
		names(isoGates) <- unique(isoWells$name)
		data@wellAnnotation$Negative.Control.Gate <- apply(data@wellAnnotation,1,function(x) {
							
			well <- subset(data@wellAnnotation,Well.Id== x[["Negative.Control"]] & plateName == x[["plateName"]],select=name)[1,][[1]]				
			
			if(length(well) && well %in% names(isoGates)) {
				isoGates[[well]][x[["Channel"]]]
			} else if (x[["name"]] %in% unique(isoWells$name)) {
				isoGates[[x[["name"]]]][x[["Channel"]]]
			} else {
				return(NA)
			}
		}) 
	} else {
		stop("gateType not supported")
	}		
	
	return(data)
})

################################################################################
## Gates/Thresholds created using setControlGates are applied to the test wells
## and the percentage of positive cells and number of positive cells are 
## reported.
################################################################################
setMethod("applyControlGates", signature("flowPlate"), function(data,gateType="Negative.Control",...) {
			
	if(gateType %in% c("Isogate","Negative.Control")) {
		## First get the control gate for each of the isotype groups.	
		wa <- data@wellAnnotation 
				
		##If any applyControlGate or summaryStats columns already exist, get rid of them
		newNames <- c('Percent.Positive','Total.Count','Positive.Count','MFI','MFI.Ratio','Predict.PP','Gate.Score')
		if(any(newNames %in% colnames(wa))) {
			keepCols <- which(!colnames(wa) %in% newNames)
			wa <- wa[,keepCols]
		}	
		
		## Calculate percent positive
		wa$Percent.Positive <- apply(data@wellAnnotation,1,function(x) {
		thresh <- as.numeric(x[["Negative.Control.Gate"]])	
		chan <- x[["Channel"]]
		frame <- data@plateSet[[x[["name"]]]]
							
		if(!is.na(thresh) && chan %in% colnames(exprs(frame)) && nrow(exprs(frame))>0 ) {
			iso <- new("rectangleGate",filterId="rectangleGate",parameters=chan,min=thresh,max=Inf)
			isoResult <- filter(frame,iso)
			return(100*(sum(isoResult@subSet)/nrow(exprs(frame))))
		} else { 
			return(NA)
		}	
						})	
 		## Get the total number of events in each well	
		wa$Total.Count <- apply(data@wellAnnotation,1,function(x) {
			thresh <- as.numeric(x[["Negative.Control.Gate"]])	
			chan <- x[["Channel"]]
			frame <- data@plateSet[[x[["name"]]]]
							
			if(!is.na(thresh) && chan %in% colnames(exprs(frame)) && nrow(exprs(frame))>0 ) {
				return(nrow(exprs(frame)))
			} else { 
				return(NA)
			}	
		})	
		
		## Calculate the number of positive events
		wa$Positive.Count <- apply(data@wellAnnotation,1,function(x) {
			thresh <- as.numeric(x[["Negative.Control.Gate"]])	
			chan <- x[["Channel"]]
			frame <- data@plateSet[[x[["name"]]]]
							
			if(!is.na(thresh) && chan %in% colnames(exprs(frame)) && nrow(exprs(frame))>0 ) {
				iso <- new("rectangleGate",filterId="rectangleGate",parameters=chan,min=thresh,max=Inf)
				isoResult <- filter(frame,iso)
				return(sum(isoResult@subSet))
			} else { 
				return(NA)
			}	
		})	
				
				
		data@wellAnnotation <- wa
		
	} else {
		stop("gateType not supported")
	}
	
	return(data)
})

################################################################################
## Some convenience functions, primarily wrappers to flowCore, Biobase
################################################################################

setMethod("phenoData","flowPlate",function(object) {
			return(phenoData(object@plateSet))
		})

setMethod("plateSet","flowPlate",function(fp,...) {
			return(fp@plateSet)
		})

setMethod("Subset","flowPlate",function(x,subset,select=NULL,...) {

			if(is.null(select)) {
				x@plateSet <- Subset(x@plateSet,subset)
			} else {
				select <- unlist(select)
				copy = select[select %in% sampleNames(x)]
				if(length(copy) != length(select)) copy <- c(copy,sampleNames(x)[pData(phenoData(x@plateSet))[,"Well.Id"] %in% select[!(select %in% sampleNames(x))]])
				
				temp <- x@plateSet@frames
				
				frames <- lapply(ls(temp),function(fileName) {
							well <- pData(phenoData(x@plateSet))[fileName,]		
							
							if(fileName %in% copy) {
								Subset(temp[[fileName]],subset)		
							} else {
								temp[[fileName]]
							}								
						})
				
				names(frames) <- ls(temp)	
				
				config <- makePlateLayout(x@wellAnnotation)
				
				plateSet <- flowPhenoMerge(as(frames,"flowSet"),config)
				
				x@plateSet <- plateSet
			}	
			
			return(x)
			
		})

setMethod("sampleNames","flowPlate",function(object) {
			sampleNames(phenoData(object@plateSet))
		})

setMethod("wellAnnotation","flowPlate",function(fp,...) {
			return(fp@wellAnnotation)
		})

################################################################################
## subsetting methods
################################################################################
setMethod("[",c("flowPlate"),function(x,i,j,...,drop=FALSE) {
			
			name <- ""
			
			if(missing(drop)) drop = FALSE
			if(missing(i)) return(x)
			
			if(is.numeric(i) || is.logical(i)) {
				copy = sampleNames(x)[i]
			} else {
				i <- unlist(i)
				copy = i[i %in% sampleNames(x)]
				if(length(copy) != length(i)) copy <- c(copy,sampleNames(x)[pData(phenoData(x@plateSet))[,"Well.Id"] %in% i[!(i %in% sampleNames(x))]])
			}
			
			x@plateSet <- x@plateSet[copy]
			x@wellAnnotation <- subset(x@wellAnnotation,name %in% copy)
			x
		})


setMethod("[[","flowPlate",function(x,i,j,...) {
			if(length(i)!=1)
				stop("subscript out of bounds (index must have length 1)")
			if(is.numeric(i)) fr <- sampleNames(x)[[i]]
			else if(i %in% sampleNames(x)) fr <- i
			else fr <- sampleNames(x)[pData(phenoData(x@plateSet))[,"Well.Id"] == i]
			
			if(is.na(fr)) stop("subscript out of bounds")
			fr <- x@plateSet[[fr]]
			
		})	
