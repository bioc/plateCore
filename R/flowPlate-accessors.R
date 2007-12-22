#########################################################################################################
##
## Filename: flowPlate-accessors.R
##
## Author: straine
##
## Data: Oct 10, 2007
##
## Description:
##
#########################################################################################################


setMethod("%on%",signature(e2="flowPlate"),function(e1,e2) {
		e2@plateSet <- fsApply(e2@plateSet,"%on%",e1=e1)	
		e2
	})

setMethod("flowPlate",signature("flowSet"),function(data,wellAnnot,...) {
			
			config <- makePlateLayout(wellAnnot)
			
			temp <- new("flowPlate")
					
			wellAnnot$Channel <- gsub("-",".",wellAnnot$Channel)
			temp@wellAnnotation <- wellAnnot
			
			## Get rid of dashes in flowSet because they're annoying for lattice
			data <- fsApply(data,function(x) {
						newNames <- gsub("-",".",colnames(exprs(x)))
						colnames(exprs(x)) <- newNames
						x
					})

			data <- flowPhenoMerge(data,config)
			
			temp@plateSet <- data
			return(temp)
		})

setMethod("setRange",signature("flowPlate","numeric","numeric","character"), function(x,minF,maxF,type="truncate") {
			if(type=="truncate") {
				x@plateSet <- fsApply(x@plateSet,function(z) {
							z@exprs[z@exprs>=maxF] <- maxF
							z@exprs[z@exprs<=minF] <- minF	
							z
						})
				
			} else if (type=="remove"){
				x@plateSet <- fsApply(x@plateSet,function(z) {
							z@exprs <- z@exprs[apply(exprs(z),1, function(y) {sum(!y<=maxF)==0}),]
							z@exprs <- z@exprs[apply(exprs(z),1, function(y) {sum(!y>=minF)==0}),]	
							z
						})
				
			}
			return(x)		
		})

setMethod("fixAutoFl",signature("flowPlate"), 
		function(fp,fsc,chanCols,unstain=NULL) {
			
			plateSet <- fp@plateSet
			## Identify the wells containing unstained samples
			if(is.null(unstain)) {
				unstainWells <- unlist(subset(pData(phenoData(plateSet)),as.logical((Sample.Type=="Unstained") ),select="name"))
			} else {
				unstainWells <- unstain
			}
			
			## Fit a linear model to the unstained data, get the slope
			unstainFits <- fsApply(plateSet[unstainWells], function(x) {
						unlist(lapply(chanCols,function(y) {
											exprData <- exprs(x)[(exprs(x)[,y]> 1),]	
											ltsReg(as.formula(paste(y," ~ ",fsc,sep="")), log(data.frame(exprData)))$coefficients[fsc]
										}))
					})	
			
			## Get a matrix of mean slopes
			coeff.mat <- matrix(apply(unstainFits,2,mean),nrow=1)
			
			## Apply the correction to the channels of interest in chanCols\
			## Truncate values less than 0 at 0
			plateSet <- fsApply(plateSet,function(x) {
						initMFIs <- log(apply(exprs(x)[,chanCols],2,median))	
						exprs(x)[,chanCols] <- exp(log(exprs(x)[,chanCols])-(log(exprs(x)[,fsc]) %*% coeff.mat)) 
						diffMFI <- matrix(rep(initMFIs - log(apply(exprs(x)[,chanCols],2,median)),nrow(exprs(x))),ncol=length(chanCols),byrow=TRUE)
						exprs(x)[,chanCols] <- exp(log(exprs(x)[,chanCols]) + diffMFI)
						x@exprs[x@exprs<0] <- 0	
						x
					})
			fp@plateSet <- plateSet
			return(fp)
		})

setMethod("compensate", signature(x="flowPlate"), function(x,spillover) {

	plateSet <- fsApply(x@plateSet,function(y) {
		fileName <- attributes(y)$descriptio[["$FIL"]]
		well <- pData(phenoData(x@plateSet))[fileName,]	
		dyeCols <- colnames(spillover)[!is.na(well[,colnames(spillover)])]
		
		if(length(dyeCols)>=2) {
			y <- compensate(y,spillover[dyeCols,dyeCols])
			
		} else {
			y
		}
	})	
	
	x@plateSet <- plateSet
	return(x)
})
	

setMethod("setContolGates", signature("flowPlate"), function(data,gateType="Negative.Control",numMads=5,...) {
			
	if(gateType=="Negative.Control") {
		## First get the control gate for each of the isotype groups.	

		isoWells <- subset(data@wellAnnotation,Sample.Type=="Isotype" & !is.na(Channel))
		
		isoGates <- lapply(unique(isoWells$Well.Id), function(x) {
			sapply(isoWells[isoWells$Well.Id==x,"Channel"], function(i) {	
				mfi <- median(exprs(data[[x]])[,i])
				mfi.mad <- mad(exprs(data[[x]])[,i])
				isoMad <- numMads
				while(quantile(exprs(data[[x]])[,i],probs=0.975)>(mfi+isoMad*mfi.mad)) {
					isoMad <- isoMad + 0.1
				}
				thresh <- mfi+isoMad*mfi.mad	
				thresh
			})		
		 })
		names(isoGates) <- unique(isoWells$Well.Id)

		data@wellAnnotation$Isogate <- apply(data@wellAnnotation,1,function(x) {
			well <- x[["Negative.Control"]]	
			if(well %in% names(isoGates)) {
				isoGates[[well]][x[["Channel"]]]
			} else if (x[["Well.Id"]] %in% unique(isoWells$Well.Id)) {
				isoGates[[x[["Well.Id"]]]][x[["Channel"]]]
			}
			else NA
		}) 
	}

	return(data)
})
		
		
setMethod("[[","flowPlate",function(x,i,j,...) {
	if(length(i)!=1)
		stop("subscript out of bounds (index must have length 1)")

	fr <- sampleNames(x@plateSet)[pData(phenoData(x@plateSet))[,"Well.Id"] == i]
	if(is.na(fr)) stop("subscript out of bounds")
	fr <- x@plateSet[[fr]]

})		

setMethod("plateSet","flowPlate",function(fp,...) {
		 return(fp@plateSet)
		})

setMethod("Subset","flowPlate",function(x,subset,...) {
			x@plateSet <- Subset(x@plateSet,subset)
			x
		})

flowPhenoMerge <- function(data, newDF) {
	
	##Assume the columns to be "merged" are the first columns in pData for each annotated data frame
	##Also assume that the ids in the first column of pData(newADF) 
	
	origDF <- pData(phenoData(data))
	
	## Make sure merging columns are unique, otherwise throw an error
	
	newIds <- newDF[,1]
	dataNames <- origDF[,1]
	
	idOrder <- sapply(newIds,function(x) {grep(x,dataNames)})
	
	## Should warn users if ids are missing
	newDF <- newDF[!is.na(idOrder>0),]
	
	## Warn if idOrder is not unique
	idOrder <- unlist(idOrder[!is.na(idOrder>0)])
	
	if(colnames(newDF)[1]=="name") { newDF <- newDF[,-1]}
	
	##Remove flowFrames from data that don't have any config information, also
	##orders the ids
	##Give users a warning that if any frames are removed
	data <- data[sampleNames(data)[idOrder]]
	
	newDF <- cbind(pData(phenoData(data)),newDF)
	
	tempMeta.df <- data.frame(colnames(newDF))
	rownames(tempMeta.df) <- colnames(newDF)
	
	tempPheno.adf <- new("AnnotatedDataFrame", data=newDF, varMetadata=tempMeta.df, dimLabels=c("rowNames", "colNames"))
	sampleNames(tempPheno.adf) <- newDF$name
	
	phenoData(data) <- tempPheno.adf
	
	return(data)
}





makePlateLayout <- function(plateDesc,abName="Ab.Name",sampleType="Sample.Type",negCon="Negative.Control",...) {
	
	##--------------------------------------------------------------------------
	##
	## Function: makePlateLayout
	##
	## Description: Create an wide data.frame for a flowSet from a tall data.frame.    
	##
	##--------------------------------------------------------------------------
	
	
	## Add some validity check function
	
	
	## Now get the unique channels, and remove dashes make them legal column names
	chans <- unique(plateDesc$Channel)
	chans <- gsub("-",".",chans[chans != ""])
	plateDesc$Channel <- gsub("-",".",plateDesc$Channel)
	
	## Make a plate layout data.frame
	plateLayout <- data.frame(Well.Id=as.character(unique(plateDesc$Well.Id)),stringsAsFactors=FALSE)
	
	for(wellChan in chans) {
		temp <- subset(plateDesc,Channel==wellChan,select=c("Well.Id",abName,sampleType,negCon))
		colnames(temp) <- c("Well.Id",wellChan,paste(sampleType,wellChan,sep="."),paste(negCon,wellChan,sep="."))
		plateLayout <- merge(plateLayout,temp,by="Well.Id",all.x=TRUE)
	}
	
	
	temp <- subset(plateDesc,Channel=="",select=c("Well.Id",sampleType))
	plateLayout <- merge(plateLayout,temp,by="Well.Id",all.x=TRUE)
	
	cbind(name=plateLayout$Well.Id,plateLayout)
	
	
}