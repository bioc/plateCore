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


setMethod("flowPlate",signature("flowSet"),function(data,config=NULL,wellAnnot=NULL,...) {
			
			temp <- new("flowPlate")
			
			## Get rid of dashes in flowSet because they're annoying for lattice
			data <- fsApply(data,function(x) {
						newNames <- gsub("-",".",colnames(exprs(x)))
						colnames(exprs(x)) <- newNames
						x
					})
			
			if(!missing(config)) {
				## Do some validity check on config
				data <- flowPhenoMerge(data,config)
			}
			
			if(!missing(wellAnnot)) {		
				temp@wellAnnotation <- wellAnnot
			}	
			
			temp@plateSet <- data
			return(temp)
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
		wellId <- pData(phenoData(x@plateSet))[fileName,"Well.Id"]
		dyeCols <- subset(x@wellAnnotation,Well.Id==wellId & !is.na(Dye),select=Channel)
		dyeCols <- dyeCols[which(dyeCols$Channel %in% colnames(spillover)),]
		if(length(dyeCols)>=2) {
			y <- compensate(y,spillover[dyeCols,dyeCols])
			
		} else {
			y
		}
	})	
	
	x@plateSet <- plateSet
	return(x)
})
	

setMethod("setContolGates", signature("flowPlate"), function(data,gateType="Isotype",numMads=5,...) {
			
	if(gateType=="Isotype") {
		## First get the control gate for each of the isotype groups.	
		isoWells <- subset(data@wellAnnotation,Sample.Type=="Isotype" & !is.na(Dye))
		
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
		
		## Get the Dye Groups for each isotype, and average any replicates.	
		controlGroups <- isoWells[isoWells$Well.Id == names(isoGates),]
		controlGroups$thresh <- apply(controlGroups,1, function(x) {
				isoGates[[x["Well.Id"]]][x["Channel"]]
			})
	}

	data@wellAnnotation$Isogate <- apply(data@wellAnnotation,1,function(x) {
		mean(subset(controlGroups,Isotype.Group==as.numeric(x[["Isotype.Group"]]) & Channel==x[["Channel"]],select=thresh))
	}) 
	
	return(data)
})
		
		
setMethod("[[","flowPlate",function(x,i,j,...) {
	if(length(i)!=1)
		stop("subscript out of bounds (index must have length 1)")

	fr <- sampleNames(data@plateSet)[pData(phenoData(data@plateSet))[,"Well.Id"] == i]
	if(is.na(fr)) stop("subscript out of bounds")
	fr <- data@plateSet[[fr]]

})		