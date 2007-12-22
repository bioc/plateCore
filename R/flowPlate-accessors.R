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
		
		if("Isogate" %in% colnames(e2@wellAnnotation)) {
			for(y in e1@transforms) {
				wellList <- which(e2@wellAnnotation$Channel %in% y@input)	
				e2@wellAnnotation[wellList,"Isogate"] <- y@f(e2@wellAnnotation[wellList,"Isogate"])			  	
			}
		}
			
		e2@plateSet <- fsApply(e2@plateSet,"%on%",e1=e1)	
		e2

	})

setMethod("flowPlate",signature("flowSet"),function(data,wellAnnot,...) {
			
			temp <- new("flowPlate")

			## Get rid of dashes in flowSet because they're annoying for lattice
			data <- fsApply(data,function(x) {
						newNames <- gsub("-",".",colnames(exprs(x)))
						colnames(exprs(x)) <- newNames
						x
					})

			config <- makePlateLayout(wellAnnot)
			data <- flowPhenoMerge(data,config)

			wellAnnot$name <- sapply(wellAnnot$Well.Id,function(x) pData(phenoData(data))[pData(phenoData(data))$Well.Id==x,"name"])
			
			wellAnnot$Channel <- gsub("-",".",wellAnnot$Channel)
			temp@wellAnnotation <- wellAnnot
			
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
		

setMethod("applyControlGates", signature("flowPlate"), function(data,gateType="Isogate",...) {
			
		if(gateType=="Isogate") {
			## First get the control gate for each of the isotype groups.	

			data@wellAnnotation$Percent.Positive <- apply(data@wellAnnotation,1,function(x) {
					thresh <- as.numeric(x[["Isogate"]])	
					chan <- x[["Channel"]]
					frame <- data@plateSet[[x[["name"]]]]

					if(!is.na(thresh) && chan %in% colnames(data@plateSet)) {
						iso <- new("rectangleGate",filterId="rectangleGate",parameters=chan,min=thresh,max=Inf)
						isoResult <- filter(frame,iso)
						return(100*(sum(isoResult@subSet)/nrow(exprs(frame))))
					} else { NA }	
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
