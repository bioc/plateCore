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


setMethod("fpbind",signature("flowPlate","flowPlate"),function(p1,p2,...) {
			
			na <- nargs()
			argl <- list(p1,p2,...)
			
			while(na > 0 && is.null(argl[[na]])) { argl <- argl[-na]; na <- na - 1 }
			argl <- argl[sapply(argl,class)=="flowPlate"]
			
			plateIds <- sapply(argl,function(x) x@plateName)
			
			if(length(plateIds) != length(unique(plateIds))) stop("Binding flowPlates with identical plateNames is not supported.")
			
			annl <- lapply(argl,function(x) data.frame(plateName=x@plateName,x@wellAnnotation))
			
			annl <- lapply(annl,function(x) { 
					x$Well.Id <- apply(x,1,function(y) {paste(y['Well.Id'],y['plateName'],sep=".")})
					x$Negative.Control <- apply(x,1,function(y) {paste(y['Negative.Control'],y['plateName'],sep=".")})
					x$name<- apply(x,1,function(y) {paste(y['Well.Id'],y['name'],sep="")}) 
					x
					})
			
			for(i in 2:length(argl)) {
				if(!all.equal(colnames(annl[[1]]),colnames(annl[[i]]))) stop("flowPllate Well Annotations must have identical column names.")
			}
					
			wellAnnotation <- annl[[1]]
			
			for(i in 2:length(annl)) {
				wellAnnotation <- rbind(wellAnnotation,annl[[i]])
			}
			
			frames <- unlist(lapply(argl,function(x) {			
					temp <- as(x@plateSet@frames,"list")
					fnames <- pData(phenoData(plateSet(x)))$name
					fIds <- pData(phenoData(plateSet(x)))$Well.Id
					names(temp) <- sapply(1:length(fnames),function(y) {paste(fIds[y],".",x@plateName,fnames[y],sep="")})
					fsList <- lapply(1:length(fnames),function(y){ 
								temp[[y]]@description[["GUID"]] <- names(temp)[y]
								temp[[y]]
					})		
					names(fsList) <- names(temp)
					fsList
				}))
			
			np <- flowPlate(as(frames,"flowSet"),wellAnnotation)

			return(np)
		})


setMethod("getGroups",signature("flowPlate"),function(data,type="Negative.Control",chan,...) {
			
	wellIds <- unique(data@wellAnnotation$Negative.Control)
	wellIds <- wellIds[wellIds %in% pData(phenoData(data@plateSet))$Well.Id]
	
	wells <- list()
	
	if(type=="Negative.Control") {
		wells <- lapply(wellIds,function(x) {
			wells <- unlist(subset(data@wellAnnotation,Channel==chan & Negative.Control==x,select=name))
			wells <- c(unlist(subset(data@wellAnnotation,Channel==chan & Well.Id==x,select=name)),wells)
			wells <- unique(wells)	
			if(!length(wells)) NA
			else wells
		})
	}
	
	wells <- wells[!is.na(wells)]
	return(wells)
})



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

setMethod("flowPlate",signature("flowSet"),function(data,wellAnnot,plateName=character(0),...) {
			
			temp <- new("flowPlate")
			
			temp@plateName=plateName
			
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
			wellAnnot <- subset(wellAnnot,name %in% sampleNames(data))
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



setMethod("summaryStats", signature("flowPlate"), function(data,...) {

	wellAnnotation <- data@wellAnnotation

	wellAnnotation$MFI <- apply(wellAnnotation,1,function(x) {
		chan <- x[["Channel"]]
		frame <- data@plateSet[[x[["name"]]]]

		if(chan %in% colnames(exprs(frame))) {
			return(median(exprs(frame)[,chan]))
		} else { NA }	
	})

	wellAnnotation$MFI.Ratio <- apply(wellAnnotation,1,function(x) {
			negWell <- x[["Negative.Control"]]
			chan <- x[["Channel"]]
			negMFI <- subset(wellAnnotation,Well.Id==negWell & Channel==chan,select=MFI)
			frame <- data@plateSet[[x[["name"]]]]
			
			if(chan %in% colnames(exprs(frame)) && negWell!="") {
				negMFI <- subset(wellAnnotation,Well.Id==negWell & Channel==chan,select=MFI)
				if(!is.na(negMFI) && !is.na(x[["MFI"]])) return(as.numeric(x[["MFI"]])/as.numeric(negMFI))
			} else { NA }	
	})	

	data@wellAnnotation <- wellAnnotation
	
	return(data)

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
			wa <- data@wellAnnotation 
			wa$Percent.Positive <- apply(data@wellAnnotation,1,function(x) {
					thresh <- as.numeric(x[["Isogate"]])	
					chan <- x[["Channel"]]
					frame <- data@plateSet[[x[["name"]]]]

					if(!is.na(thresh) && chan %in% colnames(exprs(frame))) {
						iso <- new("rectangleGate",filterId="rectangleGate",parameters=chan,min=thresh,max=Inf)
						isoResult <- filter(frame,iso)
						return(100*(sum(isoResult@subSet)/nrow(exprs(frame))))
					} else { NA }	
				})	
		
			wa$Total.Count <- apply(data@wellAnnotation,1,function(x) {
					thresh <- as.numeric(x[["Isogate"]])	
					chan <- x[["Channel"]]
					frame <- data@plateSet[[x[["name"]]]]
					
					if(!is.na(thresh) && chan %in% colnames(exprs(frame))) {
						return(nrow(exprs(frame)))
					} else { NA }	
				})	
		
			wa$Positive.Count <- apply(data@wellAnnotation,1,function(x) {
					thresh <- as.numeric(x[["Isogate"]])	
					chan <- x[["Channel"]]
					frame <- data@plateSet[[x[["name"]]]]
					
					if(!is.na(thresh) && chan %in% colnames(exprs(frame))) {
						iso <- new("rectangleGate",filterId="rectangleGate",parameters=chan,min=thresh,max=Inf)
						isoResult <- filter(frame,iso)
						return(sum(isoResult@subSet))
					} else { NA }	
				})	
		
		
			data@wellAnnotation <- wa
		}
		return(data)
	})


	

setMethod("plateSet","flowPlate",function(fp,...) {
		 return(fp@plateSet)
		})

setMethod("Subset","flowPlate",function(x,subset,...) {
			x@plateSet <- Subset(x@plateSet,subset)
			x
		})

setMethod("sampleNames","flowPlate",function(object) {
			sampleNames(phenoData(object@plateSet))
		})

## ==========================================================================
## subsetting methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

setMethod("[",c("flowPlate"),function(x,i,j,...,drop=FALSE) {
		if(missing(drop)) drop = FALSE
		if(missing(i)) return(x)

		if(is.numeric(i) || is.logical(i)) {
				copy = sampleNames(x)[i]
		} else {
			copy = i
			i    = match(i,sampleNames(x))
		}

		x@plateSet <- x@plateSet[copy]
		x@wellAnnotation <- subset(x@wellAnnotation,name %in% copy)
		x
})


setMethod("[[","flowPlate",function(x,i,j,...) {
			if(length(i)!=1)
				stop("subscript out of bounds (index must have length 1)")
			if(i %in% sampleNames(x)) fr <- i
			else fr <- sampleNames(x)[pData(phenoData(x@plateSet))[,"Well.Id"] == i]
			if(is.na(fr)) stop("subscript out of bounds")
			fr <- x@plateSet[[fr]]
			
		})	
