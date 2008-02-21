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

################################################
## This method creates a virtual plate by combining 2 or more flowPlates.
## If the sample names are not unique, then they are changed using make.unique.
## ##############################################
setMethod("fpbind",signature("flowPlate","flowPlate"),function(p1,p2,...) {
			
			## Get the number of arguments
			na <- nargs()
			argl <- list(p1,p2,...)
			
			## Get rid of anything that's not a flowPlate'
			while(na > 0 && is.null(argl[[na]])) { argl <- argl[-na]; na <- na - 1 }
			argl <- argl[sapply(argl,class)=="flowPlate"]
			
			## For now, don't bind two plates with the same name, use subsetting instead
			plateIds <- unlist(sapply(argl,function(x) unique(x@wellAnnotation$plateName))	)
			if(length(plateIds) != length(unique(plateIds))) stop("Binding flowPlates with identical plateNames is not supported.")
			
			## This method can only bind flowPlates where the well annotation data frames have the same column headers
			## I think I'll change this later, since well annotation changes over the course of an analysis (eg results columns are added)
			annl <- lapply(argl,function(x) x@wellAnnotation)	
			for(i in 2:length(argl)) {
				if(!all.equal(colnames(annl[[1]]),colnames(annl[[i]]))) stop("flowPlate Well Annotations must have identical column names.")
			}
	
			## The sample names have to be unique in order for the phenoData on the flowSet to work correctly.
			## Can't have rows with identical names
			sampNames <- unlist(lapply(argl,sampleNames))
			
			uniqNames <-make.unique(sampNames)
			
			sampLength <- lapply(argl,function(x) length(sampleNames(x)))
			
			uniqList <- vector("list",length=length(sampLength))
			index <- 1
			for(i in 1:length(sampLength)) {
				uniqList[[i]] <- uniqNames[index:(index-1+sampLength[[i]])]
				names(uniqList[[i]]) <- sampleNames(argl[[i]])
				index <- index + sampLength[[i]]		
			}

			## Remake the annotation list with the unique names
			annl <- lapply(1:length(argl),function(i) {			
						anndf <- argl[[i]]@wellAnnotation
						anndf$name <- uniqList[[i]][unlist(anndf$name)]
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
			np <- flowPlate(as(frames,"flowSet"),wellAnnotation)

			return(np)
		})

#########################################################
### Get groups, which are currently just the negative control groups.  These consist of the
### negative control well and all the associated test wells.
#########################################################
#setMethod("getGroups",signature("flowPlate"),function(data,type="Negative.Control",chan,...) {
#			
#	wellIds <- unique(data@wellAnnotation$Negative.Control)
#	wellIds <- wellIds[wellIds %in% pData(phenoData(data@plateSet))$Well.Id]
#	
#	wells <- list()
#	
#	if(type=="Negative.Control") {
#		wells <- lapply(wellIds,function(x) {
#			wells <- unlist(subset(data@wellAnnotation,Channel==chan & Negative.Control==x,select=name))
#			wells <- c(unlist(subset(data@wellAnnotation,Channel==chan & Well.Id==x,select=name)),wells)
#			wells <- unique(wells)	
#			if(!length(wells)) NA
#			else wells
#		})
#	}
#	
#	wells <- wells[!is.na(wells)]
#	return(wells)
#})


#########################################################
## This is used to transform the flowSet and  the isotype gates.  Other types
## of gates should be included later. 
#########################################################
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

######################################################
## Constructor for flowPlates
######################################################
setMethod("flowPlate",signature("flowSet"),function(data,wellAnnot,plateName="",...) {
			
			temp <- new("flowPlate")
			
			## Get rid of dashes in flowSet because they're really annoying for lattice
			data <- fsApply(data,function(x) {
						newNames <- gsub("-",".",colnames(exprs(x)))
						colnames(exprs(x)) <- newNames
						x
					})
			
			## Add plateName to well annotation
			wellAnnot <- data.frame(wellAnnot,plateName=plateName,stringsAsFactors=FALSE)		
			## Create a data.frame that corresponds to phenoData.  All info about a well is on a single row
			config <- makePlateLayout(wellAnnot)
	
			##Add the phenoData to the flowSet
			data <- flowPhenoMerge(data,config)
		
			## When a flowPlate is created from a plateLayout template, the names initially missing and need to be added.  Names correspond
			## to filenames
			if(colnames(wellAnnot)[1]=="Well.Id") {
				wellAnnot$name <- sapply(wellAnnot$Well.Id,function(x) pData(phenoData(data))[pData(phenoData(data))$Well.Id==x,"name"])
			} 
			
			## Get rid of the dashes again
			wellAnnot$Channel <- gsub("-",".",wellAnnot$Channel)
			## Only keep annotation for samples that are in  the dataset
			wellAnnot <- subset(wellAnnot,name %in% sampleNames(data))
			
			temp@plateName <- plateName
			temp@wellAnnotation <- wellAnnot
			temp@plateSet <- data
			
			return(temp)
		})

######################################################################
## There are sometimes weird fluorescence values in the flowFrame.  This function removes them in a naive fashion.
## I should change it so it only affects the channels of interest, not time, etc.
######################################################################
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

######################################################################
## Fit linear model to FSC vs channel of interest.  Correct for autofluoresence and then set the median back to it's orignal
## value.  The fitting is performed on log  transformed data.
######################################################################
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
			
			## Apply the correction to the channels of interest in chanCols
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

######################################################################
## Compensate, but only for dyes that were added.  I shouldn't need to make a copy of the flowSet, but I did in case it 
## would save me trouble with environments hanging around.
######################################################################
setMethod("compensate", signature(x="flowPlate"), function(x,spillover) {

	temp <- x@plateSet@frames
	
	frames <- lapply(ls(temp),function(fileName) {
			well <- pData(phenoData(x@plateSet))[fileName,]		
			dyeCols <- colnames(spillover)[!is.na(well[,colnames(spillover)])]
			
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


######################################################################
######################################################################
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
		
		isoGates <- lapply(unique(isoWells$name), function(x) {
			sapply(isoWells[isoWells$name==x,"Channel"], function(i) {	
				mfi <- median(exprs(data[[x]])[,i])
				mfi.mad <- mad(exprs(data[[x]])[,i])
				isoMad <- numMads
				while(quantile(exprs(data[[x]])[,i],probs=0.99)>(mfi+isoMad*mfi.mad)) {
					isoMad <- isoMad + 0.05
				}
				thresh <- mfi+isoMad*mfi.mad	
				thresh
			})		
		 })
		names(isoGates) <- unique(isoWells$name)

		data@wellAnnotation$Isogate <- apply(data@wellAnnotation,1,function(x) {

			well <- subset(data@wellAnnotation, Well.Id== x[["Negative.Control"]] & plateName == x[["plateName"]],select=name)[1,][[1]]

			if(length(well) && well %in% names(isoGates)) {
				isoGates[[well]][x[["Channel"]]]
			} else if (x[["name"]] %in% unique(isoWells$name)) {
				isoGates[[x[["name"]]]][x[["Channel"]]]
			} else NA
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

##############################
## Some convience functions
##############################

setMethod("phenoData","flowPlate",function(object) {
			return(phenoData(object@plateSet))
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

setMethod("wellAnnotation","flowPlate",function(fp,...) {
			return(fp@wellAnnotation)
		})

####################################
## subsetting methods
####################################
setMethod("[",c("flowPlate"),function(x,i,j,...,drop=FALSE) {
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
