################################################################################
##
## Filename: flowPlate-utilities.R
##
## Author: Errol Strain
##
## Date: Jan 30, 2007
##
## Description: Functions for manipulating annotated data frames that
##  are associated with flowSets, along with a 
##
################################################################################

################################################################################
## flowPhenoMerge is used to merge a data.frame with the annotated data.frame
## that describes a flowSet. Merge is on the names column, which should be 
## unique (otherwise the creation of the flowSet would fail).
################################################################################
flowPhenoMerge <- function(data, newDF) {
	
	## Get the data.frame currently associated with the flowSet
	origDF <- pData(phenoData(data))
	
	## Merging columns are in the first position
	newIds <- newDF[,1]
	dataNames <- origDF[,1]
	
	## Well Ids are equal to sample names when a flowPlate is first created.
	## If flowPhenoMerge is being used in fpbind, then the unique names
	## have already been created.
	if(identical(newDF$Well.Id,newDF$name)) {
		idOrder <- sapply(newDF$Well.Id,function(x) {grep(x,dataNames)})
	} else {
		idOrder <- sapply(newIds,function(x) {which(dataNames==x)})
	}
	
	## Remove IDs that don't match from the new data.frame
	newDF <- newDF[!is.na(idOrder>0),]
	
	## Create a numeric vector of id order
	idOrder <- unlist(idOrder[!is.na(idOrder>0)])
	
	## If the new data.frame contains a name column, which suggests that this
	## flowPhenoMerge call was the result of fpbind, drop the names as they
	## will be overwritten.
	if(colnames(newDF)[1]=="name") { newDF <- newDF[,-1]}
	
	## Put the flowSet in order
	data <- data[sampleNames(data)[idOrder]]
	
	## Create an annotated data.frame from the new data.frame (newDF)
	tempMeta.df <- rbind(varMetadata(phenoData(data)),data.frame(labelDescription=colnames(newDF)))
	rownames(tempMeta.df) <- c(rownames(varMetadata(phenoData(data))),colnames(newDF))

	## Bind the newDF to the original data.frame to add the sample names
	newDF <- cbind(pData(phenoData(data)),newDF)
	
	tempPheno.adf <- new("AnnotatedDataFrame", data=newDF, varMetadata=tempMeta.df, dimLabels=c("rowNames", "colNames"))
	sampleNames(tempPheno.adf) <- sampleNames(data)
	
	## Replace the annotated data.frame for the flowSet
	phenoData(data) <- tempPheno.adf

	return(data)
}

################################################################################
## makePlateLayout takes the "tall" data frame in wellAnnotation and creates
## a "wide" data frame used associated with the flowSet. In pData data.frame
## of a flowSet, all the info about a well is contained on one row. In the "tall"
## format, each row is one well/one channel.
################################################################################
makePlateLayout <- function(plateDesc,abName="Ab.Name",sampleType="Sample.Type",negCon="Negative.Control",...) {
		
	## Getting rid of "no visible binding errors" in CHECK
	name <- ""
	Well.Id <- ""
	plateName <- ""
	Channel <- ""

	## Add some validity check function
	if(!("name" %in% colnames(plateDesc))) 	plateDesc <- cbind(name=plateDesc$Well.Id,plateDesc)
	else plateDesc <- cbind(name=unlist(plateDesc$name),subset(plateDesc,select=-name))
	
	## Now get the unique channels, and remove dashes make them legal column names
	chans <- unique(plateDesc$Channel)
	chans <- gsub("-",".",chans[chans != ""])
	plateDesc$Channel <- gsub("-",".",plateDesc$Channel)
	
	## Make a plate layout data.frame
	plateLayout <- data.frame(name=as.character(unique(plateDesc$name)),stringsAsFactors=FALSE)
	
	## Add the well ids
	plateLayout$Well.Id <- unlist(lapply(plateLayout$name,function(x) {
					subset(plateDesc,name==x,select=Well.Id)[1,]
					}))
	
	## Get the plate name for each specific name
	plateLayout$plateName <- unlist(lapply(plateLayout$name,function(x) {
						subset(plateDesc,name==x,select=plateName)[1,]
					}))
	
	## Add the channel columns
	for(wellChan in chans) {
		temp <- subset(plateDesc,Channel==wellChan,select=c("name",abName,sampleType,negCon))
		colnames(temp) <- c("name",wellChan,paste(sampleType,wellChan,sep="."),paste(negCon,wellChan,sep="."))
		plateLayout <- merge(plateLayout,temp,by="name",all.x=TRUE)
	}
	
	
	temp <- subset(plateDesc,Channel=="",select=c("name",sampleType))
	plateLayout <- merge(plateLayout,temp,by="name",all.x=TRUE)

	## Add the row and column ids
	plateLayout$Row.Id <- sapply(plateLayout$Well.Id,function(x) {
						substr(x,1,1)
					})
	plateLayout$Column.Id <- sapply(plateLayout$Well.Id,function(x) {
						substr(x,2,3)
					})
				
	
	plateLayout <- plateLayout[order(plateLayout$Well.Id),]
	
	return(plateLayout)

}