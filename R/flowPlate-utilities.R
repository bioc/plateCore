######################################################################################################################
##
## Filename: flowPlate-functions.R
##
## Author: straine
##
## Date: Jan 30, 2007
##
## Description: Extends flowCore to work on plates.  Based on the flowSet-functions.R file.  Functions to read in plate
## 	data, find isotype gates, etc.  
##
######################################################################################################################


colorramp = function (col)  {
	coord <- as.data.frame(t(col2rgb(col))/255)
	x <- seq(0, 1, length = length(col))
	r <- approxfun(x, coord$red)
	g <- approxfun(x, coord$green)
	b <- approxfun(x, coord$blue)
	function(n) {
		x <- seq(0, 1, length = n)
		rgb(r(x), g(x), b(x))
	}
}


flowPhenoMerge <- function(data, newDF) {
	
	##Assume the columns to be "merged" are the first columns in pData for each annotated data frame
	##Also assume that the ids in the first column of pData(newADF) 
	
	origDF <- pData(phenoData(data))
#	origDF$name <- sampleNames(data)
	
	## Make sure merging columns are unique, otherwise throw an error
	
	newIds <- newDF[,1]
	dataNames <- origDF[,1]
	
	## Well Ids are equal to sample names on the first pass
	if(identical(newDF$Well.Id,newDF$name)) {
		idOrder <- sapply(newDF$Well.Id,function(x) {grep(x,dataNames)})
	} else {
		idOrder <- sapply(newIds,function(x) {grep(x,dataNames)})
	}
	
	## Should warn users if ids are missing
	newDF <- newDF[!is.na(idOrder>0),]
	
	## Warn if idOrder is not unique
	idOrder <- unlist(idOrder[!is.na(idOrder>0)])
	
	if(colnames(newDF)[1]=="name") { newDF <- newDF[,-1]}
	
	##Remove flowFrames from data that don't have any config information, also
	##orders the ids
	##Give users a warning that if any frames are removed
	data <- data[sampleNames(data)[idOrder]]
	
	tempMeta.df <- rbind(varMetadata(phenoData(data)),data.frame(labelDescription=colnames(newDF)))
	rownames(tempMeta.df) <- c(rownames(varMetadata(phenoData(data))),colnames(newDF))

	newDF <- cbind(pData(phenoData(data)),newDF)
	
	tempPheno.adf <- new("AnnotatedDataFrame", data=newDF, varMetadata=tempMeta.df, dimLabels=c("rowNames", "colNames"))
	sampleNames(tempPheno.adf) <- sampleNames(data)
	
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
	
	if(!("name" %in% colnames(plateDesc))) 	plateDesc <- cbind(name=plateDesc$Well.Id,plateDesc)
	else plateDesc <- cbind(name=unlist(plateDesc$name),subset(plateDesc,select=-name))
	
	## Now get the unique channels, and remove dashes make them legal column names
	chans <- unique(plateDesc$Channel)
	chans <- gsub("-",".",chans[chans != ""])
	plateDesc$Channel <- gsub("-",".",plateDesc$Channel)
	
	## Make a plate layout data.frame
	plateLayout <- data.frame(name=as.character(unique(plateDesc$name)),stringsAsFactors=FALSE)

	plateLayout$Well.Id <- unlist(lapply(plateLayout$name,function(x) {
					subset(plateDesc,name==x,select=Well.Id)[1,]
					}))
	
	plateLayout$plateName <- unlist(lapply(plateLayout$name,function(x) {
						subset(plateDesc,name==x,select=plateName)[1,]
					}))
	
	for(wellChan in chans) {
		temp <- subset(plateDesc,Channel==wellChan,select=c("name",abName,sampleType,negCon))
		colnames(temp) <- c("name",wellChan,paste(sampleType,wellChan,sep="."),paste(negCon,wellChan,sep="."))
		plateLayout <- merge(plateLayout,temp,by="name",all.x=TRUE)
	}
	
	
	temp <- subset(plateDesc,Channel=="",select=c("name",sampleType))
	plateLayout <- merge(plateLayout,temp,by="name",all.x=TRUE)

	
	plateLayout$Row.Id <- sapply(plateLayout$Well.Id,function(x) {
						substr(x,1,1)
					})
	plateLayout$Column.Id <- sapply(plateLayout$Well.Id,function(x) {
						substr(x,2,3)
					})
				
	
	plateLayout <- plateLayout[order(plateLayout$Well.Id),]

}