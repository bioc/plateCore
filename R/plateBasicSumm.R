#########################################################################################################
##
## Filename: plateBasicSumm.R
##
## Author: straine
##
## Data: Apr 4, 2007
##
## Description: Calculate some basic summaries for a plate, Events, MFIs, etc.
##
#########################################################################################################

######################################################################################################################
##
## Function: plateBasicSumm
##    
## Description:  Takes a plate data set (data), rectangleGates for: 1) postive cells, 2) normal morphology cells, 3) non debris
##	events, and the  channels used for analysis.  Returns a data frame summarizing the plate, with percent positive cells,
##	MFIs of positive/negative, etc.
##
######################################################################################################################

plateBasicSumm <- function(data,posCells,normCells,nonDeb,channels,fsc="FSC.A",ssc="SSC.A",bySize=0.01,minDiff=0.25,...) {

	## Column headings for the summary data frame.
	colHeadings <- c("Total.Events","NonDebris.Events","NormalMorhp.Events","FSC.A.MFI","SSC.A.MFI")
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".MFI",sep="")))))
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".pp",sep="")))))
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".posMFI",sep="")))))
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".negMFI",sep="")))))

	## Function list used to calculate values for columns listed above.
	funList <- vector("list",length=7) 
	
	## Get the number of total events
	funList[[1]] <- function(fs,normCells,posCells,nonDeb) {lapply(nonDeb,function(x) length(x@subSet))}
	## Number of non debris events
	funList[[2]] <- function(fs,normCells,posCells,nonDeb) {lapply(nonDeb,function(x) sum(x@subSet))}
	## Number of normal morphology events
	funList[[3]] <- function(fs,normCells,...) {lapply(normCells,function(x) sum(x@subSet))}
	## Median values for fsc, ssc, and other channels
	funList[[4]] <- function(fs,...) {
		lapply(c(fsc,ssc,channels),function(y) {
			fsApply(fs, function(x) { 
				 round(median(x@exprs[,y]),digits=3) 
			})
		})
	}
	## Calculate percent of cells that are positive
	funList[[5]] <- function(fs,normCells,posCells,...) {
		lapply(1:length(channels),function(y) {
			lapply(1:length(fs),function(x) {
				if(class(posCells[[x]][[y]])=="logicalFilterResult") {
					round(sum(posCells[[x]][[y]]@subSet)/sum(normCells[[x]]@subSet),digits=2)*100
				} else { NA }
			})
		})
	}	
	## MFI of positive cells
	funList[[6]] <- function(fs,normCells,posCells,...) {
		lapply(1:length(channels), function(y) {
			lapply(1:length(fs),function(x) {
				if(class(posCells[[x]][[y]])=="logicalFilterResult")  {
					round(median(fs[[x]]@exprs[posCells[[x]][[y]]@subSet,channels[y]]),digits=3)
				} else { NA }
			})
		})
 	}
 	## MFI of negative cells
	funList[[7]] <- function(fs,normCells,posCells,...) {
		lapply(1:length(channels), function(y) {
			lapply(1:length(fs),function(x) {	
				if(class(posCells[[x]][[y]])=="logicalFilterResult")  {
					round(median(fs[[x]]@exprs[!(posCells[[x]][[y]]@subSet),channels[y]]),digits=3)
				} else { NA }
			})
		})
 	}	

	## apply the function list
	plateSumm.df <- data.frame(matrix(unlist(lapply(funList,function(x) x(data,normCells,posCells,nonDeb))),
					nrow=length(data)),stringsAsFactors=FALSE)

	plateSumm.df <- cbind(pData(phenoData(data))$Well.Id,plateSumm.df)
	
	names(plateSumm.df) <- c("Well.Id",colHeadings)
	
	rownames(plateSumm.df) <- as.character(plateSumm.df$Well.Id)
		
	numPops <- callPops(data,channels,bySize=bySize,minDiff=minDiff)
	colnames(numPops) <- unlist(lapply(channels,paste,".pops",sep=""))
	
	plateSumm.df <- cbind(plateSumm.df ,numPops)
	
	name <- sampleNames(data)
	
	plateSumm.df <- cbind(plateSumm.df,name)
		
}