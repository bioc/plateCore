#########################################################################################################
##
## Filename: compFunctions.R
##
## Author: straine
##
## Data: Apr 16, 2007
##
## Description: Functions associated creating compensation matrices.  Basically these are used to preprocess
## 	the compensation flowSet before sending it to the spillover function in flowcore.
##
#########################################################################################################

#########################################################################################################
##
## Function: getCompGates
## 
## Description: Create a hash of rectangle gates representing the positive cells in each of the stained 
##	channels for a compensation tube.  These gates are "not" gates, I find the 99th percentile of the 
##	unstained cells in data (flowSet consisting of just the unstained compensation tube).  
##
#########################################################################################################
getCompGates <- function(data,dyeCols,cutOff=0.99,fsc="FSC.A") {
	
	## Create a hash for storing the rectangle gates.
	compHash <- new.env(hash=TRUE,parent=emptyenv())

	lapply(1:length(dyeCols),function(x) {
		thresh <- quantile(unlist(data@exprs[,dyeCols[x]]),probs=c(cutOff))
		rangeF <- range(unlist(data@exprs[,dyeCols[x]]))
		assign(dyeCols[[x]],!new("rectangleGate",filterId="rectangleGate",parameters=dyeCols[x],min=rangeF[1],max=thresh), env=compHash)
	})
	
	rangeF <- range(unlist(data@exprs[,fsc]))
	assign("unstain",new("rectangleGate",filterId="rectangleGate",parameters=fsc,min=rangeF[1],max=rangeF[2]), env=compHash)

	return(compHash)
}

#########################################################################################################
##
## Function: plateCompensate
##
## Description:  Compensate only for channels which have a dye
##
#########################################################################################################

plateCompensate <- function(plateSet,comp.mat,chanCols) {

	dyeCols <- unlist(lapply(chanCols,function(x) {paste(x,".dye",sep="")}))
	
	plateSet <- fsApply(plateSet,function(x) {
		fileName <- attributes(x)$descriptio[["$FIL"]]
		dyeCols <- pData(phenoData(plateSet))[fileName,dyeCols]!=""
		if(sum(dyeCols)>=2) {
			x <- compensate(x,comp.mat[chanCols[dyeCols],chanCols[dyeCols]])
			
		} else {
			x
		}
	})	
	
	return(plateSet)
}
	