#########################################################################################################
##
## Filename: AllClasses.R
##
## Author: straine
##
## Data: Apr 27, 2007
##
## Description: Classes for analyzing plate based flow cytometry experiments.
##
#########################################################################################################

setClass("flowPlate",                   
	representation(
		plateSet="flowSet",
		plateFilters="environment",
		plateConfig="data.frame",
		wellAnnotation="data.frame"),
	prototype=list(
		plateFilters=emptyenv(),
		plateConfig=data.frame(),
		wellAnnotation=data.frame())
)

setGeneric("flowPlate", function(x,...) standardGeneric("flowPlate"))

setMethod("flowPlate",signature("flowSet"),function(x) {
	temp <- new("flowPlate")
	temp@plateSet <- x
	return(temp)
})


