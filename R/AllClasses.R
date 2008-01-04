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

setClass("overlayFilterResult",
		representation(subSet="numeric"),
		prototype=list(subSet=vector("numeric")))

setClass("flowPlate",                   
	representation(
		plateName="character",
		plateSet="flowSet",
		plateFilters="filterSet",
		wellAnnotation="data.frame",
		overlayFilterResults="list",
		fjWorkspaces="environment"),
	prototype=list(
		plateName=character(0),
		plateConfig=data.frame(),
		wellAnnotation=data.frame(),
		overlayFilterResults=list(),
		fjWorkspaces=emptyenv())
)


