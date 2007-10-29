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
		plateFilters="filterSet",
		wellAnnotation="data.frame",
		fjWorkspaces="environment"),
	prototype=list(
		plateConfig=data.frame(),
		wellAnnotation=data.frame(),
		fjWorkspaces=emptyenv())
)


