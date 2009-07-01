################################################################################
##
## Filename: AllClasses.R
##
## Author: Errol Strain
##
## Data: Apr 27, 2007
##
## Description: Class for analyzing plate based flow cytometry experiments.
##
################################################################################

setClass("flowPlate",                   
	representation(
		plateName="character",
		plateSet="flowSet",
		plateFilters="filterSet",
		wellAnnotation="data.frame",
		fjWorkspaces="environment"),
	prototype=list(
		plateName=character(0),
		plateConfig=data.frame(),
		wellAnnotation=data.frame(),
		fjWorkspaces=emptyenv())
)


