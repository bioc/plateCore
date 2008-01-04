#########################################################################################################
##
## Filename: AllGeneric.R
##
## Author: straine
##
## Data: Apr 27, 2007
##
## Description:
##
#########################################################################################################


setGeneric("flowPlate", function(data,...) standardGeneric("flowPlate"))

setGeneric("setRange", function(x,minF,maxF,type) standardGeneric("setRange"))

setGeneric("fixAutoFl", function(fp,...) standardGeneric("fixAutoFl"))

setGeneric("setContolGates", function(data,...) standardGeneric("setContolGates"))

setGeneric("plateSet", function(fp,...) standardGeneric("plateSet"))

setGeneric("overlayFilterResult", function(subSet,...) standardGeneric("overlayFilterResult"))

setGeneric("overlay", function(data,...) standardGeneric("overlay"))

setGeneric("applyControlGates", function(data,...) standardGeneric("applyControlGates"))

setGeneric("summaryStats", function(data,...) standardGeneric("summaryStats"))

setGeneric("getGroups", function(data,...) standardGeneric("getGroups"))

setGeneric("pbind", function(p1,p2,...) standardGeneric("pbind"))