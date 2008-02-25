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


setGeneric("flowPlate", function(data,wellAnnotation,plateName,...) standardGeneric("flowPlate"))

setGeneric("setRange", function(x,minF,maxF,type="truncate") standardGeneric("setRange"))

setGeneric("fixAutoFl", function(fp,fsc="FSC.A",chanCols,unstain,...) standardGeneric("fixAutoFl"))

setGeneric("setControlGates", function(data,gateType="Negative.Control",numMads=5,...) standardGeneric("setControlGates"))

setGeneric("plateSet", function(fp,...) standardGeneric("plateSet"))

setGeneric("wellAnnotation", function(fp,...) standardGeneric("wellAnnotation"))

setGeneric("overlayFilterResult", function(subSet,...) standardGeneric("overlayFilterResult"))

setGeneric("overlay", function(data,...) standardGeneric("overlay"))

setGeneric("applyControlGates", function(data,gateType="Negative.Control",...) standardGeneric("applyControlGates"))

setGeneric("summaryStats", function(data,...) standardGeneric("summaryStats"))

setGeneric("getGroups", function(data,...) standardGeneric("getGroups"))

setGeneric("fpbind", function(p1,p2,...) standardGeneric("fpbind"))