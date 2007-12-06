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