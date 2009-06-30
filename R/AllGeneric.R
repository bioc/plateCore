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

##setGeneric("compensate", function(x, spillover, ...)
##			standardGeneric("compensate"))

##setGeneric("Subset",function(x,subset,select=NULL,...) standardGeneric("Subset"))

setGeneric("flowPlate", function(data,wellAnnotation,plateName,...) standardGeneric("flowPlate"))

##setGeneric("setRange", function(x,minF,maxF,type="truncate") standardGeneric("setRange"))

setGeneric("fixAutoFl", function(fp,fsc="FSC.A",chanCols,unstain=NULL,...) standardGeneric("fixAutoFl"))

setGeneric("setControlGates", function(data,gateType="Negative.Control",threshType="MAD",numMads=5,isoquantile=.995,...) standardGeneric("setControlGates"))

setGeneric("plateSet", function(fp,...) standardGeneric("plateSet"))

setGeneric("wellAnnotation", function(fp,...) standardGeneric("wellAnnotation"))

setGeneric("applyControlGates", function(data,gateType="Negative.Control",...) standardGeneric("applyControlGates"))

setGeneric("summaryStats", function(data,...) standardGeneric("summaryStats"))

setGeneric("getGroups", function(data,type="Negative.Control",chan,...) standardGeneric("getGroups"))

setGeneric("fpbind", function(p1,p2,...) standardGeneric("fpbind"))

setGeneric("plotPlate", function(fp,x=NA,method="median",main,col,values,
				width=1,na.action="zero",...) standardGeneric("plotPlate"))

setGeneric("gutterPlot", function(fp,chans=c("FSC-H","SSC-H","FL1-H","FL2-H","FL3-H","FL4-H"),...) standardGeneric("gutterPlot"))

setGeneric("setGate",function(fp,chan,wells,newGate,type="Negative.Control.Gate",...) standardGeneric("setGate"))

setGeneric("adjustGateLog10",function(fp,wells,chan,xlim,ylim,
				fsc="FSC.A",asFact="Well.Id",numEvents=250,type="Negative.Control.Gate",...) standardGeneric("adjustGateLog10"))
			
setGeneric("xyplot", function(x, data, xlab, ylab,
		as.table = TRUE,
		prepanel = prepanel.xyplot.flowPlate,
		panel = panel.xyplot.flowPlate,
		pch = ".", smooth = TRUE,
		filter = NULL,
		filterResults = NULL,
		displayFilter = TRUE,
		flowStrip=NULL,
		flowStripCex=1,
		strip=function(...,style=1) strip.default(...,style=1),
		...)
standardGeneric("xyplot"))

setGeneric("densityplot", function(x, data, xlab,
				prepanel=prepanel.densityplot.flowPlate,
				panel = panel.densityplot.flowPlate,
				as.table=TRUE, filterResult=NULL, ...)
			standardGeneric("densityplot"))
