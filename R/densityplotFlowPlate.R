##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# densityplotFlowPlate.R
# 
# Author: straine
#
# Date: Apr 30, 2008
#
# Description:
#
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

prepanel.densityplot.flowPlate <- 
		function(x, frames, channel, 
				...)
{

	xx <- eapply(frames, function(ff) {
		range(evalInFlowFrame(channel, ff))
	})

	yy <- eapply(frames, function(ff) {
			   vals <- evalInFlowFrame(channel, ff)  
			   max(do.call(density, list(x = vals))$y)
		})
	list(xlim=range(xx, finite=TRUE),ylim=c(0,max(unlist(yy))))
}



panel.densityplot.flowPlate <-
		function(x, 
				frames, channel, wellAnnotation,
				groups=NULL,
				subscripts,
				col = superpose.symbol$col,
				col.points = col,
				col.line = col,
				filterResult=NULL,
				...)
{

	
	## Getting rid of "no visible binding errors" in CHECK
	name <- ""
	Well.Id <- ""
	plateName <- ""
	Channel <- ""
	
	superpose.symbol <- trellis.par.get("superpose.symbol")

	if (is.null(groups))
	{
		ng <- length(x)+1
		col.points <- rep(col.points, length = ng)
		col.line <- rep(col.line, length = ng)
	}
	else
	{
		groups <- as.factor(groups)[subscripts]
		stopifnot(length(groups) == length(x))
		## goal: make colors etc vectors as before, but
		## associate by group
		
		ng <- nlevels(groups)
		gcode <- as.numeric(groups)
		col.points <- rep(col.points, length = ng)[gcode]
		col.line <- rep(col.line, length = ng)[gcode]
	}
	

	x <- as.character(x)


	
	for (i in seq_along(x))
	{
		nm <- x[i]
		xx <- evalInFlowFrame(channel, frames[[nm]])

		panel.densityplot(xx,data=data,plot.points=FALSE,
				col.line = col.line[i],
				col = col[i],
				...)

		if(!is.null(filterResult) && class(filterResult)=="character" && filterResult=="Negative.Control") {

			panel.abline(v=subset(wellAnnotation,name==nm & Channel==as.character(channel[[1]]))$Negative.Control.Gate)
			nc <- subset(wellAnnotation,name==nm & Channel==as.character(channel[[1]]))$Negative.Control
			ncp <- subset(wellAnnotation,name==nm & Channel==as.character(channel[[1]]))$plateName
			if(nc %in% wellAnnotation$Well.Id) {
				nc <- subset(wellAnnotation,Well.Id==nc & Channel==as.character(channel[[1]]) & plateName==ncp)$name[[1]]
				xx <- evalInFlowFrame(channel, frames[[nc]])
				
				panel.densityplot(xx,data=data,plot.points=FALSE,
						col.line = col.line[ng],
						col = col[ng],
						...)
			}
		} else if (class(filterResult)=="flowFrame") {
			panel.abline(v=subset(wellAnnotation,name==nm & Channel==as.character(channel[[1]]))$Negative.Control.Gate)
			xx <- evalInFlowFrame(channel, filterResult)
	
			panel.densityplot(xx,data=data,plot.points=FALSE,
				col.line = col.line[ng],
				col = col[ng],
				...)		
		}
	}
}





setMethod("densityplot",
		signature(x = "formula", data = "flowPlate"),
 	function(x, data, xlab,
				prepanel = prepanel.densityplot.flowPlate,
				panel = panel.densityplot.flowPlate,
				as.table = TRUE,
				filterResult=NULL,
				...)
		{
			
			flowData <- plateSet(data)
			pd <- pData(phenoData(flowData))
			ocall <- sys.call(sys.parent())
			ccall <- match.call(expand.dots = FALSE)
			ccall <- manipulate.call(ocall, ccall)
			uniq.name <- createUniqueColumnName(pd)
			## ugly hack to suppress warnings about coercion introducing NAs
			pd[[uniq.name]] <- factor(sampleNames(data),
					levels=unique(sampleNames(data))) 
			channel <- x[[2]]
			if (length(channel) == 3)
			{
				channel <- channel[[2]]
				x[[2]][[2]] <- as.name(uniq.name)
			} else x[[2]] <- as.name(uniq.name)
			channel.name <- expr2char(channel)
			channel <- as.expression(channel)
			if (missing(xlab)) xlab <- channel.name

			
			ccall$x <- x
			ccall$data <- pd
			ccall$wellAnnotation <- data@wellAnnotation
			ccall$prepanel <- prepanel
			ccall$panel <- panel
			ccall$as.table <- as.table
			ccall$xlab <- xlab
			ccall$frames <- flowData@frames
			ccall$channel <- channel
			ccall[[1]] <- quote(lattice::densityplot)
			ans <- eval.parent(ccall)
			ans$call <- ocall
			ans
		}
)
