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

prepanel.density.flowPlate <- 
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



panel.density.flowPlate <-
		function(x, 
				frames, channel,
				groups=NULL,
				subscripts,
#				darg,
#				col = superpose.symbol$col,
#				col.points = col,
#				pch = superpose.symbol$pch,
#				cex = superpose.symbol$cex,
#				col.line = col,
#				lty = superpose.line$lty,
#				lwd = superpose.line$lwd,
				...)
{

	
#	if (is.null(groups))
#	{
#		nx <- length(x)
#		col.points <- rep(col.points, length = nx)
#		col.line <- rep(col.line, length = nx)
#		pch <- rep(pch, length = nx)
#		cex <- rep(cex, length = nx)
#		lty <- rep(lty, length = nx)
#		lwd <- rep(lwd, length = nx)
#		alpha <- rep(alpha, length = nx)
#	}
#	else
#	{
#		groups <- as.factor(groups)[subscripts]
#		stopifnot(length(groups) == length(x))
#		## goal: make colors etc vectors as before, but
#		## associate by group
#		
#		ng <- nlevels(groups)
#		gcode <- as.numeric(groups)
#		pch <- rep(pch, length = ng)[gcode]
#		cex <- rep(cex, length = ng)[gcode]
#		lty <- rep(lty, length = ng)[gcode]
#		lwd <- rep(lwd, length = ng)[gcode]
#	}
#	browser()
	x <- as.character(x)
	for (i in seq_along(x))
	{
		nm <- x[i]
		xx <- evalInFlowFrame(channel, frames[[nm]])

		panel.densityplot(xx,plot.points=FALSE,
#				darg=darg,
#				col = col[i],
#				cex = cex[i],
#				pch = pch[i],
#				lty = lty[i],
#				lwd = lwd[i],
				...)
	}
}





setMethod("densityplot",
		signature(x = "formula", data = "flowPlate"),
 	function(x, data, xlab,
				prepanel = prepanel.density.flowPlate,
				panel = panel.density.flowPlate,
				as.table = TRUE,
				...)
		{
			
			flowData <- plateSet(data)
			pd <- pData(phenoData(flowData))
			ocall <- sys.call(sys.parent())
			ccall <- match.call(expand.dots = FALSE)
			ccall <- manipulate.call(ocall, ccall)
			uniq.name <- createUniqueColumnName(pd)
			## ugly hack to suppress warnings about coercion introducing NAs
			pd[[uniq.name]] <- factor(sampleNames(flowData))
			channel <- x[[2]]
			if (length(channel) == 3)
			{
				channel <- channel[[2]]
				x[[2]][[2]] <- as.name(uniq.name)
			}
			else x[[2]] <- as.name(uniq.name)
			channel.name <- expr2char(channel)
			channel <- as.expression(channel)
			if (missing(xlab)) xlab <- channel.name
			ccall$x <- x
			ccall$data <- pd
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