# xyplotFlowPlate.R
# 
# TODO: Add comment
#
# Author: straine
###############################################################################

prepanel.xyplot.flowplate <- 
		function(x, 
				frames, channel.x, channel.y,
				...)
{
	if (length(nm <- as.character(x)) > 1)
		stop("must have only one flow frame per panel")
	if (length(nm) == 1)
	{
		xx <- evalInFlowFrame(channel.x, frames[[nm]])
		yy <- evalInFlowFrame(channel.y, frames[[nm]])
		list(xlim = range(xx, finite = TRUE),
				ylim = range(yy, finite = TRUE),
				dx = diff(xx), dy = diff(yy))
	}
	else list()
}

panel.xyplot.flowplate <-
		function(x, 
				frames,
				channel.x, channel.y,
				channel.x.name, channel.y.name, 
				filter = NULL,
				filterResults = NULL,
				displayFilter = TRUE,
				pch, smooth,
				...)
{
	x <- as.character(x)
	if (length(x) > 1) stop("must have only one flow frame per panel")
	if (length(x) < 1) return()
	
	nm <- x
	xx <- evalInFlowFrame(channel.x, frames[[nm]])
	yy <- evalInFlowFrame(channel.y, frames[[nm]])
	
	this.filter.result <- NULL
	
	groups <- 
			if (!is.null(filterResults))
			{
				filterResults[[nm]]@subSet
			}
			else if (!is.null(filter) && !is.null(filter[[nm]]))
			{
				this.filter.result <- filter(frames[[nm]], filter[[nm]])
				this.filter.result@subSet
			}
			else NULL
	
	if (smooth) {
		panel.smoothScatter(xx, yy, ...)
	}
	else panel.xyplot(xx, yy, pch = pch,
				groups = groups,
				subscripts = seq_along(groups),
				...)
	
	if ((!is.null(filter) && !is.null(filter[[nm]])) && (is.list(displayFilter) || displayFilter))
	{

		display.pars <- list(border = TRUE)
		filter.boundary <-
				filterBoundary(filter = filter[[nm]],
						parameters = c(channel.x.name, channel.y.name),
						frame = frames[[nm]],
						result = this.filter.result)
		do.call(panel.polygon,
				c(filter.boundary, display.pars))
	}
	
	
}

setMethod("xyplot",
		signature(x = "formula", data = "flowPlate"),
		function(x, data, xlab, ylab,
				as.table = TRUE,
				prepanel = prepanel.xyplot.flowplate,
				panel = panel.xyplot.flowplate,
				pch = ".", smooth = TRUE,
				filter = NULL,
				filterResults = NULL,
				displayFilter = TRUE,
				flowStrip=NULL,
				flowStripCex=1,
				strip=function(...,style=1) strip.default(...,style=1),
				...)
		{

			flowData <- plateSet(data)
			pd <- pData(phenoData(flowData))
			uniq.name <- createUniqueColumnName(pd)
			## ugly hack to suppress warnings about coercion introducing
			## NAs (needs to be `undone' inside prepanel and panel
			## functions):
			pd[[uniq.name]] <- factor(sampleNames(flowData)) 
			channel.y <- x[[2]]
			channel.x <- x[[3]]
			if (length(channel.x) == 3)
			{
				channel.x <- channel.x[[2]]
				x[[3]][[2]] <- as.name(uniq.name)
				x[[2]] <- NULL
			}
			else
			{
				x[[3]] <- as.name(uniq.name)
				x[[2]] <- NULL
			}
			channel.x.name <- expr2char(channel.x)
			channel.y.name <- expr2char(channel.y)
			channel.x <- as.expression(channel.x)
			channel.y <- as.expression(channel.y)
			
			if (missing(xlab)) xlab <- channel.x.name
			if (missing(ylab)) ylab <- channel.y.name
			
			if(!missing(filterResults) && is.character(filterResults)) {
				if(filterResults=="Negative.Control") {		
					data <- overlay(data,type="Negative.Control",channel=channel.y.name)
					flowData <- plateSet(data)
					filterResults <- data@overlayFilterResults
				}
			}

			if(is.character(filter) && (filter=="Isogate" || filter=="Negative.Control")) {

				isoGates <- subset(data@wellAnnotation,Channel==channel.y.name)

				filter <- lapply(pd[[uniq.name]],function(x) {
						thresh <- subset(isoGates,Well.Id==pd[x,"Well.Id"] & plateName==pd[x,"plateName"],select=Negative.Control.Gate)[[1]]
						if(!(pd[x,"Well.Id"] %in% isoGates$Well.Id)) {
							return(NULL)
						} else {
							xx <- evalInFlowFrame(channel.x, (flowData@frames)[[as.character(x)]])
							xx <- range(xx, finite = TRUE)
							yy <- evalInFlowFrame(channel.y, (flowData@frames)[[as.character(x)]])
							yy <- range(yy, finite = TRUE)
							yy[[1]] <- thresh
							if(thresh>yy[[2]]) yy[[2]] <- thresh
							minrect <- c(xx[1],yy[1])
							maxrect <- c(xx[2],yy[2])
							names(minrect) <- names(maxrect) <- c(channel.x.name,channel.y.name)
							return(new("rectangleGate",filterId="rectangleGate",
								parameters=c(channel.x.name,channel.y.name),min=c(minrect[1],minrect[2]),max=c(maxrect[1],maxrect[2])))
						}
					})
				names(filter) <- pd[[uniq.name]]
				
			} else if(!is.null(filter)) {
				filter <- lapply(pd[[uniq.name]],function(x) filter)
				names(filter) <- pd[[uniq.name]]
			}

			if(!is.null(flowStrip)) {
				
				if(!("Well.Id" %in% flowStrip)) flowStrip <- c("Well.Id",flowStrip)
				
				## Assumes as.factor in formula x is Well.id
				labels <- subset(data@wellAnnotation,Channel==channel.y.name,select=flowStrip)
				
				if("MFI" %in% flowStrip) labels$MFI <- round(as.numeric(labels$MFI),digits=2)
				if("MFI.Ratio" %in% flowStrip) labels$MFI.Ratio <- round(as.numeric(labels$MFI.Ratio),digits=2)
				if("Percent.Positive" %in% flowStrip) {
					labels$Percent.Positive <- round(as.numeric(labels$Percent.Positive),digits=0)
					labels$Percent.Positive <- sapply(labels$Percent.Positive,function(x) paste(x,"%",sep=""))
				}
				
				labels <- apply(labels,1,function(x) {
					temp <- x[[1]]
					if(length(x)>1) for(i in 2:length(x)) {temp <- paste(temp,x[[i]],sep=" : ")}
					temp
				})
				
				strip=strip.custom(factor.levels=labels,par.strip.text=list(cex=flowStripCex))
			}
			
			
			densityplot(x, data = pd, 
					
					prepanel = prepanel,
					panel = panel,
					
					frames = flowData@frames,
					channel.x = channel.x,
					channel.y = channel.y,
					channel.x.name = channel.x.name,
					channel.y.name = channel.y.name,
					
					filter = filter,
					filterResults = filterResults,
					displayFilter = displayFilter,
					as.table = as.table,
					
					xlab = xlab,
					ylab = ylab,
					pch = pch, smooth = smooth,
					strip=strip,
					...)
		})





