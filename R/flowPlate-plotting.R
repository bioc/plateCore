# TODO: Add comment
#
# Author: straine
###############################################################################

setMethod("plotPlate",signature("flowPlate"),function(fp,x=NA,method="median",main,col,values,
						width=2,na.action="zero",...) {
					
	ncol = length(unique((pData(phenoData(fp)))[,"Column.Id"]))		
	nrow = length(unique((pData(phenoData(fp)))[,"Row.Id"]))		
					
	nrwell <- ncol*nrow

	if(missing(values) & !is.na(x) & method=="mahalanobis" & all(x %in% plateSet(fp)@colnames)){
		mat <- fsApply(plateSet(fp), each_col, median)
		mat <- mat[,x]
		mat.cov <- cov.rob(mat)
		mat.mean <- apply(mat, 2, mean)
		values <- mahalanobis(mat, mat.cov$center, mat.cov$cov, method="mcd")
	} else if(missing(values) & !is.na(x) & length(x) == 1 & x %in% plateSet(fp)@colnames) {
		values <- fsApply(plateSet(fp),function(ff) {
			 eval(parse(text=paste(method,"(exprs(ff)[,\"",x,"\"])",sep="")))
		})[,1]
	} else if (missing(values) & !is.na(x) & x == "events"){
		values <- fsApply(plateSet(fp),function(ff) {
					nrow(exprs(ff))
				})[,1]
	} else if(!missing(values)) {
		if(!is.numeric(values) || !is.vector(values) || length(values)!=nrwell)
			stop("'values' must be a numeric vector of length 'ncol*nrow'")		
	} else {
		stop("x is not valid")
	}
	
	## Put the data in order and check for missing values
	colIds <- unique((pData(phenoData(fp)))[,"Column.Id"])
	rowIds <- unique((pData(phenoData(fp)))[,"Row.Id"])
	
	sampNames <- unlist(lapply(rowIds,function(row) {
		lapply(colIds,function(col){
			subset(pData(phenoData(fp)),Column.Id==col & Row.Id==row, select=name)[1,1]		
		})
	}))
	
	tempValues <- rep(NA,length(sampNames))
	names(tempValues) <- sampNames
	
	tempValues[names(values)] <- values
	
	values <- tempValues

	
	valuesRange=range(values, na.rm=TRUE)
	
		
	## user coordinates: x=(-0.5...13.5), y=(-0.5...9.5)
	xlim   = c(-0.5, ncol+1.5)
	ylim   = c(-0.5, nrow+1.5)
	colbarwid = 0.3
	fw     = diff(xlim)+colbarwid
	fh     = diff(ylim)
	
	height <- width/fw * fh
	args <- list(width=width, height=height)
			
	layout(matrix(1:2, ncol=2), widths=c(diff(xlim), colbarwid), heights=1)
	
	## device coordinates
	u2px = function(x) (x-xlim[1]) / fw * width
	u2py = function(y) (y-ylim[1]) / fh * height
	
	par(mai=c(0,0,0,0))
	cex = 1.5
	plot(x=0, y=0, type="n", bty="n", xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlim=xlim, ylim=ylim)
	graphics::text((1:ncol), 0, unique(pData(phenoData(fp))[,"Column.Id"]), adj=c(0.5,0), cex=cex)
	graphics::text(0, (nrow:1), unique(pData(phenoData(fp))[,"Row.Id"]), adj=c(0, 0.5), cex=cex)
	if(!missing(main))
		graphics::text((ncol+1)/2, nrow+1, main, adj=c(0.5, 1), cex=cex)
	
	nrcolors   = 256
	thepalette = colorramp(col)(nrcolors)
	
	# the mapping from values to color indices
	z2icol <- function(z) {
		res = round((z-valuesRange[1])/diff(valuesRange)*(nrcolors-1))+1
		res[res > nrcolors] = nrcolors
		res[res < 1       ] = 1
		return(res)
	}
	icol2z <- function(i) {
		(i-1)/(nrcolors-1)*diff(valuesRange)+valuesRange[1]
	}
	stopifnot(all(z2icol(icol2z(1:nrcolors))==1:nrcolors))
	circcol <- thepalette[z2icol(values)]
	
	## circles
	radius = 0.45
	xc = radius*cos(seq(0, 2*pi, len=73))
	yc = radius*sin(seq(0, 2*pi, len=73))
	x0 = 1     + (0:(nrwell-1)) %% ncol
	y0 = nrow - (0:(nrwell-1)) %/% ncol
	
	switch(na.action,
			zero = {
				circcol[is.na(circcol)] <- thepalette[z2icol(0)]
				wh <- 1:nrwell
			}, 
			omit = {
				wh <- which(!is.na(circcol))
			},
			stop(paste("Invalid value of 'na.action':", na.action))
	)
	
	for (i in wh) {
		polygon(x = x0[i]+xc,
				y = y0[i]+yc,
				col=circcol[i])
	}
	
	xmin = 0.5
	xmax = ncol + 0.5
	ymin = 0.5
	ymax = nrow + 0.5
	polygon(c(xmin, xmax, xmax, xmin, xmin), c(ymin, ymin, ymax, ymax, ymin))
	
	par(mai=0.5*c(1,0,1,0))
	image(x=0, y=icol2z(1:nrcolors), z=matrix(1:nrcolors, nrow=1),
			col = thepalette, xaxt="n")   
	
	x0 = 1 + (wh-1) %%  ncol
	y0 = 1 + (wh-1) %/% ncol
	dx = dy = 0.4
	x1 = u2px(x0-dx)
	x2 = u2px(x0+dx)
	y1 = u2py(y0-dy)
	y2 = u2py(y0+dy)
	
	res <- list(which=wh, coord=floor(cbind(x1, y1, x2, y2) + 0.5))
	invisible(res)
})

setMethod("gutterPlot",signature("flowPlate"),function(fp,chans=c("FSC.H","SSC.H","FL1.H","FL2.H","FL3.H","FL4.H"),...) {

	resultMat <- fsApply(plateSet(fp),function(x) {
		apply(exprs(x)[,chans],2,function(y) {
			if(min(y) < max(y)) {
				(sum(y==max(y)) + sum(y=min(y)))/length(y)
			} else {
				sum(y=min(y))/length(y)
			}		
		})
	})

	# Create the main plot window
	plot(1, 0.5, type="n", xlim=c(1,nrow(resultMat)), ylim=c(0,1),
			xaxt="n",
			xlab="Well ID",
			ylab="% Events Pegged Full or Min Scale", 
			main=fp@plateName, ...)
	
	for (i in chans) {
		points(1:nrow(resultMat), resultMat[,i], type="b", pch=which(chans==i), col=which(chans==i))
	}
	
	axis(1, at=1:nrow(resultMat), labels=pData(phenoData(fp))$Well.Id);
	legend(x="topleft", legend=chans, cex=1, 
			bty="n", pch=1:length(chans), col=1:length(chans));
	
})

setMethod("timePlateMap",signature("flowPlate"),function(fp,parameter=1,trim=c(0.05, 0.95),...) {
	
	ncol = length(unique((pData(phenoData(fp)))[,"Column.Id"]))		
	nrow = length(unique((pData(phenoData(fp)))[,"Row.Id"]))				

	nrwell <- ncol*nrow
						
#	def.par <- par(no.readonly = TRUE) # save default, for resetting...
	# Create a 8 X 12 layout of plots
#	mat = textPlateMap(1:nrwell)
	nf  = layout(matrix(1:nrwell,nrow=nrow,ncol=ncol, byrow=TRUE))
	layout.show(nf)
	par(mar=c(0,0,0,0))
	
#	plateMap   = setUpPlateMap(aFlowSet)
	rangeMat   = fsApply(plateSet(fp), each_col, range)
	quantRange = fsApply(plateSet(fp), each_col, quantile, probs=trim)

	# Really, just take the middle of the 5th and 95th percentiles
	yLim     = quantile(quantRange[,parameter])[c(2,4)]
	xLim     = as.numeric(range(rangeMat[,"Time"]))
	
	goodIndex = 1
	for (i in 1:length(plateMap$wellMapVec)) {
		if (plateMap$wellMapVec[i] == TRUE) {
			
			# Note this bit of indirection is necessary as the flowSet wells may
			# not be in the "correct" order or even of the length of a full plate.  
			# However, the wellMapVec will be in correct order.
			j = plateMap$plotOrder[goodIndex]      
			x = as.vector(aFlowSet[[j]]@exprs[, "Time"])
			y = as.vector(aFlowSet[[j]]@exprs[, parameter])
			goodIndex = goodIndex + 1
		} else {
			x = mean(xLim)
			y = mean(yLim)
		}
		
		plot(x, y, pch=".", xlim=xLim, ylim=yLim,
				xaxt="n", yaxt="n",col="blue")
		# And plot a smoothed line over the data
		loSmooth = lowess(x,y, f=0.66)
		lines(loSmooth, col="black", lwd=2)
	}
	
	# Finish up the plate figure
#	par(def.par)
	plotTitle   = fp@plateName
	if (is.numeric(parameter)) {
		plotTitle = paste(plotTitle, "by", names(aFlowSet[[1]])[parameter])
	}
	
	mtext(plotTitle, side=3)
	mtext("Time", side=1, cex=0.75)
})

