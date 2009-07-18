################################################################################
##
## Specialized functions for plotting flowPlates. The motivation for including
##  these two methods is that the flow cytometry group at Amgen has found it
##  useful to visualize data using these approaches. plotPlate makes simple row
##  by column plot of the plate with colored wells. gutterPlot shows how many
##  events in each well are at either their min or max values (too many events
##  at the min/max may indicate a problem).
##
## Author: Errol Strain, Jon Gosink (and Florian Hahne for plotPlate from prada)
##
################################################################################

################################################################################
## colorramp is used by the plotPlate function to get the color gradient
################################################################################
colorramp <- function (col)  {
	coord <- as.data.frame(t(col2rgb(col))/255)
	x <- seq(0, 1, length = length(col))
	r <- approxfun(x, coord$red)
	g <- approxfun(x, coord$green)
	b <- approxfun(x, coord$blue)
	function(n) {
		x <- seq(0, 1, length = n)
		rgb(r(x), g(x), b(x))
	}
}

################################################################################
## plotPlate creates a row-column plot with circles representing wells. Wells are
## colored according to some value of choice (MFI, etc). This code was copied
## from prada (Florian Hahne) and modified for flowPlates.
################################################################################
setMethod("plotPlate",signature("flowPlate"),function(fp,x=NA,method="median",main,col,values,
						width=2,na.action="zero",...) {
	
	## Getting rid of "no visible binding errors" in CHECK
	name <- ""
	Row.Id <- ""
	Column.Id <- ""
					
	ncol = length(unique((pData(phenoData(fp)))[,"Column.Id"]))		
	nrow = length(unique((pData(phenoData(fp)))[,"Row.Id"]))		
					
	nrwell <- ncol*nrow

	if(missing(values) & !is.na(x) & method=="mahalanobis" & all(x %in% plateSet(fp)@colnames)){
		mat <- fsApply(plateSet(fp), each_col, median)
		mat <- mat[,x]
		mat.cov <- cov.rob(mat)
		mat.mean <- apply(mat, 2, mean)
		values <- mahalanobis(mat, mat.cov$center, mat.cov$cov, method="mcd")
	} else if(missing(values) & !is.na(x) & length(x) == 1 & x %in% plateSet(fp)@colnames & method %in% colnames(fp@wellAnnotation)) {
		temp <- fp@wellAnnotation[fp@wellAnnotation$Channel==x,c("name",method)]
		values <- temp[,2]
		names(values) <- temp[,1]		
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
			rownames(subset(pData(phenoData(fp)),Column.Id==col & Row.Id==row))[1]		
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


################################################################################
# gutterPlot creates a plot showing what proportion of events in a well are at
# either their minimum or maximum values (i.e. "in the gutter"). 
################################################################################
setMethod("gutterPlot",signature("flowPlate"),function(fp,chans=c("FSC-H","SSC-H","FL1-H","FL2-H","FL3-H","FL4-H"),...) {

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


################################################################################
# mfiPlot creates a plot showing the mfi ratio versus the percent positive cells
################################################################################
setMethod("mfiPlot",signature("flowPlate"),function(fp,thresh=2,Sample.Type="Test",...) {
			
			## Declaring variables for R CMD check
			Gate.Score = ""
			
			mfiDf <- subset(wellAnnotation(fp),Sample.Type==Sample.Type)
			mfiDf2 <- subset(mfiDf,abs(Gate.Score)<thresh)
			mfiDf3 <- subset(mfiDf,abs(Gate.Score)>=thresh)
			mfiDf$LogMFI.Ratio = log10(mfiDf$MFI.Ratio)
			mfiDf$PosCount <- round(mfiDf$Percent.Positive,0)
			mfiDf$NegCount <- 100-mfiDf$PosCount 
			robMFI <- glmrob(cbind(PosCount,NegCount) ~ LogMFI.Ratio, data=mfiDf,
					family=binomial(link="logit"))
			
			
			mfiRange = range(log10(mfiDf2$MFI.Ratio),na.rm=TRUE)
			
			x=seq(mfiRange[1],mfiRange[2],0.1)
			
			x2 <- -robMFI$coefficients[[1]]-robMFI$coefficients[[2]]*x
			
			
			plot(mfiDf2$MFI.Ratio,mfiDf2$Percent.Positive, log="x",...)
			points(mfiDf3$MFI.Ratio,mfiDf3$Percent.Positive,col="red",...)
			lines(10^x,100/(1 + exp(x2)),lwd=2,col='red')
					
			
		})