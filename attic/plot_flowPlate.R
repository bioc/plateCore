# plot_flowPlate.R
# 
# TODO: Add comment
#
# Author: straine
###############################################################################


plotPlate2  = function(x, nrow=8, ncol=12, main, col,channel,
	 width, na.action="zero") {
	
    data <- subset(x@wellAnnotation,Channel==channel)
 
 
	nrow <- length(unique(pData(phenoData(plateSet(x)))[,"Row.Id"]))
	ncol <- length(unique(pData(phenoData(plateSet(x)))[,"Column.Id"]))
	
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
	graphics::text((1:ncol), 0, paste(1:ncol), adj=c(0.5,0), cex=cex)
	graphics::text(0, (nrow:1), c("A","B","C","D","E","F","G","H"), adj=c(0, 0.5), cex=cex)
	if(!missing(main))
		graphics::text((ncol+1)/2, nrow+1, main, adj=c(0.5, 1), cex=cex)
	
	nrcolors   = 256
	thepalette = colorramp(col)(nrcolors)
	
	# the mapping from values to color indices
	z2icol <- function(z) {
		res = round((z-xrange[1])/diff(xrange)*(nrcolors-1))+1
		res[res > nrcolors] = nrcolors
		res[res < 1       ] = 1
		return(res)
	}
	icol2z <- function(i) {
		(i-1)/(nrcolors-1)*diff(xrange)+xrange[1]
	}
	stopifnot(all(z2icol(icol2z(1:nrcolors))==1:nrcolors))
	circcol <- thepalette[z2icol(x)]
	
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
	
	if(!missing(device) && device %in% c("pdf", "png", "jpeg"))
		dev.off()	   
	
	x0 = 1 + (wh-1) %%  ncol
	y0 = 1 + (wh-1) %/% ncol
	dx = dy = 0.4
	x1 = u2px(x0-dx)
	x2 = u2px(x0+dx)
	y1 = u2py(y0-dy)
	y2 = u2py(y0+dy)
	
	res <- list(which=wh, coord=floor(cbind(x1, y1, x2, y2) + 0.5))
	invisible(res)
}
