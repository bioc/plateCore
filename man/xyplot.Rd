\name{xyplot}
\alias{xyplot}

\title{ Scatter plots (dotplots) for flowPlates. }

\description{
	A function to create dotplots, and smoothed scatter plots,
	from flowPlates. This function is a slightly modified
	version of xyplot from flowViz.  The flowPlate xyplot allows users to 
	overlay plots of test samples versus controls, and makes creating
	informative flowStrips easier. Refer to the documentation for
	xyplot from flowViz and lattice for more detailed information.
}

\usage{
xyplot(x, data, xlab, ylab,
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
}
\arguments{
\item{x}{A formula describing the layout of the plots.}
\item{data}{A flowPlate.}
\item{filter}{A flowCore filter to apply to each flowFrame. Alternatively,
				filter can take the character string "Isogate" of "Negative.Control"
				when creating overlay plots.}
Other arguments are identical to xyplot from flowViz.
}

\value{
  Returns a \code{flowPlate}
}
\author{ Errol Strain }
\examples{
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

## Create a rectangle filter
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))

xyplot(`SSC-H` ~ `FSC-H` | as.factor(name), 
	fp[1], smooth=FALSE, filter=rectGate, displayFilter=FALSE)
}
\keyword{methods}
\keyword{dplot}
