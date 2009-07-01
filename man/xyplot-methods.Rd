\name{xyplot}
\alias{xyplot}
\alias{xyplot-method}
\alias{xyplot,formula,flowPlate-method}

\docType{methods}

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
\S4method{xyplot}{formula,flowPlate}(x, data, xlab, ylab,
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
\item{x}{A formula describing the layout of the plots. Plots for flowPlates
	usually condition on either \code{as.factor(name)} or
	\code{as.factor(Well.Id)} since only one flowFrame can be shown on each
	panel (with the exception of Negative.Control overlays).}
\item{data}{A flowPlate.}
\item{xlab}{Label for x-axis.}
\item{ylab}{Label for y-axis.}
\item{as.table}{Defaults to table layout.}
\item{prepanel}{Lattice-flowViz prepanel function.}
\item{panel}{Lattice-flowViz panel function.}
\item{pch}{Plotting character.}
\item{smooth}{Plot a smoothed scatterplot by default.}
\item{filter}{A flowCore filter to apply to each flowFrame.}
\item{filterResults}{If filterResults="Negative.Control", the negative control
	wells corresponding to a test well are overlayed in the test well plots.}
\item{displayFilter}{Defaults to displaying filter on the plot.}
\item{flowStrip}{Character vector indicating additional information to be printed
	on the strip. Values can include any combination of "Well.Id","MFI","MFI.Ratio", and
	"Percent.Positive".}
\item{flowStripCex}{Font size for the flowStrip.}
\item{strip}{Lattice strip function.}
\item{\dots}{Optional arguments}
}
\seealso{
\code{\link[flowViz:xyplot]{flowViz::xyplot}}
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
