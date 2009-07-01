\name{densityplot}
\alias{densityplot}
\alias{densityplot-method}
\alias{densityplot,formula,flowPlate-method}
\alias{densityplot,formula,flowPlate,ANY,ANY,ANY,missing-method}

\docType{methods}

\title{ One-dimensional density plots for flowPlates }

\description{
	This function is a modified version of densityplot from
	the flowViz package that allows for multiple flowFrames
	per panel. The flowViz densityplot plots the density curves in a 
	one per panel style, while the flowPlate densityplot overlays
	the densities.
}

\usage{
\S4method{densityplot}{formula,flowPlate}(x, data, xlab,
    prepanel=prepanel.densityplot.flowPlate,
    panel = panel.densityplot.flowPlate,
    as.table=TRUE,
    filterResult=NULL,
    \dots)
}

\arguments{
  \item{x}{A formula describing the layout of the plots.}
  \item{data}{A flowPlate.}
  \item{xlab}{Label for the x-axis}
  \item{prepanel}{Lattice-flowViz prepanel function.}
  \item{panel}{Lattice-flowViz panel function.}
  \item{as.table}{Defaults to table layout.}
  \item{filterResult}{ filterResult can either take the character
  	string "Negative.Control" and have the negative control wells added
  	to the panels, or if filterResult is a flowFrame then the density curve
  	for the flowFrame will be added to each panel.}
  \item{\dots}{optional arguments}
Other arguments are identical to densityPlot from flowViz.
}
\seealso{
\code{\link[flowViz:densityplot]{flowViz::densityplot}}
}
\examples{
# Load the plateCore package and data
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# Overlay the first 3 flowFrames.  If the groups argument was
# omitted, then the flowFrames would be combined into a single
# density curve.
densityplot(~ `FSC-H`, fp[1:3], groups=name, auto.key=TRUE)
}

\keyword{methods}
\keyword{dplot}
