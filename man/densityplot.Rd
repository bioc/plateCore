\name{densityplot}

\alias{densityplot}

\title{ One-dimensional density plots for flowPlates }



\description{
	This function is a modified version of densityplot from
	the flowViz package that allows for multiple flowFrames
	per panel. The flowViz densityplot plots the density curves in a 
	shingle layout, while the flowPlate densityplot overlays
	the densities.
}

\usage{
densityplot(formula,flowSet)(
    x,
    data, 
    xlab,
    prepanel=prepanel.densityplot.flowPlate,
    panel = panel.densityplot.flowPlate,
    filter=NULL,
    ...)
}

\arguments{
\item{x} { A formula describing the layout of the plots.}
\item{data} { A flowPlate. }
\item{filter} { A flowCore filter to apply to each flowFrame. Alternatively,
				filter can take the character string "Isogate" of "Negative.Control"
				when creating overlay plots.}
Other arguments are identical to densityPlot from flowViz.
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
