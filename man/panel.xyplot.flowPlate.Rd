\name{panel.xyplot.flowPlate}
\alias{panel.xyplot.flowPlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lattice-flowViz style panel function for flowPlate xyplot. }
\description{
   This function should not be called directly, use \code{\link{xyplot}}.
}
\usage{
panel.xyplot.flowPlate(x, 
				frames,
				channel.x, channel.y,
				channel.x.name, channel.y.name, 
				filter = NULL,
				filterResults = NULL,
				displayFilter = TRUE,
				pch, smooth,
				wellAnnotation = NULL,
				col = superpose.symbol$col,
				...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ character }
  \item{frames}{ flowFrames }
  \item{channel.x}{}
  \item{channel.y}{}
  \item{channel.x.name}{}
  \item{channel.y.name}{}
  \item{filter}{}
  \item{filterResults}{}
  \item{displayFilter}{}
  \item{pch}{} 
  \item{smooth}{}
  \item{wellAnnotation}{}
  \item{col}{}
  \item{\dots}{optional arguments}
}

\author{ Errol Strain }

\seealso{ See Also \code{\link{xyplot}} }
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

