\name{panel.xyplot.flowPlate}
\alias{panel.xyplot.flowPlate}

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

\arguments{
  \item{x}{ character }
  \item{frames}{ flowFrames }
  \item{channel.x}{ xyplot channel.x parameter }
  \item{channel.y}{ xyplot channel.y parameter }
  \item{channel.x.name}{ xyplot channel.x.name parameter }
  \item{channel.y.name}{ xyplot channel.y.name parameter }
  \item{filter}{ xyplot filter parameter }
  \item{filterResults}{ xyplot filterResults parameter }
  \item{displayFilter}{ xyplot displayFilter parameter }
  \item{pch}{ xyplot pch parameter } 
  \item{smooth}{ xyplot smooth parameter }
  \item{wellAnnotation}{ wellAnnotation data.frame }
  \item{col}{ xyplot col parameter }
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

\keyword{ methods }

