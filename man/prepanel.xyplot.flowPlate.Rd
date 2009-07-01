\name{prepanel.xyplot.flowPlate}
\alias{prepanel.xyplot.flowPlate}

\title{ Lattice-flowViz style panel function for flowPlate xyplot.  }
\description{
     This function should not be called directly, use \code{\link{xyplot}}. 
}
\usage{
prepanel.xyplot.flowPlate(x, frames, channel.x, channel.y,\dots)
}

\arguments{
  \item{x}{ Character }
  \item{frames}{flowFrames }
  \item{channel.x}{Character string for channel name.}
  \item{channel.y}{Character string for channel name.}
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

