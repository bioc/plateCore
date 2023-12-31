\name{prepanel.densityplot.flowPlate}
\alias{prepanel.densityplot.flowPlate}

\title{ Lattice-flowViz style panel function for flowPlate densityplot. }
\description{
  This function should not be called directly, use \code{\link{densityplot}}.
}
\usage{
prepanel.densityplot.flowPlate(x, frames, channel,\dots)
}

\arguments{
  \item{x}{ Character }
  \item{frames}{ flowFrames }
  \item{channel}{Character string for channel name.}
  \item{\dots}{optional arguments}
}

\author{ Errol Strain }

\seealso{ See Also \code{\link{densityplot}} }
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

\keyword{ methods }

