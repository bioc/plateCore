\name{gutterPlot}
\alias{gutterPlot}
\alias{gutterPlot-method}
\alias{gutterPlot,flowPlate-method}

\docType{methods}

\title{ gutterPlot }
\description{
  A Quality Control plot to check the number of events in each channel that are at 
  either their minimum or maximum value.  A large number of these events may indicate
  a problem with the sample.
}
\usage{
gutterPlot(fp, chans = c("FSC-H", "SSC-H", "FL1-H", "FL2-H", "FL3-H", "FL4-H"), \dots)
}

\arguments{
  \item{fp}{ A \code{flowPlate}. }
  \item{chans}{ Channels of interest to show on the \code{gutterPlot}. }
  \item{\dots}{optional arguments}
}

\value{
  Creates a plot where the x-axis is the different wells in a flowPlate, and the
  y-axis is the fraction of events at the boundary.
}

\author{ Jon Gosink and Errol Strain }

\examples{
library(plateCore)
data(plateCore)

## Create a flowPlate
fp <- flowPlate(pbmcPlate,wellAnnotation,"p1001")

gutterPlot(fp,chans=c("FSC-H","SSC-H","FL1-H","FL2-H","FL3-H","FL4-H"))
}

\keyword{ methods }

