\name{Subset}
\docType{methods}
\alias{Subset}
\alias{Subset-method}
\alias{Subset,flowPlate-method}
\alias{Subset,flowPlate,ANY-method}

\title{ Subset }
\description{
  Select a subset of events in a flowPlate. If a flowPlate
  and filter are supplied, then this function calls the 
  Subset function from flowCore. Additionally, the plateCore
  version of Subset also makes it easy to filter individual
  flowFrames and keep the flowPlate structure.

}
\usage{
\S4method{Subset}{flowPlate,ANY}(x, subset, select=NULL,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{flowPlate} }
  \item{subset}{A filter object }
  \item{select}{An optional vector of either sample names or well ids.}
  \item{\dots}{optional arguments}
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

## Apply the filter only to sample A01. The other flowFrames
## are not filtered.
fp <- Subset(fp,rectGate,'A01')
}
\keyword{ methods }
