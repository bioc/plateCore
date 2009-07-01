\name{setControlGates}
\alias{setControlGates}
\alias{setControlGates-method}
\alias{setControlGates,flowPlate-method}

\docType{methods}

\title{ Create control gates for a flowPlate }
\description{
  A function to estimate the threshold between positive and negative cells. This threshold corresponds to
  a one-dimensional gate, and cells above the gate are considered positive.
  The default value of \code{numMads=5} generally works well on the 
  linear scale, but will need to be adjusted for transformed data. If each well
  contains a large number of events for the cell type of interest (>1000), then
  using the 99.5th quantile usually gives similar values.
}
\usage{
setControlGates(data, gateType, threshType="MAD", numMads=5, isoquantile=.995, \dots)
}

\arguments{
  \item{data}{A \code{flowPlate} }
  \item{gateType}{The type of gate to be set. Currently only "Negative.Control" gates are supported. }
  \item{threshType}{Values can be either "MAD", for median absolute deviation
    based gating, or "isoQuant" for quantile based gating. }
  \item{numMads}{Number of median absolute deviations above the median to set the initial gate. }
  \item{isoquantile}{Quantile setting for "isoQuant" threshType. }
  \item{\dots}{optional arguments.}
}

\value{
 Returns a \code{flowPlate}
}

\author{ Errol Strain }

\examples{
library(plateCore)
data(plateCore)

## Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
pbmcPlate <- Subset(pbmcPlate, rectGate)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# Create a set of negative control gates and then apply them
fp <- setControlGates(fp,gateType="Negative.Control")

# There should now be a Negative.Control.Gate column in wellAnnotation
head(wellAnnotation(fp))
}

\keyword{ methods }
