\name{setControlGates}
\alias{setControlGates}
\alias{setControlGates-method}
\alias{setControlGates,flowPlate-method}

\docType{methods}

\title{ Create control gates for a flowPlate }
\description{
  A function to estimate the threshold between positive and negative cells. This threshold corresponds to
  a one-dimensional gate, and cells above the gate are considered positive. The gate is initially set by 
  calculating 5 median absolute deviations above the median. If the value is less than the 99th quantile,
  then the gate is iteratively moved up. The default value of \code{numMads=5} often works well on the 
  linear scale, but will need to be adjusted for transformed data.
}
\usage{
setControlGates(data, gateType, threshType="MAD", numMads=5, isoquantile=.995, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A \code{flowPlate} }
  \item{gateType}{ The type of gate to be set. Currently only "Negative.Control" gates are supported. }
  \item{threshType}{}
  \item{numMads}{ Number of median absolute deviations above the median to set the initial gate. }
  \item{isoquantile}{}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
