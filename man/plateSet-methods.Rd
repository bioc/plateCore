\name{plateSet}
\alias{plateSet}
\alias{plateSet-method}
\alias{plateSet,flowPlate-method}

\docType{methods}

\title{ Get the flowSet from a flowPlate object }
\description{
  A function to retrieve the \code{flowSet} from a \code{flowPlate}.
}
\usage{
plateSet(fp, \dots)
}
\arguments{
  \item{fp}{ A \code{flowPlate} }
  \item{\dots}{optional arguments}
}

\value{
 Returns a \code{flowSet}
}
\author{ Errol Strain }

\examples{
library(plateCore)
data(plateCore)

## Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
pbmcPlate <- Subset(pbmcPlate, rectGate)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

## Retrieves the flowSet
fs <- plateSet(fp)
}

\keyword{ methods }
