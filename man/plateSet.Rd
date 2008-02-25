\name{plateSet}
\alias{plateSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get the flowSet from a flowPlate object }
\description{
  A function to retrieve the \code{flowSet} from a \code{flowPlate}.
}
\usage{
plateSet(fp, ...)
}
\arguments{
  \item{fp}{ A \code{flowPlate} }
}

\value{
 Returns a \code{flowSet}
}
\author{ Errol Strain }

\examples{
library(plateCore)
data(plateCore)

## Get the lymphocytes
rectGate <- rectangleGate("FSC.H"=c(300,700),"SSC.H"=c(50,400))
normGate <- norm2Filter("SSC.H","FSC.H",scale.factor=1.5)
pbmcPlate <- Subset(pbmcPlate, rectGate & normGate)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

## Retrieves the flowSet
fs <- plateSet(fp)
}

\keyword{ methods }
