\name{flowPlate}
\alias{flowPlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a flowPlate }
\description{
  Constructor for a \code{\link{flowPlate-class}} object. \code{sampleNames} for the \code{flowSet} 
  should match the Well.Id column of \code{wellAnnotation}. Well.Ids must be unique to \code{sampleNames},
  which is usually ensured by using the 3 character designations for wells (e.g. "A01","A02",...,"H12").
}
\usage{
flowPlate(data,wellAnnotation,plateName,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{\code{flowSet} object to be made into a \code{flowPlate}}
  \item{wellAnnotation}{\code{data.frame} describing the layout and contents of the \code{flowPlate}}
  \item{plateName}{Name of the \code{flowPate}, should be unique within the set \code{flowPlates} under consideration}
  \item{\dots}{optional arguments}
}

\value{
  Returns a \code{flowPate} object.
}
\author{ Errol Strain }
\examples{
library(plateCore)
data(plateCore)

## Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
normGate <- norm2Filter("SSC-H","FSC-H",scale.factor=1.5)
pbmcPlate <- Subset(pbmcPlate, rectGate & normGate)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")
}
\keyword{ methods }
