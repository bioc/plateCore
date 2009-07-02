\name{flowPlate-class}
\docType{class}
\alias{flowPlate}
\alias{flowPlate-class}
\alias{flowPlate,flowSet-method}
\alias{[,flowPlate,ANY-method}
\alias{[[,flowPlate,ANY-method}
\alias{sampleNames,flowPlate-method}
\alias{phenoData,flowPlate-method}

\title{ Create a flowPlate }
\description{
  Constructor for a \code{\link{flowPlate}} object. \code{sampleNames} for the \code{flowSet} 
  should match the Well.Id column of \code{wellAnnotation}. Well.Ids must be unique to \code{sampleNames},
  which is usually ensured by using the 3 character designations for wells (e.g. "A01","A02",...,"H12").
}
\usage{
flowPlate(data,wellAnnotation,plateName,\dots)
}

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
pbmcPlate <- Subset(pbmcPlate, rectGate)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

## Subset the flowPlate, creating another flowPlate
fpSmall <- fp["A01"]

## Extract a flowFrame from a flowPlate
ff <-fp[["A01"]]

## Retrieve sample names from flowPlate
sampNames <- sampleNames(fp)

## Retrieve the annotatedDataFrame describing the flowPlate
adf <- phenoData(fp)
}
\keyword{ methods }
