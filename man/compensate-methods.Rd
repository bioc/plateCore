\name{compensate}
\alias{compensate}
\alias{compensate-method}
\alias{compensate,flowPlate,ANY-method}

\docType{methods}

\title{ Compensate a flowPlate to correct for the effects of spillover. }
\description{
  Flow samples are often stained with multiple types of fluorophores. 
  Unfortunately, the emission spectra for these different fluorophores
  often overlap, and the signals must be corrected before proceeding with
  the analysis. Compensate adjusts for spillover using the 
  method implemented in the package flowCore. Unlike flowCore, compensate 
  only adjusts for the dyes/fluorophores listed in wellAnnotation.
  }
\usage{
\S4method{compensate}{flowPlate,ANY}(x, spillover)
}

\arguments{
  \item{x}{A flowPlate }
  \item{spillover}{ The compensation matrix where the row and column names match the 
  fluorescence channels of the flowPlate. }
}
\value{
  Returns a compensated \code{flowPlate}.

}
\author{ Errol Strain }
\seealso{ See Also \code{\link[flowCore::compensation-class]{compensation-class}} }
\examples{
library(plateCore)
data(plateCore)

# Create the compensation matrix
comp.mat <- spillover(x=compensationSet,unstained=sampleNames(compensationSet)[5],
patt=".*H",fsc="FSC-H",ssc="SSC-H",method="median")

## Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
pbmcPlate <- Subset(pbmcPlate, rectGate)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# apply the compensation matrix
fp <- compensate(fp,comp.mat)
}

\keyword{ methods }
