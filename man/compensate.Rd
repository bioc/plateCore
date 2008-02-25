\name{compensate}
\alias{compensate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compensate a flowPlate to correct for the effects of spillover. }
\description{
  Flow samples are often stained with multiple types fluorophores, and their
  emission spectra can overlap. Compensate adjusts for spillover using the 
  method implemented in the package flowCore. Unlike flowCore, compensate 
  only adjusts for the dyes/fluorophores listed in wellAnnotation.
  }
\usage{
compensate(x, spillover)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ flowPlate }
  \item{spillover}{ The compensation matrix where the row and column names match the 
  fluorescence channels of the flowPlate. }
}
\value{
  Returns a compensated \code{flowPlate}.

}
\author{ Errol Strain }
\seealso{ See Also \code{\link[flowCore:flowCore-package]{compensate}} }
\examples{
library(plateCore)
data(plateCore)

# Create the compensation matrix
comp.mat <- spillover(x=compensationSet,unstained=sampleNames(compensationSet)[5],
patt=".*H",fsc="FSC.H",ssc="SSC.H",method="median")

## Get the lymphocytes
rectGate <- rectangleGate("FSC.H"=c(300,700),"SSC.H"=c(50,400))
normGate <- norm2Filter("SSC.H","FSC.H",scale.factor=1.5)
pbmcPlate <- Subset(pbmcPlate, rectGate & normGate)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# apply the compensation matrix
fp <- compensate(fp,comp.mat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
