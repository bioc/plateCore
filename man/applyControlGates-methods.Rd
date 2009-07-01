\name{applyControlGates}
\alias{applyControlGates}
\alias{applyControlGates-methods}
\alias{applyControlGates,flowPlate-method}

\docType{methods}

\title{ Apply control gates to a flowPlate }
\description{
  Once the function setControlGates has been used to create gates for a flowPlate object, then apply the gates to test
samples using this function. The applyControlGates function is separated from setControlGates since gates may need be adjusted
  outside of setControlGates.  
}
\usage{
applyControlGates(data, gateType="Negative.Control", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A flowPlate dataset. }
  \item{gateType}{ The type of gate to be applied to the flowPlate. Currently only "Negative.Control" gates are supported. }
  \item{\dots}{optional arguments}
}
\value{
  Returns a \code{flowPlate} where the wellAnnotation now contains additional 
  columns corresponding to the percentage of cells above background (Percent.Positive) and the number
  of positive cells (Positive.Count). 
}
\author{ Errol Strain }

\seealso{ See Also \code{\link{setControlGates}} }
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
fp <- applyControlGates(fp,gateType="Negative.Control")

# Percent Positive and Counts columns are now in the wellAnnotation
head(wellAnnotation(fp))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

