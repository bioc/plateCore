\name{summaryStats}
\alias{summaryStats}
\alias{summaryStats-method}
\alias{summaryStats,flowPlate-method}

\docType{methods}

\title{ Compute summary statistics on a flowPlate }
\description{
  This function computes the median fluorescence intensity (MFI) and the
  MFI ratio (ratio of test well MFI to negative control MFI) for each
  well/channel in a \code{flowPlate}. The predicted percent
  positive (Predict.PP) and gate score (Gate.Score) come from a robust
  logistic regression of the MFI ratio to either the percentage of positive cells
  or the actual count of positive cells.
  Predict.PP is the estimated percent positive based on the MFI ratio and
  Gate.Score is the number of standardized residuals the sample data point
  is away from the best fit line. The glmrob function from the robustbase
  package is used for the regression. Results from summaryStats are stored in the
  \code{wellAnnotation data.fame}.
}
\usage{
summaryStats(data, Events="Percentage", \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A \code{flowPlate} }
  \item{Events}{ The robust logistic regression can be performed using either 
  the percentage of events above the negative control gate ("Percentage") or the
  actual number of events above the gate ("Actual"). }
  \item{\dots}{optional arguments}
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
fp <- applyControlGates(fp,gateType="Negative.Control")

# Compute summary statistics
fp <- summaryStats(fp)

# There should now be MFI and MFI.ratio columns in the wellAnnotation
head(wellAnnotation(fp))

}

\keyword{ methods }
