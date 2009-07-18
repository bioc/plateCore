\name{mfiPlot}
\alias{mfiPlot}
\alias{mfiPlot-method}
\alias{mfiPlot,flowPlate-method}

\docType{methods}

\title{ mfiPlot }
\description{
  A Quality Control plot that shows the MFI Ratio versus the percentage of 
  positive cells in a flowPlate. The robust logistic regression is performed
  using gmlrob from the robustbase package.
}
\usage{
mfiPlot(fp, thresh=2, Sample.Type="Test",Events="Percentage", \dots)
}

\arguments{
  \item{fp}{ A \code{flowPlate}. }
  \item{thresh}{ Points more than "thresh" number of standard deviations away from the
  best fit line will be colored red. }
  \item{Sample.Type}{ Type of sample to show on plot. Defaults to "Test" } 
  \item{Events}{ The robust logistic regression can be performed using either 
  the percentage of events above the negative control gate ("Percentage") or the
  actual number of events above the gate ("Actual"). }
  \item{\dots}{optional arguments to plot and points.}
}

\value{
  Creates a plot where the x-axis is MFI Ratio and the
  y-axis is the percentage of cells above the negative control gate.
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

## Create an MFI plot
mfiPlot(fp,thresh=2.5,xlab="MFI Ratio (Test MFI / Isotype MFI)",xlim=c(0.1,250),
		ylab="Percentage of cells above the isotype gate",pch=23)

}

\keyword{ methods }

