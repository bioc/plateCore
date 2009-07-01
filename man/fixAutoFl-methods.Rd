\name{fixAutoFl}
\alias{fixAutoFl}
\alias{fixAutoFl-method}
\alias{fixAutoFl,flowPlate-method}

\docType{methods}

\title{Correct for the effects of cell size (FSC) on autofluorescence}
\description{
  The \code{fixAutoFl} function uses the method of Hahne et al. 2006 (Genome Biology) to fit
  a robust regression to the fluorescence channel of interest versus forward scatter (FSC). 
  The current implementation scales the corrected data so the median fluorescence intensity
  (MFI) is the same before and after \code{fixAutoFl}.
}
\usage{
fixAutoFl(fp,fsc="FSC.A",chanCols,unstain,\dots) 
}

\arguments{
  \item{fp}{ A \code{flowPlate} }
  \item{chanCols}{ Selected channels to correct for autofluorescence. }
  \item{unstain}{ Name(s) of the unstained samples. The function will try to find samples with
  Sample.Type="Unstained" if no names are given. If there are multiple unstained samples the
  function will average the slopes. }
  \item{fsc}{ Name of the FSC parameter. }
  \item{\dots}{optional arguments}
}

\value{
	Returns a flowPlate with autofluorescence due to cell size (FSC) corrected.
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

## Correct for autofluorescence in FL1.H-FL4.H
fp <- fixAutoFl(fp,fsc="FSC-H",chanCols=c("FL1-H","FL2-H","FL3-H","FL4-H")) 
}
\keyword{methods}
