\name{fixAutoFl}
\alias{fixAutoFl}

\title{Correct for the effects of cell size (FSC) on autofluorescence.}
\description{
  The 'fixAutoFl' function uses the method of Hahne et al. 2006 (Genome Biology) to fit
  a robust regression to the fluorescence channel of interest versus forward scatter (FSC). 
  The current implementation scales the corrected data so the median fluorescence intensity
  (MFI) is the same before and after fixAutoFl().
}
\usage{
fixAutoFl(flowPlate,fsc="FSC.A",chanCols,unstain) 
}

\arguments{
  \item{flowSet}{ A flowSet object from flowCore. }
  \item{chanCols}{ Selected channels to correct for autofluorescence. }
  \item{unstain}{ Name(s) of the unstained samples. The function will try to find samples with
  Sample.Type="Unstained" if no names are given. If there are multiple unstained samples the
  function will average the slopes. }
  \item{fsc}{ Name of the FSC parameter. }

}

\value{
	Returns a flowPlate with autofluorescence due to cell size (FSC) corrected.
	}

\author{ Errol Strain }
\examples{
## Load the plateCore sample data
data(plateCore)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(platePBMC,wellAnnotation,plateName="plate1")

## Correct for autofluorescence in FL1.H-FL4.H
fp <- fixAutoFl(fp,chanCols=c("FL1.H","FL2.H","FL3.","FL4.H")) 
}
\keyword{methods}
