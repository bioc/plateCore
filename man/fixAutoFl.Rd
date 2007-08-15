\name{fixAutoFl}
\alias{fixAutoFl}

\title{Correct for the effects of cell size (FSC) on autofluorescence.}
\description{
  The 'fixAutoFl' function uses the method of Hahne et al. 2006 (Genome Biology) to fit
  a robust regression to the fluorescence channel of interest on forward scatter (FSC). 
}
\usage{
fixAutoFl(flowSet,chanCols,unstain,fsc="FSC.A") 
}

\arguments{
  \item{flowSet}{ A flowSet object from flowCore. }
  \item{chanCols}{ Selected channels to correct for autofluorescence. }
  \item{unstain}{ Name(s) of the unstained samples. }
  \item{fsc}{ Name of the FSC parameter. }

}

\value{
	Returns a flowSet with autofluorescence due to cell size (FSC) removed.
	}

\author{ Errol Strain }
\examples{
fixAutoFl(flowSet,chanCols=c("FITC.A","PE.A"),unstain='Sample_A01.fcs') 
}
\keyword{methods}
