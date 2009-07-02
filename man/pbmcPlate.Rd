\name{pbmcPlate}
\alias{pbmcPlate}
\docType{data}
\title{ pbmcPlate Data Set}
\description{
  One 96-well plate from a BD FACS CAP analysis of Peripheral Blood Mononucleocyte (PMBC) cells.
}
\usage{data(plateCore)}
\format{
  The format is an object of class \code{flowSet} composed of 96
  \code{flowFrames}. Each \code{flowFrame} corresponds to one well from the plate.

}
\details{
  BD FACS CAP (\url{http://www.bd.com/technologies/discovery_platform/BD_FACS_CAP.asp})
  is a platform for screening a large number of antibodies (200+) on human samples.
  Antibodies are arrayed 3-per well on a 96-well plate, along with the appropriate
  controls. In this experiment, an early version of FACS CAP was used to screen
  PBMCs from 2 donors for 189 different human cell surface markers. The complete
  data set is available from the FICCS site shown below. The \code{pbmcPlate} include
  with plateCore is a lymphocyte enriched subset of one of the replicate plates 
  for donor 1. 
}
\source{
  Complete dataset available at
  \url{http://www.ficcs.org/software.html#Data_Files}, the Flow Informatics and Computational
  Cytometry Society website (FICCS)  
}
\references{
  Add reference for plateCore paper (when/if published).
}

\keyword{datasets}
