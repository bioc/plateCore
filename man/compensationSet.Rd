\name{compensationSet}
\alias{compensationSet}
\docType{data}
\title{ Sample Compensation Data Set}
\description{
  Sample Compensation Data Set
}
\usage{data(plateCore)}
\format{
  The format is an object of class \code{flowSet} composed of 5
  \code{flowFrames}. The \code{flowSet} consists of 4 stained 
  and one unstained \code{flowFrames}. Peripheral Blood Mononucleocytes
  (PBMCs) were stained with FITC (Fluorescein isothiocyanate), PE (phycoerythrin),
  PerCp (Peridinin-chlorophyll), and APC (Allophycocyanin).
}

\source{
  Sample data set from BD FACS CAP analysis.
}
\author{ Errol Strain }
\seealso{ See Also \code{\link[flowCore:compensation-class]{compensation-class}}} 
\examples{
library(plateCore)
data(plateCore)

# Create the compensation matrix
comp.mat <- spillover(x=compensationSet,unstained=sampleNames(compensationSet)[5],
patt=".*H",fsc="FSC-H",ssc="SSC-H",method="median")

}

\keyword{datasets}
