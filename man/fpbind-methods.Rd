\name{fpbind}
\alias{fpbind}
\alias{fpbind-method}
\alias{fpbind,flowPlate,flowPlate-method}

\docType{methods}

\title{ Merge multiple flowPlates into a single virtual flowPlate }
\description{
  A function to combine multiple \code{flowPlates} into a single \code{flowPlate} object. The 
  \code{plateName} identifiers for the \code{flowPlates} must be unique within the set to be
  bound, otherwise the bind will fail.
}
\usage{
fpbind(p1, p2, ...)
}
\arguments{
  \item{p1}{ First \code{flowPlate} }
  \item{p2}{ Second \code{flowPlate} }
  \item{\dots}{ Additional \code{flowPlates} }
}

\value{
 Returns a \code{flowPlate}
}
\author{ Errol Strain }
\examples{
library(plateCore)
data(plateCore)

# Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
pbmcPlate <- Subset(pbmcPlate, rectGate)

# Create a flowPlate object from the platePBMC and the wellAnnotation
fp1 <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")
fp2 <- flowPlate(pbmcPlate,wellAnnotation,plateName="P2")

# Combine the plates.
virtPlate <- fpbind(fp1,fp2)

}

\keyword{ methods }
