\name{wellAnnotation}
\alias{wellAnnotation}
\alias{wellAnnotation-method}
\alias{wellAnnotation,flowPlate-method}

\docType{methods}

\title{ Retrieve a data.frame describing the content of a flowPlate }
\description{
  wellAnnotation returns the tall \code{data.frame} describing the layout of a 
  \code{flowPate}, where each row corresponds to one well-channel. 
}
\usage{
wellAnnotation(fp, \dots)
}

\arguments{
  \item{fp}{ A flowPlate dataset. }
  \item{\dots}{optional arguments}
}
\value{
  Returns a \code{data.frame}.  
}
\author{ Errol Strain }

\examples{
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# Look at the top of wellAnnotation
head(wellAnnotation(fp))

}
\keyword{ methods }

