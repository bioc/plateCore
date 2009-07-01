\name{\%on\%}
\docType{methods}
\alias{\%on\%}
\alias{\%on\%-methods}
\alias{\%on\%,ANY,flowPlate-method}


\title{Methods for Function \%on\% in Package `plateCore'}
\description{
This operator is used to construct a \code{transformFilter} that first applies a \code{transformList} to the data before
applying the \code{filter} operation to a \code{flowPlate}. 
}
\author{Errol Strain}
\examples{
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

## Create a rectangle filter
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))

xyplot(`FL1-H` ~ `FSC-H` | as.factor(name), 
	transform("FL1-H"=log10) \%on\% fp, smooth=FALSE, filter=rectGate, displayFilter=FALSE)

}
\keyword{methods}

