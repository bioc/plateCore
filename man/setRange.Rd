\name{setRange}
\alias{setRange}

\title{ Truncate the signal values on a \code{flowPlate} }
\description{
  Cytometers occasionally output values that are outside the range of the instrument.
  These values are either truncated or removed at the limits set by minF and maxF.
  Currently, this truncates/removes everthing in the \code{exprs} matrix, including
  the "Time" column. Future versions will only truncate a user-defined set of columns.
}
\usage{
setRange(x, minF, maxF, type="truncate")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{flowPlate} }
  \item{minF}{ Minimum signal value }
  \item{maxF}{ Maximum signal value }
  \item{type}{ Either "truncate" to truncate the signals, or "remove" to delete them}
}

\value{
  Returns a \code{flowPlate}
}
\author{ Errol Strain }
\examples{
library(plateCore)
data(plateCore)

## Get the lymphocytes
rectGate <- rectangleGate("FSC.H"=c(300,700),"SSC.H"=c(50,400))
normGate <- norm2Filter("SSC.H","FSC.H",scale.factor=1.5)
pbmcPlate <- Subset(pbmcPlate, rectGate & normGate)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

fp <- setRange(fp,minF=1,maxF=10000)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
