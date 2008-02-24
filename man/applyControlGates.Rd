\name{applyControlGates}
\alias{applyControlGates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Once the function setControlGates() has been used to create gates for a flowPlate object, then apply the gates to test
samples using this function. }
\description{
  The applyControlGates() function is separate from setControlGates since the gates may need be adjusted
  outside of setControlGates.  
}
\usage{
applyControlGates(data, gateType="Negative.Control", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A flowPlate dataset. }
  \item{gateType}{ ~~Describe \code{\dots} here~~ }
}
\value{
  \item{data }{A flowPlate where the wellAnnotation now contains additional 
  columns correspondin to }
}
\author{ Errol Strain }

\seealso{ See Also \code{\link{setControlGates}}, ~~~ }
\examples{
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

fp <- setControlGates(fp,gateType="Negative.Control")
fp <- applyControlGates(fp,gateType="Negative.Control")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ flowPlate }

