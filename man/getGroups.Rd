\name{getGroups}
\alias{getGroups}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve the negative control groups from a \code{\link{flowPlate-class}} }
\description{
   
}
\usage{
getGroups(data,type="Negative.Control",chan, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A flowPlate dataset. }
  \item{type} {Currently only Negative.Control groups are supported.}
  \item{chan} {Fluorescence channel of interest ("FL1-H", "PE-H", etc.)}
  \item{\dots}{optional arguments}
}
\value{
  Returns a list of groups, where each group contains a single negative
  control well and the associated test well for a particular channel. 
}
\author{ Errol Strain }

\examples{
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# Create a set of negative control gates and then apply them
negCon <- getGroups(fp,chan="FL1-H")

negCon[1:2]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

