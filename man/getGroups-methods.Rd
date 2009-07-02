\name{getGroups}
\alias{getGroups}
\alias{getGroups-method}
\alias{getGroups,flowPlate-method}

\docType{methods}

\title{ Retrieve Negative control groups from a flowPlate }
\description{
   Retrieve a list of negative control-based groups from a flowPlate, based
   on the information in wellAnnotation.
}
\usage{
getGroups(data,type="Negative.Control",chan, \dots)
}

\arguments{
  \item{data}{ A flowPlate dataset. }
  \item{type}{Currently only Negative.Control groups are supported.}
  \item{chan}{Fluorescence channel of interest ("FL1-H", "PE-H", etc.)}
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

\keyword{ methods }

