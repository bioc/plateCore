\name{replaceColHeadings}
\alias{replaceColHeadings}

\title{Replace column headings in the exprs() matrix of a flowSet.}
\description{
	The 'replaceColHeadings' function is used substitute the channel names (columns) in the 
	flowFrame elements of a flowSet.  
}
\usage{
replaceColHeadings(flowSet,original=c("FL1-H","FL2-H","FL3-H","FL4-H"),replace=c("FITC.A","PE.A","PerCP.Cy5.5.A","APC.A"))
}

\arguments{
  \item{flowSet}{ A flowSet object from flowCore. }
  \item{original}{ Original channel names. }
  \item{replace}{ Replacement channel names. }

}

\value{
	Returns a flowSet with the new channel names.
}

\author{ Errol Strain }
\examples{
replaceColHeadings(flowSet,"FL1-H","FITC.A")
}

\keyword{methods}
