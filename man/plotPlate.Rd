\name{plotPlate}

\title{ plotPlate }
\description{
  Make a row vs. column plot of a plate, where the wells are colored according to some value 
  of choice (e.g. number of events, median signal intensity, percent positive, etc.)
}
\usage{
plotPlate(fp, x = NA, method = "median", main, col, values, width = 1, na.action = "zero", ...)
}

\arguments{
  \item{fp}{ A \code{flowPlate}. }
  \item{x}{ A \code{character} indicating the variable of interest.  Valid choices are 
  "events", any single channel name (e.g. FSC.H, SSC.H, FL1.H, etc.), or vector of channel names
   if the method is mahalanobis. }
  \item{method}{ Valid choices are mean, median, sd, mad, or mahalanobis. }
  \item{main}{ Main text for  the plot }
  \item{col}{ Character vector of two colors. }
  \item{values}{ Optional list of values, with names corresponding to sampleNames, that will be 
  used for plotting. }
  \item{width}{ Width of the well. }
  \item{na.action}{ Handling of NA values, either "zero" or "omit". }

}

\value{
  Plots the plate to the standard output.
}
\references{ The original version of this plot came from the \code{prada} pacakge. }
\author{ Jon Gosink and Errol Strain }

\examples{
library(plateCore)
data(plateCore)

## Create a flowPlate
fp <- flowPlate(pbmcPlate,wellAnnotation,"p1001")

plotPlate(fp,x="events",method="median",col=c("yellow", "darkblue"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
