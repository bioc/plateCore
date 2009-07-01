\name{panel.densityplot.flowPlate}
\alias{panel.densityplot.flowPlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lattice-flowViz style panel function for flowPlate densityplot. }
\description{
  This function should not be called directly, use \code{\link{densityplot}}.
}
\usage{
panel.densityplot.flowPlate(x, 
				frames, channel, wellAnnotation,
				groups=NULL,
				subscripts,
				col = superpose.symbol$col,
				col.points = col,
				col.line = col,
				filterResult=NULL,
				\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ character }
  \item{frames}{ flowFrames }
  \item{channel}{}
  \item{wellAnnotation}{wellAnnotation data.frame}
  \item{groups}{} 
  \item{subscripts}{}
  \item{col}{}
  \item{col.points}{}
  \item{col.line}{}
  \item{filterResult}{}
  \item{\dots}{optional arguments}
}

\author{ Errol Strain }

\seealso{ See Also \code{\link{densityplot}} }
\examples{
# Load the plateCore package and data
library(plateCore)
data(plateCore)

# Create a flowPlate from the sample data in plateCore
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

# Overlay the first 3 flowFrames.  If the groups argument was
# omitted, then the flowFrames would be combined into a single
# density curve.
densityplot(~ `FSC-H`, fp[1:3], groups=name, auto.key=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

