\name{summaryStats}
\alias{summaryStats}

\title{ Compute summary statistics on a flowPlate }
\description{
  This function computes the median fluorescence intensity (MFI) and the
  MFI ratio (ratio of test well MFI to negative control MFI) for each
  well/channel in a \code{flowPlate}. These results are stored in the
  \code{wellAnnotation data.fame}.
}
\usage{
summaryStats(data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A \code{flowPlate} }
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

# Create a set of negative control gates and then apply them
fp <- setControlGates(fp,gateType="Negative.Control")
fp <- applyControlGates(fp,gateType="Negative.Control")

# Compute summary statistics
fp <- summaryStats(fp)

# There should now be MFI and MFI.ratio columns in the wellAnnotation
head(wellAnnotation(fp))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
