\name{plateCore-package}
\alias{plateCore-package}
\alias{plateCore}
\docType{package}
\title{ plateCore: A Bioconductor package for high throughput analysis of flow cytometry data
}
\description{
\code{plateCore} is a Bioconductor packaged created to make processing and analysis of large, complex flow datasets
in R easier. High throughput flow studies are often run in a 96 or 384-well plate format, with a number of different samples, 
controls, and antibodies-dye conjugates present on the plate. Analyzing the output from the cytometer requires keeping track of the contents
of each well, matching sample wells with control wells, gating each well/channel separately, making the appropriate plots, and
summarizing the results. \code{plateCore} extends the \code{flowCore} and \code{flowViz} packages to work on
\code{flowPlate} objects that represent these large flow datasets. For those familiar with \Rpackage{flowCore} and \code{flowViz}, 
the gating (filtering), transformation, and other data manipulations for \code{flowPlates} are very similar to \code{flowSets}. 
}
\details{
\tabular{ll}{
Package: \tab plateCore\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-06-29\cr
}

}
\author{ Errol Strain, Florian Hahne, and Perry Haaland
Maintainer: Errol Strain <estrain@gmail.com>
}
\references{
Insert flowCore and flowViz publications.
}

\keyword{ package }

\examples{
library(plateCore)
data(plateCore)

## Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
normGate <- norm2Filter("SSC-H","FSC-H",scale.factor=1.5)
pbmcPlate <- Subset(pbmcPlate, rectGate & normGate)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

}
