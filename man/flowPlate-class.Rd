\name{flowPlate-class}
\docType{class}
\alias{flowPlate-class}
\alias{[,flowPlate,ANY-method}
\alias{[[,flowPlate,ANY-method}


\title{The flowPlate class.}
\description{
	\code{flowPlates} are the 
}

\section{Slots}{
  
  \describe{
    
    \item{\code{plateName}:}{A \code{character} string containing the name of the plate.
      }

    \item{\code{plateSet}:}{A \code{\link[flowCore::flowSet-class]{flowSet-class}} containing FCS
	  event data. Prior to creating a \code{flowPlate}, the FCS files are first read
	  into a \code{\link[flowCore::flowSet-class]{flowSet-class}} using \code{\link[flowCore::read.flowSet]{read.flowSet}}.
      }

    \item{\code{wellAnnotation}:}{A \code{data.frame} describing the layout of the
	  plate. Each row describes one channel for a well.
	  }
    
  }
}



\section{Methods}{
  
  \describe{

    \item{[, [[}{Subsetting. \code{x[i]} where \code{i} is either a scalar or
      \code{character} corresponding to a sample name,
      returns a \code{flowPlate} object, and \code{x[[i]]} a
      \code{\linkS4class{flowFrame}} object. 
      
      \emph{Usage:}

      \code{   flowSet[i]}
      
      \code{   flowSet[[i]]}
      
    }

  }
}


\author{ Errol Strain }
\examples{
library(plateCore)
data(plateCore)

##Look at the wellAnnotation
wellAnnotation[1:4,]

## Get the lymphocytes
rectGate <- rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400))
pbmcPlate <- Subset(pbmcPlate, rectGate)

## Create a flowPlate object from the platePBMC and the wellAnnotation
fp <- flowPlate(pbmcPlate,wellAnnotation,plateName="P1")

## Subset the flowPlate, creating another flowPlate
fpSmall <- fp["A01"]

## Extract a flowFrame from a flowPlate
ff <-fp[["A01"]]

}
\keyword{ methods }
