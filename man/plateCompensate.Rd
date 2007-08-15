\name{plateCompensate}
\alias{plateCompensate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to compensate flow plate data. }
\description{
  Flow plate data may contain wells with different combinations
  of fluorochromes.  This function identifies wells containing
  two or more fluorochromes and corrects for the spillover.
}
\usage{
plateCompensate(plateSet, comp.mat, chanCols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plateSet}{ A flowSet described using plate phenoData. }
  \item{comp.mat}{ Compensation matrix where the row/column names are the fluoresence channels. }
  \item{chanCols}{ Character vector of fluorescence channels. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{ Returns a compensated flowSet.}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ Errol Strain }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(plateSet,comp.mat,chanCols) {

	dyeCols <- unlist(lapply(chanCols,function(x) {paste(x,".dye",sep="")}))
	
	plateSet <- fsApply(plateSet,function(x) {
		fileName <- attributes(x)$descriptio[["$FIL"]]
		dyeCols <- pData(phenoData(plateSet))[fileName,dyeCols]!=""
		if(sum(dyeCols)>=2) {
			x <- compensate(x,comp.mat[chanCols[dyeCols],chanCols[dyeCols]])
			
		} else {
			x
		}
	})	
	
	return(plateSet)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
