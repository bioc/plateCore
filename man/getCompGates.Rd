\name{getCompGates}
\alias{getCompGates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Retrieve rectangleGates of the unstained cells in compensation samples. }
\description{
  Create rectangleGates for removing the negative cells from compensation samples
  containing a mixture of postive and negative cells.  
}
\usage{
getCompGates(data, dyeCols, cutOff = 0.99, fsc = "FSC.A")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ flowSet of compensation data. }
  \item{dyeCols}{ Channels of interest for compensastion }
  \item{cutOff}{ Quantile cutoff for positive/negative cells }
  \item{fsc}{ Forward scatter channel (FSC). }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{rectangleGates in channel of interest plus FSC.}
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
function(data,dyeCols,cutOff=0.99,fsc="FSC.A") {
	
	## Create a hash for storing the rectangle gates.
	compHash <- new.env(hash=TRUE,parent=emptyenv())

	lapply(1:length(dyeCols),function(x) {
		thresh <- quantile(unlist(data@exprs[,dyeCols[x]]),probs=c(cutOff))
		rangeF <- range(unlist(data@exprs[,dyeCols[x]]))
		assign(dyeCols[[x]],!new("rectangleGate",filterId="rectangleGate",parameters=dyeCols[x],min=rangeF[1],max=thresh), env=compHash)
	})
	
	rangeF <- range(unlist(data@exprs[,fsc]))
	assign("unstain",new("rectangleGate",filterId="rectangleGate",parameters=fsc,min=rangeF[1],max=rangeF[2]), env=compHash)

	return(compHash)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
