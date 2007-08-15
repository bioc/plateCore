\name{applyIsoGates}
\alias{applyIsoGates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply the isotype rectangle gates to a flowSet. }
\description{
  Isotype gates are used to determine the boundary between positive and 
  negative cells.  For each isotype control, a rectangleGate is created
  in the channel of interest and foward scatter (FSC).  The data is then
  filtered using isoGates to get the number of postive cells in each well.
}
\usage{
applyIsoGates(data, isoGates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A flowSet object which is annotated using plate phenoData.}
  \item{isoGates}{ A list or environment of isotype rectangleGates. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{A list of filterResults.}
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
function(data,isoGates) {

	## Get the names of the columns with dyes
    dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	
	apply(pData(phenoData(data)),1,function(x) {
		lapply(dyeCols,function(y) {
			if((x[y]!="") & (x[y] \%in\% ls(isoGates))) {
				tempFilt <- get(x[y],env=isoGates)
				filter(data[[x["name"]]], !tempFilt)
			} 
			else { NA }
		})
	})
		
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
