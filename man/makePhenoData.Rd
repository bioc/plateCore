\name{makePhenoData}
\alias{makePhenoData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create an annotated data frame for a flowSet. }
\description{
  This function takes a data frame describing the plate layout, and returns an annotated data frame.
}
\usage{
makePhenoData(plateDesc.df, tempMeta.df = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plateDesc.df}{ Data frame describing the plate laytout. }
  \item{tempMeta.df}{ Optional meta data describing variables in plateDesc.Df. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Annotated data frame describing the layout of a flowSet. }
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
function(plateDesc.df,tempMeta.df=NULL) {

	if(is.null(tempMeta.df)) {
		tempMeta.df <- data.frame(colnames(plateDesc.df))
		rownames(tempMeta.df) <- colnames(plateDesc.df)
	}
	
	tempPheno.adf <- new("AnnotatedDataFrame", data=plateDesc.df, varMetadata=tempMeta.df, dimLabels=c("rowNames", "colNames"))
    sampleNames(tempPheno.adf) <- plateDesc.df$name
	return(tempPheno.adf)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
