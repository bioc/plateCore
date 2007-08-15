\name{read.flowPlate}
\alias{read.flowPlate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to read in a directory containing flow plate data. }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
read.flowPlate(path = ".", pattern = NULL, wellAnnotation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{ Relative path to the directory containing fcs files. }
  \item{pattern}{ Location of the 3-letter well identifier in the file name, number of characters from the end.
  	Defaults to c(6,4), indicating the well identifier is between 6 and 4 character from the end of the filename
  	string. }
  \item{wellAnnotation}{ Annotated data frame describing the plate layout. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{ Returns a flowSet with plate phenoData. }
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
function(path=".",pattern=NULL,wellAnnotation) {
	
	files = list.files(path) 

	## Now get the well ids from the sample filenames
	if(is.null(pattern)) {
		Well.Id <- unlist(lapply(files, function(x) substr(x, start=nchar(x)-2,stop=nchar(x))))
	} else {
		Well.Id <- unlist(lapply(files, function(x) substr(x, start=nchar(x)-pattern[1],stop=nchar(x)-pattern[2])))
	}

	## Create a data frame with the sample file names and corresponding well id
	name = as.character(list.files(path))
	tempPlate.df <- data.frame(cbind(name,Well.Id),stringsAsFactors=FALSE)

	## Merge the tempPlate.df with a subset of the current phenoData (excluding name column from phenoData)
	## name gets overwritten with the sample file names
	tempPlate2.df <- merge(tempPlate.df,subset(pData(wellAnnotation),select=-name),by="Well.Id")

	## The name column is usually assumed to be in the first position for the other flowCore stuff
	name = as.character(tempPlate2.df$name)
	tempPlate2.df <- subset(tempPlate2.df,select=-name)
	tempPlate2.df <- data.frame(cbind(name,tempPlate2.df),stringsAsFactors=FALSE)

	## Now replace the phenoData in the wellAnnotation.adf with the new data frame
	row.names(tempPlate2.df) <- name
	pData(wellAnnotation) <- tempPlate2.df
	
	## It's still converting strings to factors for some reason, fixing it here for now.
	pData(wellAnnotation)$name <- as.character(pData(wellAnnotation)$name)
		
	## Call the read.flowSet function to actually read in the plates
	newFlowSet <- read.flowSet(path=path)

	## Copy over the annotated data frame describing the experiment
    phenoData(newFlowSet) <- wellAnnotation
    
    ## Get rid of the dashes
    newFlowSet <- fsApply(newFlowSet,function(x) {
    	newNames <- gsub("-",".",colnames(exprs(x)))
    	colnames(exprs(x)) <- newNames
    	x
    })
    
    gc()
    
	return(newFlowSet)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
