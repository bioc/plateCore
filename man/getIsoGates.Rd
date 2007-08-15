\name{getIsoGates}
\alias{getIsoGates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create rectangleGates for each isotype control.  }
\description{
  Isotype controls are used to define the boundary between positive/negative cells.
}
\usage{
getIsoGates(data, cutOff = 0.99, unstainDelta = 0.024)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ flowSet decribed using a plate annotated data frame.  }
  \item{cutOff}{ Quantile for cutoff between positive/negative cells. }
  \item{unstainDelta}{ Constant added to the quantile cutoff. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Returns an environment containing isotype gates, with 
  each gate named according to it's associated fluorochrome. }
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
function(data,cutOff=0.99,unstainDelta=0.024) {

	## Get the names of the columns with dyes
    dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	
	## Make a data frame containing the dye/well info for the isotype controls
	iso.df <- subset(pData(phenoData(data)),as.logical(Sample.Type=="Isotype"),select=c("Well.Number",dyeCols))

	## Get the list of dyes, wells, and column names corresponding to the dyes.  Note that this assumes one dye 
	## per well.
	isoDyes <- apply(iso.df,1,function(x) x[dyeCols[x[dyeCols]!=""]] )
	isoWells <- names(isoDyes)
	isoCols <- apply(iso.df,1,function(x) gsub("\\.","-",sub("\\.dye","",dyeCols[x[dyeCols]!=""]) ))
	
	## Create a hash for storing the rectangle gates.
	isoHash <- new.env(hash=TRUE,parent=emptyenv())

	## For each Isotype control, create a gate based on the channel of interest and foward scatter
	lapply(1:length(isoCols),function(x) {
		thresh <- quantile(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]]),probs=c(cutOff)) + unstainDelta 
		rangeF <- ((range(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]])))[1])	
		if(rangeF>0) {rangeF <- 0}
		names(thresh) <- isoCols[x]			
		names(rangeF) <- isoCols[x]
		assign(isoDyes[[x]],new("rectangleGate",filterId="rectangleGate",parameters=isoCols[x],min=rangeF,max=thresh), env=isoHash)
	})

	return(isoHash)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
