\name{getLinIsoGates}
\alias{getLinIsoGates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create rectangleGates for each isotype control. }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getLinIsoGates(data, numMads = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ flowSet decribed using a plate annotated data frame. }
  \item{numMads}{ Number of median absolute deviations above the median used to set the 
  boundary between positive/negative cells. }
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
function(data,numMads=3) {

	## Get the names of the columns with dyes
    dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]
	
	## Make a data frame containing the dye/well info for the isotype controls
	iso.df <- subset(pData(phenoData(data)),as.logical(Sample.Type=="Isotype"),select=c("Well.Number",dyeCols))

	## Get the list of dyes, wells, and column names corresponding to the dyes.  Note that this assumes one dye 
	## per well.
	isoDyes <- apply(iso.df,1,function(x) x[dyeCols[x[dyeCols]!=""]] )
	isoWells <- names(isoDyes)
	isoCols <- apply(iso.df,1,function(x) sub("\\.dye","",dyeCols[x[dyeCols]!=""]))
	
	## Create a hash for storing the rectangle gates.
	isoHash <- new.env(hash=TRUE,parent=emptyenv())

	## For each Isotype control, create a gate based on the channel of interest and foward scatter
	lapply(1:length(isoCols),function(x) {
		mfi <- median(10^(data[[isoWells[x]]]@exprs[,isoCols[x]]))
		mfi.mad <- mad(10^(data[[isoWells[x]]]@exprs[,isoCols[x]]))
		thresh <- c(1,log10(mfi+numMads*mfi.mad)) 
		rangeF <- (range(unlist(data[[isoWells[x]]]@exprs[,isoCols[x]])))[1]	
		if(rangeF>0) {rangeF <- 0}
		rangeF <- c(0,rangeF)
		names(thresh) <- c("FSC.A",isoCols[x])			
		names(rangeF) <- c("FSC.A",isoCols[x])
		assign(isoDyes[[x]],new("rectangleGate",filterId="rectangleGate",parameters=c("FSC.A",isoCols[x]),min=rangeF,max=thresh), env=isoHash)
	})

	return(isoHash)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
