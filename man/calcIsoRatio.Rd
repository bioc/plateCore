\name{calcIsoRatio}
\alias{calcIsoRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A function to calculate the MFI ratio of test wells to control wells. }
\description{
  
}
\usage{
calcIsoRatio(data, summ.df, channels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A flowSet with plate phenoData. }
  \item{summ.df}{ Summary data frame from the plateBasicSumm function.}
  \item{channels}{ Character vector of fluorescence channels. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Returns a data frame, with the original summ.df bound to the MFI ratios.}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(data,summ.df,channels) {
	
	mfiCols <- unlist(lapply(channels,function(x) paste(x,".ratioMFI",sep="")))
	posCols <- unlist(lapply(channels,function(x) paste(x,".MFIratio.pos",sep="")))
	negCols <- unlist(lapply(channels,function(x) paste(x,".MFIratio.neg",sep="")))
	
	temp.df <- data.frame(matrix(0,nrow=nrow(summ.df),ncol=length(mfiCols)*3))
	colnames(temp.df) <- c(mfiCols,posCols,negCols)

	## Get the names of the columns with dyes
    dyeCols <- colnames(pData(phenoData(data)))[grep(".*\\.dye",colnames(pData(phenoData(data))))]	


	## Make a data frame containing the dye/well info for the isotype controls
	iso.df <- subset(pData(phenoData(data)),as.logical(Sample.Type=="Isotype"),select=c("Well.Number",dyeCols))

	## Get the list of dyes, wells, and column names corresponding to the dyes.  Note that this assumes one dye 
	## per well.
	isoDyes <- apply(iso.df,1,function(x) x[dyeCols[x[dyeCols]!=""]] )
	isoWells <- as.numeric(iso.df[,1])
	isoCols <- apply(iso.df,1,function(x) sub("dye","MFI",dyeCols[x[dyeCols]!=""]) )
	
	## Create a hash for storing the rectangle gates.
	isoHash <- new.env(hash=TRUE,parent=emptyenv())

	lapply(1:length(isoCols),function(x) {
		assign(isoDyes[[x]],summ.df[isoWells[x],isoCols[[x]]],env=isoHash)
	})


	iso.list <- apply(pData(phenoData(data))[,dyeCols],1, function(x) {
			lapply(x,function(y) {
				if(y!="") { 
					get(y,env=isoHash) 
				}
				else{ NA } 
			})
		})
	
	iso.df <- data.frame(matrix(unlist(iso.list),ncol=length(channels),byrow=TRUE))	
	
	temp.df[,mfiCols] <- round(summ.df[,unlist(lapply(channels,function(x) gsub("-","\\.",paste(x,".MFI",sep=""))))]/iso.df,digits=3)
	temp.df[,posCols] <- round(summ.df[,grep(".*posMFI$",colnames(summ.df))]/iso.df,digits=3)
	temp.df[,negCols] <- round(summ.df[,grep(".*negMFI$",colnames(summ.df))]/iso.df,digits=3)
	
	names(temp.df) <- sub("-","\\.",names(temp.df))
	cbind(summ.df,temp.df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
