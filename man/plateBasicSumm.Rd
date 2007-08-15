\name{plateBasicSumm}
\alias{plateBasicSumm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create a summary data frame for a plate flowSet. }
\description{
  Given a plate flowSet, calculate the MFIs, percent postives, MFI ratios, etc. for
  each well. 
}
\usage{
plateBasicSumm(data, posCells, normCells, nonDeb, channels, fsc = "FSC.A", ssc = "SSC.A")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ flowSet described using an annotated data frame. }
  \item{posCells}{ List of filterResults for positive cells.. }
  \item{normCells}{ List of filterReults for normal morphology cells. }
  \item{nonDeb}{ List of filterResutls for non-debris events. }
  \item{channels}{ Channels of interest. }
  \item{fsc}{ Name of forward scatter channel. }
  \item{ssc}{ Name of side scatter channel. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Data frame of MFIs, percent postive, MFI ratios, for each well.}
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
function(data,posCells,normCells,nonDeb,channels,fsc="FSC.A",ssc="SSC.A") {

	## Column headings for the summary data frame.
	colHeadings <- c("Total.Events","NonDebris.Events","NormalMorhp.Events","FSC.A.MFI","SSC.A.MFI")
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".MFI",sep="")))))
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".pp",sep="")))))
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".posMFI",sep="")))))
	colHeadings <- c(colHeadings,unlist(lapply(channels,function(x) gsub("-",'\\.',paste(x,".negMFI",sep="")))))

	## Function list used to calculate values for columns listed above.
	funList <- vector("list",length=7) 
	
	## Get the number of total events
	funList[[1]] <- function(fs,normCells,posCells,nonDeb) {lapply(nonDeb,function(x) length(x@subSet))}
	## Number of non debris events
	funList[[2]] <- function(fs,normCells,posCells,nonDeb) {lapply(nonDeb,function(x) sum(x@subSet))}
	## Number of normal morphology events
	funList[[3]] <- function(fs,normCells,...) {lapply(normCells,function(x) sum(x@subSet))}
	## Median values for fsc, ssc, and other channels
	funList[[4]] <- function(fs,...) {
		lapply(c(fsc,ssc,channels),function(y) {
			fsApply(fs, function(x) { 
				 round(median(x@exprs[,y]),digits=3) 
			})
		})
	}
	## Calculate percent of cells that are positive
	funList[[5]] <- function(fs,normCells,posCells,...) {
		lapply(1:length(channels),function(y) {
			lapply(1:length(fs),function(x) {
				if(class(posCells[[x]][[y]])=="logicalFilterResult") {
					round(sum(posCells[[x]][[y]]@subSet)/sum(normCells[[x]]@subSet),digits=2)*100
				} else { NA }
			})
		})
	}	
	## MFI of positive cells
	funList[[6]] <- function(fs,normCells,posCells,...) {
		lapply(1:length(channels), function(y) {
			lapply(1:length(fs),function(x) {
				if(class(posCells[[x]][[y]])=="logicalFilterResult")  {
					round(median(fs[[x]]@exprs[posCells[[x]][[y]]@subSet,channels[y]]),digits=3)
				} else { NA }
			})
		})
 	}
 	## MFI of negative cells
	funList[[7]] <- function(fs,normCells,posCells,...) {
		lapply(1:length(channels), function(y) {
			lapply(1:length(fs),function(x) {	
				if(class(posCells[[x]][[y]])=="logicalFilterResult")  {
					round(median(fs[[x]]@exprs[!(posCells[[x]][[y]]@subSet),channels[y]]),digits=3)
				} else { NA }
			})
		})
 	}	

	## apply the function list
	plateSumm.df <- data.frame(matrix(unlist(lapply(funList,function(x) x(data,normCells,posCells,nonDeb))),
					nrow=length(data)),stringsAsFactors=FALSE)

	plateSumm.df <- cbind(pData(phenoData(data))$Well.Id,plateSumm.df)
	
	names(plateSumm.df) <- c("Well.Id",colHeadings)
	
	plateSumm.df
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
