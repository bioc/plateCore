\name{getMarkers}
\alias{getMarkers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get flow plate data output by marker. }
\description{
  Given a flowSet, summary data frame, and the channels of interest retrun a list of 
  markers with the appropriate results (MFIs, percent positive, etc. )
}
\usage{
getMarkers(data, summ.df, chanCols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ flowSet decribed using a plate annotated data frame. }
  \item{summ.df}{ Data frame summarizing the flow plate data. }
  \item{chanCols}{ Fluoresence channels of interest. }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{ Data frame of markers with MFIs, percent postives, MFI ratios, etc.}
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
function(data,summ.df,chanCols) {

	## Markers (genes) are only in the "Test" samples
	genes.df <- subset(pData(phenoData(data))[,c("Well.Id",chanCols,"Sample.Type")],Sample.Type=="Test",select=-Sample.Type)

	## Make a list of data frames showing the well/channel for each marker 	
	genes2.df <- lapply(chanCols,function(x) {
		 temp.df <- cbind(genes.df[,c(x,"Well.Id")],x)
		 colnames(temp.df) <- c("Antigen","Well.Id","Channel")
		 temp.df 
		})
	
	## Now bind the list of dataframes into a single data frame
	temp.df <- data.frame(matrix(nrow=0,ncol=3))
	colnames(temp.df) <- c("Antigen","Well.Id","Channel")
	for(i in 1:length(genes2.df)) {
		temp.df <- rbind(temp.df,genes2.df[[i]])
	}
	geneMap.df <- temp.df[temp.df[,"Antigen"]!="",]
	
	rownames(summ.df) <- as.character(summ.df$Well.Id)
	
	## Now get the marker MFI, etc. values from the plate data
	out.mat <- apply(geneMap.df,1,function(x) {
		matrix( c(x["Antigen"],summ.df[x[["Well.Id"]],c(grep(x[["Channel"]],colnames(summ.df)))]) ,nrow=1 )
	})
	
	cols <- c("Antigen",sub(paste(chanCols[1],"\\.",sep=""),"",colnames(summ.df)[grep(chanCols[1],colnames(summ.df))])) 
	
	out.mat <- matrix(unlist(out.mat,recursive=TRUE),ncol=length(cols),byrow=TRUE)
	colnames(out.mat) <- cols
	data.frame(out.mat,stringsAsFactors=FALSE)
	
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
