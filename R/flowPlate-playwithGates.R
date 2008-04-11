###############################################################################
#
# Filename: flowPlate-playwithGates.R
#
# Author: straine
#
# Description: Using the playwith package to have some rudimentary way of adjusting
#	the negative control gates visually.  Unfortunately, I can't return just a numeric
#   value, I have to return the environment (the function returns with the playwith
#    graphics window is still open, and the modal=TRUE option to playwith doesn't seem
#    to actually stop the console.)
#
# Date : Apr 6, 2008
###############################################################################


setGate <- function(fp,chan,wells,newGate,type="Negative.Control.Gate") {
	fp@wellAnnotation[plate@wellAnnotation$Channel==chan & plate@wellAnnotation$Well.Id %in% wells,type] <- newGate
	return(fp)
}

adjustGateLog10 <- function(fp,wells,chan,xlim,ylim,fsc="FSC-A",asFact="Well.Id",numEvents=250,type="Negative.Control.Gate") {
	
	require(playwith)
	
	dataDf <- fsApply(plateSet(fp[wells]),function(x) {
				data.frame(exprs(x)[sample(1:nrow(exprs(x)),min(nrow(exprs(x)),numEvents)),c(fsc,chan)])		
			})
	
	groupDf <- cbind(dataDf[[1]],Well.Id=wells[1])
	
	for(i in 2:length(wells)) {
		groupDf <- rbind(groupDf,cbind(dataDf[[i]],Well.Id=wells[i]))
	}
	
	NegGate <- median(subset(fp@wellAnnotation,Channel==chan & Well.Id %in% wells,select=type)[[1]])
	
	gateTool <- function(playState) {
		gateScrollBox <- gtkHBox(show=TRUE)
		vbox <- gtkVBox()
		vbox$packStart(gateScrollBox, expand=FALSE)
		gateScrollBox$packStart(gtkLabel("Gate (Linear) :"), expand=FALSE)
		gateEntry <- gtkEntry()
		gtkEntrySetText(gateEntry,as.character(playState$env$NegGate))
		gateEntry["width-chars"] <- 30
		gSignalConnect(gateEntry, "activate", 
				gate_handler, data=playState)
		gateScrollBox$packStart(gateEntry, expand=FALSE)
		foo <- gtkToolItem()
		foo$add(vbox)
		foo
	}
	
	## this is called when the button is clicked
	gate_handler <- function(widget, playState) {
		playState$env$NegGate <-  as.numeric(widget["text"])
		playReplot(playState)
	}	
	
	latForm <- as.formula(paste("log10(",gsub("-",".",chan),") ~ ",gsub("-",".",fsc)," | as.factor(",asFact,")",sep=""))
	
	playwith(xyplot(latForm, 
					groupDf,xlim=xlim,ylim=ylim,NegGate=NegGate,
					panel=function(x,y,...){
						panel.xyplot(x,y,pch=19,cex=0.5)
						panel.abline(h=log10(NegGate))
						pp <- round(100*sum(y>log10(NegGate))/length(y),digits=0)
						panel.text(0.9*xlim[2],0.9*ylim[2],paste(pp,"%",sep=""),cex=0.75)
					}),modal=TRUE,bottom=list(gateTool)
	)
	
	return(playDevCur()[['env']])
}


