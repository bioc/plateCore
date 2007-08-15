###############################################################################
## plotWell.R
## 
## given an flowFrame and parent and child filters (optional) this function calls
## a 2x2 layout of basic plot with FSCxSSC, and then FL1, FL2, and FL3 in all
## combinations
##
## Also maps gene names to wells.
##
# Author: haaland, straine
###############################################################################
plotWell = function(x=comp1,child,well,main,parentFilter,wellMaps.adf,facsTOgene.df,chanCols) {
	if(class(x)!="flowFrame"){
		stop("plate must be a list of flowFrame objects.")
	}
	if(!missing(parentFilter)){
		if(class(parentFilter)!="logicalFilterResult") {
			stop("parent must be a logicalFilterResult.")
		}
	}	

	#########################################################
	## This if/else block creates the labels containing gene
	## names and  channels.
	#########################################################
	plotPs = vector("character",length=0)
	labels = vector("character",length=0)
	channels <- pData(wellMaps.adf)[well,]
	labels2 = vector("character",length=0)
	
	chanLab <- names(channels)[grep("\\.A$",names(channels))]
	chanLab2 <- names(channels)[grep("\\.dye$",names(channels))]

	if(channels[chanLab2[1]]!="") { 
		if(channels[chanLab[1]]!="") {
			labels <- c(unlist(labels),channels[chanLab[1]])
			labels2 <- c(unlist(labels2),channels[chanLab2[1]])
		} else {
			labels <- c(unlist(labels),paste(chanCols[1],"-",channels[chanLab[1]],sep=""))
		}
		plotPs <- c(unlist(plotPs),chanCols[1])
	}
	if(channels[chanLab2[2]]!="") { 
		if(channels[chanLab[2]]!="") {
			labels <- c(unlist(labels),channels[chanLab[2]])
			labels2 <- c(unlist(labels2),channels[chanLab2[2]])
		} else {
			labels <- c(unlist(labels),paste(chanCols[2],"-",channels[chanLab[2]],sep=""))
		}
		plotPs <- c(unlist(plotPs),chanCols[2])
	}
	if(channels[chanLab2[3]]!="") { 
		if(channels[chanLab[3]]!="") {
			labels <- c(unlist(labels),channels[chanLab[3]])
			labels2 <- c(unlist(labels2),channels[chanLab2[3]])
		} else {
			labels <- c(unlist(labels),paste(chanCols[3],"-",channels[chanLab[3]],sep=""))
		}
		plotPs <- c(unlist(plotPs),chanCols[3])
	}	
	if(channels[chanLab2[4]]!="") { 			
		if(channels[chanLab[4]]!="") {
			labels <- c(unlist(labels),channels[chanLab[4]])
			labels2 <- c(unlist(labels2),channels[chanLab2[4]])
		} else {
			labels <- c(unlist(labels),paste(chanCols[4],"-",channels[chanLab[4]],sep=""))
		}
		plotPs <- c(unlist(plotPs),chanCols[4])
	}	

	
	if(!missing(facsTOgene.df) && length(labels2)>0) {
		for(i in 1:length(labels)) {
			if(!(labels[i] %in% c(chanLab,NA))) {
				if(facsTOgene.df[facsTOgene.df[,"Antigen"]==labels[i],"Name"]!="") {
					labels[i] <- paste(labels[i],":",labels2[i]," (",facsTOgene.df[facsTOgene.df[,"Antigen"]==labels[i],"Name"][1],")",sep="")
				}
			}
		}
		labels <- unlist(labels)
	}
	
	par(mfrow=c(2,2))

	if(missing(parentFilter)) {
		plot(x,c("FSC-H","SSC-H"),transformation=function(x) x, xlim=c(0,1),ylim=c(0,1),
				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))	
		if(length(labels2)==0) {
			plot(x,c(chanCols[1],chanCols[2]),transformation=function(x) x,xlab=chanCols[1],ylab=chanCols[2],xlim=c(0,1),ylim=c(0,1),
				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
			plot(x,c(chanCols[2],chanCols[3]),transformation=function(x) x,xlab=chanCols[2],ylab=chanCols[3],xlim=c(0,1),ylim=c(0,1),
				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
			plot(x,c(chanCols[3],chanCols[4]),transformation=function(x) x,xlab=chanCols[3],ylab=chanCols[4],xlim=c(0,1),ylim=c(0,1),
				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
		}

		else if(length(labels)>2) { 
			plot(x,c(plotPs[[1]],plotPs[[2]]),transformation=function(x) x,xlab=labels[[1]],ylab=labels[[2]],xlim=c(0,1),ylim=c(0,1),
				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
   			plot(x,c(plotPs[[1]],plotPs[[3]]),transformation=function(x) x,xlab=labels[[1]],ylab=labels[[3]],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
   			plot(x,c(plotPs[[2]],plotPs[[3]]),transformation=function(x) x,xlab=labels[[2]],ylab=labels[[3]],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))	
		}
		else{
			plot(x,c("FSC-H",plotPs[1]),xlab="FSC-H",transformation=function(x) x,ylab=labels[1],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))	
			plot(x,c("FSC-H",plotPs[2]),xlab="FSC-H",transformation=function(x) x,ylab=labels[2],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))	
		}
	}
	else {
		plot(x[parentFilter],c("FSC-H","SSC-H"),transformation=function(x) x,xlim=c(0,1),ylim=c(0,1),colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))		

		if(length(labels2)==0) {
			plot(x[parentFilter],c(chanCols[1],chanCols[2]),transformation=function(x) x,xlab=chanCols[1],ylab=chanCols[2],
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
			plot(x[parentFilter],c(chanCols[2],chanCols[3]),transformation=function(x) x,xlab=chanCols[2],ylab=chanCols[3],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
			plot(x[parentFilter],c(chanCols[3],chanCols[4]),transformation=function(x) x,xlab=chanCols[3],ylab=chanCols[4],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
		}
		else if(length(labels)==1) {
			plot(x[parentFilter],c("SSC-H",plotPs[1]),transformation=function(x) x,xlab="SSC-H",ylab=labels[1],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
		}		
		else if(length(labels)>2) { 
			plot(x[parentFilter],c(plotPs[[1]],plotPs[[2]]),transformation=function(x) x,xlab=labels[1],ylab=labels[2],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
   			plot(x[parentFilter],c(plotPs[[1]],plotPs[[3]]),transformation=function(x) x,xlab=labels[1],ylab=labels[3],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
   			plot(x[parentFilter],c(plotPs[[2]],plotPs[[3]]),transformation=function(x) x,xlab=labels[2],ylab=labels[3],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))	
		}
		else{
			plot(x[parentFilter],c(plotPs[1],plotPs[2]),transformation=function(x) x,xlab=labels[1],ylab=labels[2],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
			plot(x[parentFilter],c("SSC-H",plotPs[2]),transformation=function(x) x,xlab="SSC-H",ylab=labels[1],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
			plot(x[parentFilter],c("SSC-H",plotPs[3]),transformation=function(x) x,xlab="SSC-H",ylab=labels[2],xlim=c(0,1),ylim=c(0,1),
   				colramp = colorRampPalette(c("white",rev(rainbow(14))[-(1:3)],"#FF0000FF")))
		}
	}
	mtext(outer=TRUE,line=-2,text=main,side=3)
	par(mfrow=c(1,1))
	invisible()

	
}
