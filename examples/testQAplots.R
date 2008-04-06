# TODO: Add comment
#
# Author: straine
###############################################################################

library(plateCore)
data(plateCore)

fpRaw <- flowPlate(pbmcPlate,wellAnnotation,"p1001")

rectGate <- rectangleGate("FSC.H"=c(300,700),"SSC.H"=c(50,400))
normGate <- norm2Filter("SSC.H","FSC.H",scale.factor=1.5)
fp <- Subset(fpRaw, rectGate & normGate)

fp <- fixAutoFl(fp,fsc="FSC.H",chanCols=c("FL1.H","FL2.H","FL3.H","FL4.H"))

fp <- setControlGates(fp,gateType="Negative.Control",numMads=5)
fp <- applyControlGates(fp)

fp <- summaryStats(fp)



win.metafile("plotPlate_TotalEvents.emf")
plotPlate(fpRaw,x="events",method="median",col=c("yellow", "darkblue"))
dev.off()

win.metafile("plotPlate_FL1median.emf")
plotPlate(fpRaw,x="FL1.H",method="median",col=c("yellow", "darkblue"))
dev.off()

win.metafile("plotPlate_FL1mad.emf")
plotPlate(fpRaw,x="FL1.H",method="mad",col=c("yellow", "darkblue"))
dev.off()

win.metafile("plotPlate_Mahalanobis.emf")
plotPlate(fp,x=c("FSC.H","SSC.H"),method="mahalanobis",col=c("yellow", "darkblue"),
		main=paste("Mahalanobis Distance For:", fp@plateName))
dev.off()

win.metafile("medianFSCvsSSC.emf")
data <- fsApply(plateSet(fp),function(x) {apply(exprs(x)[,c("FSC.H","SSC.H")],2,median)})
data <- cbind(Well.Id=pData(phenoData(fp)[,"Well.Id"]),data)

plot(data[,2],data[,3],type="n",xlab="Median FSC",ylab="Median SSC")
text(data[,2],data[,3],data[,1])
dev.off()

win.metafile("gutterPlot.emf")
gutterPlot(fp,chans=c("FSC.H","SSC.H","FL1.H","FL2.H","FL3.H","FL4.H"))
dev.off()

jpeg("timePlot.jpeg")
xyplot(FSC.H ~ Time | as.factor(Well.Id),fpRaw,layout=c(12,8),
		strip=strip.custom(par.strip.text=list(cex=0.7)))
dev.off()

win.metafile("plotPlate_PositiveCount.emf")
plotPlate(fp,x="FL2.H",method="Positive.Count",col=c("yellow", "darkblue"))
dev.off()

win.metafile("plotPlate_PercentPositive.emf")
plotPlate(fp,x="FL2.H",method="Percent.Positive",col=c("yellow", "darkblue"))
dev.off()