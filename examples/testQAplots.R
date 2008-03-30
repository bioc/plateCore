# TODO: Add comment
#
# Author: straine
###############################################################################

library(plateCore)
data(plateCore)

fp <- flowPlate(pbmcPlate,wellAnnotation,"p1001")

plotPlate(fp,x="events",method="median",col=c("yellow", "darkblue"))

plotPlate(fp,x="FL1.H",method="median",col=c("yellow", "darkblue"))

plotPlate(fp,x="FL1.H",method="mad",col=c("yellow", "darkblue"))

#values=c(rep(1,90),rep(10,6))
#names(values) <- sampleNames(fp)
#
#plotPlate(fp,col=c("yellow", "darkblue"),values=values)


plotPlate(fp,x=c("FSC.H","SSC.H"),method="mahalanobis",col=c("yellow", "darkblue"),main=paste("Mahalanobis Distance For:", fp@plateName))



data <- fsApply(plateSet(fp),function(x) {apply(exprs(x)[,c("FSC.H","SSC.H")],2,median)})
data <- cbind(Well.Id=pData(phenoData(fp)[,"Well.Id"]),data)

plot(data[,2],data[,3],type="n",xlab="Median FSC",ylab="Median SSC")
text(data[,2],data[,3],data[,1])


gutterPlot(fp,chans=c("FSC.H","SSC.H","FL1.H","FL2.H","FL3.H","FL4.H"))


xyplot(FSC.H ~ Time | as.factor(Well.Id),fp,layout=c(12,8),smooth=FALSE)



rectGate <- rectangleGate("FSC.H"=c(300,700),"SSC.H"=c(50,400))
normGate <- norm2Filter("SSC.H","FSC.H",scale.factor=1.5)
fp.lymph <- Subset(fp, rectGate & normGate)

fp.lymph <- fixAutoFl(fp.lymph,fsc="FSC.H",chanCols=c("FL1.H","FL2.H","FL3.H","FL4.H"))

fp.lymph <- setControlGates(fp.lymph,gateType="Negative.Control",numMads=5)
fp.lymph <- applyControlGates(fp.lymph)

fp.lymph <- summaryStats(fp.lymph)

plotPlate(fp,x="FL2.H",method="Positive.Count",col=c("yellow", "darkblue"),na.action="zero")
plotPlate(fp,x="FL2.H",method="Percent.Positive",col=c("yellow", "darkblue"),na.action="omit")