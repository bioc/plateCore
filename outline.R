###################################################
### chunk number 1:  eval=FALSE
###################################################
## platePBMC <- Subset(platePBMCraw, rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400)))	


###################################################
### chunk number 2: morphGate
###################################################
library(plateCore)
data(plateCore)
pbmcFP <- flowPlate(pbmcPlate,wellAnnotation,plateName="PBMC.001")

#densityplot(~ `FL1-H` | as.factor(Well.Id),pbmcFP[1:4])

#densityplot(~ `FL1-H` | as.factor(Well.Id),transform("FL1-H"=log10) %on% pbmcFP[1:4])

print(xyplot(`SSC-H` ~ `FSC-H` | as.factor(name), pbmcFP[1], smooth=TRUE,
		filter=rectangleGate("FSC-H"=c(300,600),"SSC-H"=c(50,400)) ))


###################################################
### chunk number 3: fluidic
###################################################
pbmcFP2 <- Subset(pbmcFP, rectangleGate("FSC-H"=c(300,700),"SSC-H"=c(50,400)))	

col2 = c("darkseagreen","red","royalblue","brown","orange","turquoise", "orchid","yellow","black","green","darkred","lightblue")
groups2 = c(1:12)
print(ecdfplot(~`FSC-H`|as.factor(Row.Id),plateSet(pbmcFP),col=col2,lty=c(rep(1,6),rep(2,6)),
				key = list(lines=list(col=col2), lty=c(rep(1,6),rep(2,6)), cex=1, text =list(as.character(groups2)), columns = 4, title="")))	



###################################################
### chunk number 4:  eval=FALSE
###################################################
## platePBMC <- setControlGates(platePBMC,gateType="Negative.Control",numMads=6)


###################################################
### chunk number 5:  eval=FALSE
###################################################
## platePBMC <- applyControlGates(platePBMC)
## platePBMC <- summaryStats(platePBMC)


###################################################
### chunk number 6: 
###################################################
pbmcFP2 <- setControlGates(pbmcFP2,gateType="Negative.Control",numMads=5)
pbmcFP2 <- applyControlGates(pbmcFP2)
platePBMC <- summaryStats(pbmcFP2)


###################################################
### chunk number 7: isoGate
###################################################
wells <- unique(subset(platePBMC@wellAnnotation,Negative.Control=="A03")$Well.Id)
densityplot(~ `FL1-H` | as.factor(Well.Id),transform("FL1-H"=log10) %on% pbmcFP2[wells],filterResult="Negative.Control")
xyplot(`FL1-H` ~ `FSC-H` | as.factor(Well.Id),transform("FL1-H"=log10) %on% plateSet(pbmcFP2[wells]),filterResult="Negative.Control")


print(xyplot(`SSC-H` ~ `FSC-H` | as.factor(name), plateSet(platePBMC)[1:4], smooth=TRUE,
				filter=rectangleGate("FSC-H"=c(300,600),"SSC-H"=c(50,400))))

###################################################
### chunk number 8: wellAnnotation
###################################################
head(wellAnnotation(platePBMC))


###################################################
### chunk number 9:  eval=FALSE
###################################################
## densityplot(~ `FL2-H` | as.factor(plateName), virtPlate,filterResult="Negative.Control")


###################################################
### chunk number 10: pbmcHeat
###################################################
library(gplots)
temp <- read.delim("ppExp.csv",stringsAsFactors=FALSE,header=TRUE)
rownames(temp) <- temp[,1]
temp <- temp[,-1]

mapcol <- rev(heat.colors(10))
b <- seq(0,100,by=10)
#heatmap.2(as.matrix(temp),col=mapcol,breaks=b,scale='none',trace='none',cexRow=0.3)
heatmap.2(as.matrix(temp),col=redgreen(10),breaks=b,scale='none',trace='none',cexRow=0.3)


###################################################
### chunk number 11: pbmcCDbd69
###################################################

fileNames <- list.files("../../publicationPlateCore/pbmcRData",full.names=TRUE)

plates <- lapply(fileNames,function(x){
			load(x)
			platePBMC[c("A03","B06")]
		})

virtPlate <- fpbind(plates[[1]],plates[[2]],plates[[3]],plates[[4]],plates[[5]])


print(densityplot(~ `FL2-H` | as.factor(plateName), transform("FL2-H"=log10) %on% virtPlate,filterResult="Negative.Control"))



