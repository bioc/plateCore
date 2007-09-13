# Count the number of gutter events in a flowSet

gutterPlot <- function (aFlowSet, ...) {

  setLength = length(aFlowSet);
  numParams = length(getParameterNames(aFlowSet));
  resultMat = matrix(nrow = setLength, ncol = numParams);
  pd        = phenoData(aFlowSet);
  plateName = as.vector(unique(pd@data$PlateBarcode));
  plateMap  = setUpPlateMap(aFlowSet);

  for (aFile in 1:setLength) {
#    index     = plateMap$plotOrder[position];
    fcsObj    = aFlowSet[[aFile]];
    numEvents = dim(fcsObj@exprs)[1];
    goodMat   = matrix(nrow=numEvents, ncol=numParams);
    
    stats      = summary(fcsObj);
    paramMins  = stats["Min.", ];
    paramMaxes = stats["Max.", ];
    for (i in 1:dim(fcsObj@exprs)[2]) {
      vec = ((fcsObj@exprs[, i] > paramMins[i]) & (fcsObj@exprs[,i] < paramMaxes[i]))
      goodMat[, i] = vec;
    }
    
    # But I want to capture the events at the edge
    goodMat = (goodMat == FALSE);
    resultMat[aFile, ] = apply(goodMat, 2, sum)/numEvents;
  }
  
  resultMat     = resultMat[plateMap$plotOrder,];
  orderedLabels = plateMap$wellNames[plateMap$plotOrder];
  
  # Create the main plot window
  plot(1, 0.5, type="n", xlim=c(1,setLength), ylim=c(0,1),
       xaxt="n",
       xlab="Well ID",
       ylab="% Events Pegged Full or Min Scale", 
       main=plateName, ...);
       
  for (i in 1:numParams) {
    points(1:setLength, resultMat[,i], type="b", pch=i, col=i)
  }
  
  axis(1, at=1:setLength, labels=orderedLabels);
  legend(x="topleft", legend=names(aFullPlate[[1]]), cex=1, 
         bty="n", pch=1:numParams, col=1:numParams);

}


