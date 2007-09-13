`parallelPlateMap` <-
function (aFlowSet, numEvents, plateSize=96) {

  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  # Create a 8 X 12 layout of plots
  mat = textPlateMap(1:plateSize);
  nf  = layout(mat, respect=TRUE);
  layout.show(nf);
  par(mar=c(0,0,0,0));

  plateMap = setUpPlateMap(aFlowSet);

  goodIndex   = 1;
  timeChannel = guessTimeChannel(aFlowSet);
  paramNames  = getParameterNames(aFlowSet);
  for (i in 1:length(plateMap$wellMapVec)) {
    plot(1,1, pch=".", xlim=c(1,length(paramNames)),
         ylim=c(0,1), xaxt="n", yaxt="n", type="n");

    if (plateMap$wellMapVec[i] == TRUE) {
      j   = plateMap$plotOrder[goodIndex];
      mat = aFlowSet[[j]]@exprs[1:numEvents, -timeChannel];
      mat = normalizeMatrix(mat);
      for (j in 1:numEvents) {
        lines(mat[j,], col=palette(rainbow(12))[j]);
      }
      goodIndex = goodIndex + 1;
    }
  }

  par(def.par); # reset default paramters
  pd          = phenoData(aFlowSet);
  plateName   = as.vector(unique(pd@data$PlateBarcode)); mtext(plateName, side=3);
  mtext(paste(paramNames[-timeChannel], collapse=" | "), side=1, cex=0.75);
}

