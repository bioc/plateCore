`timePlateMap2` <-
function (aFlowSet, plateSize=96, parameter=1, trim=c(0.05, 0.95)) {

  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  # Create a 8 X 12 layout of plots
  mat = textPlateMap(1:plateSize);
  nf  = layout(mat, respect=TRUE);
  layout.show(nf);
  par(mar=c(0,0,0,0));

  plateMap   = setUpPlateMap(aFlowSet);
  rangeMat   = fsApply(aFlowSet, each_col, range);
  quantRange = fsApply(aFlowSet, each_col, quantile, probs=trim);
#  yLim     = as.numeric(range(rangeMat[,parameter]));

  # Really, just take the middle of the 5th and 95th percentiles
  yLim     = quantile(quantRange[,parameter])[c(2,4)];
  xLim     = as.numeric(range(rangeMat[,"Time"]));

  goodIndex    = 1;
  for (i in 1:length(plateMap$wellMapVec)) {
    if (plateMap$wellMapVec[i] == TRUE) {

      # Note this bit of indirection is necessary as the flowSet wells may
      # not be in the "correct" order or even of the length of a full plate.  
      # However, the wellMapVec will be in correct order.
      j = plateMap$plotOrder[goodIndex];      
      x = as.vector(aFlowSet[[j]]@exprs[, "Time"]);
      y = as.vector(aFlowSet[[j]]@exprs[, parameter]);
      goodIndex = goodIndex + 1;
    } else {
      x = mean(xLim);
      y = mean(yLim);
    }

    plot(x, y, pch=".", xlim=xLim, ylim=yLim,
        xaxt="n", yaxt="n",col="blue");
    # And plot a smoothed line over the data
    loSmooth = lowess(x,y, f=0.66);
    lines(loSmooth, col="black", lwd=2);
  }

  # Finish up the plate figure
  par(def.par);
  pd          = phenoData(aFlowSet);
  plateName   = as.vector(unique(pd@data$PlateBarcode));
  plotTitle   = plateName;
  if (is.numeric(parameter)) {
     plotTitle = paste(plotTitle, "by", names(aFlowSet[[1]])[parameter]);
  }
  
  mtext(plotTitle, side=3);
  mtext("Time", side=1, cex=0.75);
}

