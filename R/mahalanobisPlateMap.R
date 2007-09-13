`mahalanobisPlateMap` <-
function (aFlowSet, ...) {

  mahalList = flowSetMahalanobis(aFlowSet, ...);
  pd        = phenoData(aFlowSet);
  plateName = as.vector(unique(pd@data$PlateBarcode));
  plateMap  = setUpPlateMap(aFlowSet);

  # Note, I log the values here..  Log scale is more sensitive to minor problems.
  # Linear scale just shows the major problems.
  plotPlate(mahalList[plateMap$plotOrder], ind=plateMap$present,
            main=paste("Mahalanobis Distance For:", plateName),
            col=c("yellow", "darkblue"), cex.main=1, cex.desc=1,
            desc=c("Weird", "Normal"));
}

