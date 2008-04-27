`mahalanobisPlateHistogram` <-
function (aFlowSet, ...) {

  mahalList = flowSetMahalanobis(aFlowSet, ...);
  pd        = phenoData(aFlowSet);
  plateName = as.vector(unique(pd@data$PlateBarcode));

  plot(density(mahalList), xlab="Mahalanobis Distance",
       ylab="Number of Wells", main=plateName);
}

