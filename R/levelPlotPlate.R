`levelPlotPlate` <-
function (aFlowSet) {

  colorVals   <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404");  # Yellow -> Brown
  colorScheme <- colorRampPalette(colorVals)(20);
  pd          = phenoData(aFlowSet);
  plateName   = as.vector(unique(pd@data$PlateBarcode));

  levelplot(`SSC-A` ~ `FSC-A` | Column + RowSymbol, aFlowSet,
            n=15, exclude.time=TRUE, labels=FALSE, contour=FALSE,
            col.regions=colorScheme, colorkey=FALSE, strip=FALSE,
            main=plateName);
}

