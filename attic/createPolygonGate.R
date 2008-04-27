`createPolygonGate` <-
function (aFlowFrame, params=c("FSC-A", "SSC-A"),
                               gateName=NA, ...) {

  if (class(aFlowFrame) != "flowFrame") {
    stop("You can only draw gates on a flowFrame, I give up!");
  }

  index      = vector();
  dyeNames   = as.vector(colnames(aFlowFrame@exprs));
  paramNames = as.vector(getParameterNames(aFlowFrame));
  
  # Accept either the name of the fluor, or the antibody target
  for (i in 1:2) {
    value = grep(params[i], dyeNames, ignore.case=TRUE);
    if (length(value) == 0) {
      value = grep(params[i], paramNames, ignore.case=TRUE);
    }
    index[i] = value;
  }
  
  if (is.na(gateName)) {
    gateName = paste("Polygon on", params[1], "and", params[2]);
  }
  
  plot(aFlowFrame@exprs[,index[1]], aFlowFrame@exprs[,index[2]],
        xlab=params[1], ylab=params[2], pch=".",
        main=aFlowFrame@description$'$FIL', ...);

  coords = locator(type="l");

  polyCoords <- matrix(ncol=2, nrow=length(coords$x));
  polyCoords[,1] = coords$x;
  polyCoords[,2] = coords$y;
  colnames(polyCoords) = dyeNames[index];
  myGate = polygonGate(filterId=gateName, boundaries=polyCoords);

  return(myGate);
}

