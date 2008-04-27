`readGatedFCS` <-
function(filePath, compMatrix, filterObj) {

  aFlowFrame = read.FCS(filePath);
  if (!is.null(compMatrix)) {
    aFlowFrame = compensate(aFlowFrame, compMatrix);
  }
  aFlowFrame = Subset(aFlowFrame, filterObj);
  return(aFlowFrame);
}

