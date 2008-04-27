`setUpPlateMap` <-
function (aFlowSet, ...) {
  
  retList       = list();
  pd            = phenoData(aFlowSet);
  cols          = as.vector(pd@data$Column);
  rows          = as.vector(pd@data$RowSymbol);
  plotOrder     = order(rows, cols);
  wellNames     = paste(rows, cols, sep="");
  wellMapVec    = mapPresentWells(wells=wellNames, returnType="vector", plateSize=96);
  wellMapMat    = mapPresentWells(wells=wellNames, returnType="matrix", plateSize=96);
  presentWells  = grep(TRUE, wellMapVec);

  retList[["columns"]]    = cols;
  retList[["rows"]]       = rows;
  retList[["plotOrder"]]  = plotOrder;
  retList[["wellNames"]]  = wellNames;
  retList[["wellMapVec"]] = wellMapVec;
  retList[["wellMapMat"]] = wellMapMat;
  retList[["present"]]    = presentWells;
  
  invisible(retList);
}

