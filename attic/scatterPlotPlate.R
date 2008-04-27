`scatterPlotPlate` <-
function (aFlowSet, parameters=c("FSC-A","SSC-A"),
                              identifier="wellId", method="median") {
  pd          = phenoData(aFlowSet);
  plateName   = as.vector(unique(pd@data$PlateBarcode));
  pointLabels = vector();
  mat         = matrix();
  
  # HERE I GIVE TWO CHOICES, WHETHER YOU WANT TO LABLE THE GRAPH
  # BY THE wellId OR BY THE sampleName OF THE FCS FILE IN THE WELL.
  if (identifier == "wellId") {
    cols        = as.vector(pd@data$Column);
    rows        = as.vector(pd@data$RowSymbol);
    pointLabels = paste(rows, cols, sep="");

  } else if (identifier == "sampleName") {
    pointLabels = as.factor(phenoData(aFlowSet)@data$sampleNames);
    plateName   = paste(plateName[1], "...", plateName[length(plateName)]);
  }
  
  # CREATE A MATRIX WITH A VALUE FOR EACH STAIN PARAMETER FOR EACH WELL
  # FOR EACH PLATE.  RIGHT NOW THAT IS EITHER BY MEDIAN OR MODE INTENSITY.
  mat = summarizePlate.FCS(aFlowSet, discardBlankCols=TRUE,
                           descriptiveRowNames=FALSE, method=method);
  
  rowNames = row.names(mat);
  xLim     = range(mat[,parameters[1]]);
  yLim     = range(mat[,parameters[2]]);
  
  # Set up an empty plot
  plot(mean(xLim), mean(yLim), xlim=xLim, ylim=yLim, type="n",
       xlab=parameters[1], ylab=parameters[2],
       main=paste(method, "Fluorescent Intensity on ", plateName));

  # Then fill it with points, one per each well
  for (i in 1:dim(mat)[1]) {
    points(mat[i,parameters[1]], mat[i,parameters[2]], type="p", pch=".");
    
    # Note, I go through coniptions to get the well ID from the phenoData
    # and not just directly from the plate information.  Is this more reliable?
    # aFullPlate[["20070221c_H1_H01.fcs"]]@description$'WELL ID'
    index = pd@data$sampleNames == rowNames[i];   # roundabout way to do this
    label = pointLabels[index];
    text(mat[i,parameters[1]], mat[i,parameters[2]], labels=c(label));
  }
}

