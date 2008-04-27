`summarizePlate.FCS` <-
function (aFlowSet, channelNames="target", discardBlankCols=TRUE,
                                descriptiveRowNames=FALSE, method="median") {

  pd           = phenoData(aFlowSet);
  factors      = pd@data;
  fileNames    = factors$name;
  mat          = matrix();
  descVec      = vector();

  # CREATE A MATRIX WITH A VALUE FOR EACH STAIN PARAMETER FOR EACH WELL
  # FOR EACH PLATE.  RIGHT NOW THAT IS EITHER BY MEDIAN OR MODE INTENSITY.
  if (method == "median") {
    mat = fsApply(aFlowSet, each_col, median);

  } else if (method == "mode") {
    mat = fsApply(aFlowSet, mode.FCS);
  }

  # I PREFER TO USE THE NAME OF THE TARGET MOLECULE, NOT THE FLUOR.
  if (channelNames == "target") {
    dimnames(mat)[[2]] = getParameterNames(aFlowSet);
  }

  # THEN TACK THE SUMMARIZED DATA AND THE FACTORS TOGETHER.
  for (i in 1:dim(factors)[2]) {
    if (class(factors[,i]) == "character") {
      factors[,i] = as.numeric(as.factor(factors[,i]));
    } else {
      factors[,i] = as.numeric(factors[,i]);
    }
  }

  summaryArray = cbind(mat, factors);



  # ACCUMULATE A SHORT DESCRIPTION OF EACH flowFrame (WHICH MAY OR MAY NOT
  # BE APPLIED AS THE ROWNAMES FOR THE OUTGOING DATA.
  if (descriptiveRowNames == TRUE) {
    for (fileName in fileNames) {
      temp = pd@data[fileName,];
      description = paste(temp$RowSymbol, temp$Column, sep="");
      description = paste(description, temp$SampleBarcode,
                        temp$Treatment1Name, temp$Treatment1Concentration,
                        sep=":");

      descVec   = append(descVec, description);
    }
    dimnames(summaryArray)[[1]] = descVec;
  }


  # Discard any columns where there are only NA's...
  if (discardBlankCols == TRUE) {
    keepColumns = as.vector(rep(TRUE, dim(summaryArray)[2]));
    numWells    = dim(summaryArray)[1];
    for (i in 1:length(keepColumns)) {
      index = is.na(summaryArray[,i]);
      if (sum(index) == numWells) {
        keepColumns[i] = FALSE;
        next;

      # ...and convert NA's to the minimum value for that column
      } else if (sum(index) > 0) {
        replacement = min(summaryArray[,i], na.rm=TRUE);
        summaryArray[index,i] = replacement;
      }
    }
    summaryArray = summaryArray[,keepColumns];
  }

  return(as.matrix(summaryArray));
}

