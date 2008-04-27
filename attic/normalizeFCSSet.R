`normalizeFCSSet` <-
function (fcsSet, params=list(), method="multiplicative",
                             normalizeTime=FALSE) {

  if (normalizeTime == TRUE) {
      cat("\n\n    OK, I'm normalizing the 'Time' column -- that's pretty weird\n\n");
  }

  rightChanneOrder = getParameterNames(fcsSet[[1]]);
  modeArray        = array(dim=c(length(fcsSet), length(rightChanneOrder)),
                           dimnames=list(NULL, rightChanneOrder));
                    
  # Loop over the set the first time, figuring out what all the medians are
  for (i in 1:length(fcsSet)) {
    thisFile = fcsSet[[i]];

    # They should all have the same list of parameters, in the same order
    if (all.equal(getParameterNames(thisFile), rightChanneOrder) != TRUE) {
      stop(paste("Error!   File", i,
                 "doesn't have the same list of parameters as the first file"));
    }

    # loop over every channel and record the modes
    values        = mode.FCS(thisFile);
    modeArray[i,] = values;
  }

  # Now find the average of all these modes for each parameter
  normVector = as.vector(apply(modeArray, 2, mean));


  # NOTE!!! ONLY THE MULTIPLICATIVE METHOD IS CURRENTLY IMPLEMENTED!!!
  # Now, loop over the set again applying the needed corrections
  # Don't normalize the "Time" channel
  for (i in 1:length(fcsSet)) {
  
    multiFactor   = normVector / modeArray[i,];  # is this right??
#    addFactor     = normVector - modeArray[i,];    # double check this too.
#    print(paste("File", i));
#    print(multiFactor);
    thisFile      = fcsSet[[i]];
    thisData      = thisFile@data;
    for (j in 1:dim(thisData)[1]) {
      thisData[j,] = round(thisData[j,] * multiFactor);
#      thisData[j,] = round(thisData[j,] + addFactor);          # a different normalization
    }

    # HAH, NO NEED TO NORMALIZE THE "TIME" COLUMN
    if (normalizeTime == FALSE) {
      index = guessTimeChannel(thisFile);
      if (length(index) > 0) {
        thisData[,index] = thisFile@data[,index];
      } else {
        cat(paste("\n\n     I tried NOT to normalize time, but couldn't figure",
                   " out which column to avoid.  So I normalized all columns.\n\n"));
      }
    }
        
    thisFile@data = thisData;
    fcsSet[i]     = thisFile;   # need to do explicit replacement?
  }

#  return(modeArray);
  return(fcsSet);
}

