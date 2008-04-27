`evenFCSSet` <-
function (fcsSet, trim=TRUE) {

  firstChannelSet = getParameterNames(fcsSet[[1]]);
  minObsArray     = array(dim=dim(fcsSet)[1], dimnames=list(NULL, firstChannelSet));
  maxObsArray     = minObsArray;
  
  for (i in 1:dim(fcsSet)[2]) {
    thisFile = fcsSet[[i]];
    
        # They should all have the same list of parameters, in the same order
    if (all.equal(getParameterNames(thisFile), rightChanneOrder) != TRUE) {
      stop(paste("Error!   File", i,
                 "doesn't have the same list of parameters as the first file"));
    }

    # loop over every channel and record the min and max possible
    minObsArray[i,] = apply(thisFile@data, 2, min);
    maxObsArray[i,] = apply(thisFile@data, 2, max);
  }
  minVec = apply(minObsArray, 2, max);   # max of the mins
  maxVec = apply(maxObsArray, 2, min);   # min of the maxes

  # HERE, I SHOULD EITHER TRIM THEM, OR REPORT DISCREPANCIES
  # FOR EACH CHANNEL, DELETE ALL OBSERVATIONS THAT ARE:
  #   - LESS THAN THE BIGGEST MINIMUM
  #   - GREATER THAN THE SMALLEST MAXIMUM
  for (i in 1:dim(fcsSet)[2]) {
    thisFile = fcsSet[[i]];
    thisData = thisFile@data;
    
    for (i in 1:dim(thisData)[1]) {
      keeperIndex[i] = all((thisData[i,] < maxVec) &
                           (thisData[i,] > minVec));
    }
    
    thisFile@data = thisData[keeperIndex,];
    fcsSet[i] = thisFile;
  }
  
  return(fcsSet);
}

