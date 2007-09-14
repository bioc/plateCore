`guessTimeChannel` <-
function (fcsObj) {

  timeChannel = NULL;
  if (is(fcsObj, "flowSet") == TRUE) {
    return(guessTimeChannel(fcsObj[[1]]));
  }

  if (is(fcsObj, "flowFrame") == FALSE) {
    stop("Hey, that isn't an FCS object you are trying to find the time channel on.");
    return(timeChannel);
  }

  paramNames  = getParameterNames(fcsObj);
  timeChannel = grep("Time", paramNames, ignore.case=TRUE);

  # Last ditch, find the channel that is always increasing
  if (length(timeChannel) < 1) {
    maxRow = min(30, dim(fcsObj@exprs)[1]);
    dat    = fcsObj@exprs[1:maxRow,];  # the FCS data

    # Check each channel, stop when you get to the one always increasing
    for (column in 1:dim(dat)[2]) {
      # diffVec should always be >= 0 in the time column
      diffVec = dat[2:maxRow,column] - dat[1:(maxRow-1),column];
      if (sum(diffVec >= 0) == length(diffVec)) {
        timeChannel = column;
        return(timeChannel);
      }
    }
  }
  return(timeChannel);
}

