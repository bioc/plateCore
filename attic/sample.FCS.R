`sample.FCS` <-
function(fcsObj, size=1, bootstrap=TRUE, change.time=FALSE, ...) {

  if (class(fcsObj) != "flowFrame") {
    print("Hey, thats not a flowFrame object, I can't sample it.");
    return();
  }

  numCells    = dim(fcsObj@exprs)[1];
  selection   = vector();
  params      = getParameterNames(fcsObj);
  timeChannel = guessTimeChannel(fcsObj);

  if (bootstrap == TRUE) {
    selection   = round((numCells - 1) * runif(numCells*size)) + 1;

  } else {
    selection = sample(1:numCells, size=size, ...);
    if (size >= numCells) {
      print(cat("\nWarning, this FCS object has ", numCells, " cell events. \n",
                "You are trying to select ", size, " cells.  \n",
                "Shouldn't you do a boostrap instead (ie set bootstrap==TRUE).\n\n", sep=""));
    }
  }


  # SOMETHING FUNNY HAPPENS IF YOU TRY A SAMPLE OF JUST ONE.
  if (length(selection) == 1) {
    resampled = t(fcsObj@exprs[selection,]);
  } else {
    resampled   = fcsObj@exprs[selection,];
  }

  if (change.time == TRUE) {
    timePoints              = as.vector(resampled[,timeChannel]);
    resampled[,timeChannel] = sample(timePoints, ...);
    resampled               = resampled[order(resampled[,timeChannel]),];
  }

  dimnames(resampled)[[1]] = seq(1:dim(resampled)[1]);
  fcsObj@exprs              = resampled;

  return(fcsObj);
}

