`deltaTime.FCS` <-
function(aFlowFrame, orderTime=TRUE) {

  if (is(aFlowFrame, "flowFrame") == FALSE) {
    stop("Hey, I need a flowCore style flowFrame object for this to work!");
    return(NULL);
  }

  timeChannel  = guessTimeChannel(aFlowFrame);
  timePoints   = as.vector(aFlowFrame@exprs[,timeChannel]);

  # Make sure the data points are in the correct time order.  Note, certain
  # sampling steps may have scrambled this.  Be aware, though, if you resort
  # stuff here, then these values won't be in the same order as the data in
  # the input flowFrame...
  if (orderTime == TRUE) {
    timePoints = timePoints[order(timePoints)];
  }

  deltaTime = c(timePoints[1], diff(timePoints));

  return(deltaTime);
}

