`stripGutters` <-
function (fcsObj, gutters=list(), method="gutters", quantiles=c(0.05, 0.95)) {

  # Strip the extremes from a flowFrame.  
  # method = "gutters", "ltZero", "quantile"

  if (class(fcsObj) != "flowFrame") {
    stop("I was expecting a flowFrame object.");
  }

  keeperIndex = vector();
  goodMat     = matrix(nrow=dim(fcsObj@exprs)[1], ncol=dim(fcsObj@exprs)[2]);
  
  # Get rid of just the most extreme observed values (min and max)
  if (method == "gutters") {
     stats       = summary(fcsObj);
     paramMins   = stats["Min.",];  # + .Machine$double.eps;
     paramMaxes  = stats["Max.",];  # - .Machine$double.eps;
  
     for (i in 1:dim(fcsObj@exprs)[2]) {
        vec = ((fcsObj@exprs[,i] > paramMins[i]) &
            (fcsObj@exprs[,i] < paramMaxes[i]));
        goodMat[,i] = vec;
     }

  # Just get rid of points that are less than zero in any of the channels
  } else if (method == "ltZero") {
     for (i in 1:dim(fcsObj@exprs)[2]) {
        vec = (fcsObj@exprs[,i] > 0);
        goodMat[,i] = vec;
     }

  # Strip values out by their quantiles.  Default is to strip out everything
  # outside of the 5th and 95th percentiles.  
  } else if (method == "quantile") {  
     for (i in 1:dim(fcsObj@exprs)[2]) {
        limits = quantile(fcsObj@exprs[,i], probs=quantiles);
        vec = ((fcsObj@exprs[,i] > limits[1]) &
            (fcsObj@exprs[,i] < limits[2]));
        goodMat[,i] = vec;
     }
  }

  # KEEP ONLY THE CELL EVENTS (ROWS) WHERE ALL OF ITS PARAMETERS ARE
  # INSIDE OF THE MIN AND MAX LIMITS.
  keeperIndex = apply(goodMat, 1, all);

  retFcsObj        = fcsObj;
  retFcsObj@exprs  = fcsObj@exprs[keeperIndex,];
  gc();
  return(retFcsObj);
}

