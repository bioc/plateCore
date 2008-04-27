`mode.FCS` <-
function (fcsObj) {

  if (class(fcsObj) != "flowFrame") {
    stop("Sorry, I can only get the mode on a flowFrame object.");
  }
  
  modeList    = vector();
  densityList = apply(fcsObj@exprs, 2, density);

  for (i in 1:length(densityList)) {
    d           = densityList[[i]];
    modeList[i] = d$x[d$y == max(d$y)];
  }
  
  names(modeList) = as.vector(unlist(colnames(fcsObj)));
  return(modeList);
}

