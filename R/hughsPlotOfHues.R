`hughsPlotOfHues` <-
function (fcsObj, chanOne=1, chanTwo=2) {

  if (is(fcsObj, "FCS") == FALSE) {
    stop("Hey, that isn't an FCS object you are trying to plot!");
    return(NULL);
  }

  paramNames  = getParameterNames(fcsObj);
  timeChannel = guessTimeChannel(fcsObj);
  dat         = fcsObj@data[,-c(chanOne, chanTwo, timeChannel)];
  dat.cov     = cov.rob(dat);     # this "robust" method may run slowly
  dat.mean    = apply(dat, 2, mean);
  dat.mahal   = mahalanobis(dat, dat.mean, dat.cov);
  
  plot(density(log(dat.mahal), bw=0.5),
       main="Log Mahalanobis distances");
  rug(log(dat.mahal));
  plot(fcsObj@data[,chanOne], fcsObj@data[,chanTwo],
       xlab=paramNames[chanOne], ylab=paramNames[chanTwo],
       pch=".", col=log(dat.mahal));

  return(dat.mahal);
}

