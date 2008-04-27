`flowSetMahalanobis` <-
function (aFlowSet, ignore.time=TRUE, method="median") {

  if (is(aFlowSet, "flowSet") == FALSE) {
    stop("Hey, I need a flowCore style flowSet object for this to work!");
    return(NULL);
  }

  pd          = phenoData(aFlowSet);
  paramNames  = as.vector(colnames(aFlowSet));
  mat         = matrix(nrow=length(aFlowSet), ncol=length(paramNames),
                       dimnames=list(pd@data$name, paramNames));
  
  if (method == "median") {
    mat = fsApply(aFlowSet, each_col, median);
    
  } else if (method == "mode") {
    mat = fsApply(aFullPlate, mode.FCS)
  }
  
  if (ignore.time == TRUE) {
    timeChannel = guessTimeChannel(aFlowSet);
    mat         = mat[,-timeChannel];
  }
  
  mat.cov     = cov.rob(mat);     # this "robust" method may run slowly
  mat.mean    = apply(mat, 2, mean);
#  mat.mahal   = mahalanobis(mat, mat.mean, mat.cov$cov);
  mat.mahal   = mahalanobis(mat, mat.cov$center, mat.cov$cov, method="mcd");

  return(mat.mahal);
}

