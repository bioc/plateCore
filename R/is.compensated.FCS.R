`is.compensated.FCS` <-
function(fcsObj) {

  if (class(fcsObj) == "FCS") {
    fcsObj = list(fcsObj);
  }

  answerMat = matrix(nrow=length(fcsObj), ncol=dim(fcsObj[[1]])[2]);

  for (i in 1:length(fcsObj)) {
    answerMat[i,] = apply(fcsObj[[i]]@data, 2, is.integer);
  }
  return(answerMat);
}

