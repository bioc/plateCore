`rotateMatrix` <-
function(mat) {

  orgRows = dim(mat)[1];
  orgCols = dim(mat)[2];
  newMat  = matrix(nrow=orgCols, ncol=orgRows);
  j = 0;
  for (i in orgRows:1) {
    j = j+1;
    newMat[,j] = mat[i,];
  }
  dimnames(newMat) = list(dimnames(mat)[[2]], rev(dimnames(mat)[[1]]));
  return(newMat);
}

