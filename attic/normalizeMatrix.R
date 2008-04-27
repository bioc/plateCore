`normalizeMatrix` <-
function (mat) {
  for (i in 1:dim(mat)[2]) {
    colRange = range(mat[,i]);
    colScale = colRange[2] - colRange[1];
    mat[,i] = (mat[,i] - colRange[1]) / colScale;
  }
  return(mat);
}

