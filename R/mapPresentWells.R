`mapPresentWells` <-
function (wells, plateSize=96, inputPresent=TRUE,
                             returnType="vector") {

  wells      = toupper(wells);
  sentRows   = substr(wells, 1,1);
  sentCols   = as.numeric(sub("[A-Z]", "", wells));

  # Set up a matrix representing the plate
  rowNames   = vector();
  colNumbers = vector();
  numRows    = 8;

  if (plateSize == 96) {
    rowNames   = LETTERS[1:8];
    colNumbers = 1:12;
  }

  if (plateSize == 384) {
    rowNames   = LETTERS[1:16];
    colNumbers = 1:24;
    numRows    = 16;
  }

  # Create a return matrix representing the plate
  mat = matrix(data=rep(FALSE, plateSize), nrow=numRows, byrow=TRUE,
               dimnames=list(rowNames, colNumbers));

  # Alternatively, create a vector representing the plate
  vec = vector();
  for (r in rowNames) {
    for (cN in colNumbers) {
      index = paste(r, cN, sep="");
      vec[index] = FALSE;
    }
  }

  # Then "populate" the plate with the observed data
  for (i in 1:length(wells)) {
    index = paste(sentRows[i], sentCols[i], sep="");
    vec[index]                    = TRUE;
    mat[sentRows[i], sentCols[i]] = TRUE;
  }

  if (returnType == "matrix") { return(mat); }
  if (returnType == "vector") { return(vec); }
}

