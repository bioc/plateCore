`textPlateMap` <-
function (myData = vector(), sigFigs = 3) {

  # Guess 96 or 384.
  # Add row and column headers
  # Round as appropriate
  # Deal with NAs appropriately.

  if (class(myData) == "numeric") {
    myData  = signif(myData, digits=sigFigs);
  } else if (class(myData) == "character") {
    myData  = substr(myData, start=1, stop=sigFigs);
  }
  
  # Set up defaults
  rowNames   = vector();
  colNumbers = vector();
  
  if (length(myData) == 96) {
    rowNames   = LETTERS[1:8];
    colNumbers = 1:12;
    numRows    = 8;

  } else if (length(myData) == 384) {
    rowNames   = LETTERS[1:16];
    colNumbers = 1:24;
    numRows    = 16;

  } else {
    stop("I don't know how to handle that data.");
  }

  mat = matrix(myData, nrow=numRows, byrow=TRUE,
               dimnames=list(rowNames, colNumbers));
  return(mat);
}

