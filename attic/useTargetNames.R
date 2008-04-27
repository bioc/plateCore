`useTargetNames` <-
function(fcsObj) {

  if (as.vector(class(fcsObj)) == "flowFrame") {
    dimnames(fcsObj@exprs)[[2]] <- getParameterNames(fcsObj);
    return(fcsObj);
  }

  if (as.vector(class(fcsObj)) == "flowSet") {
    fcsObj@colnames  <- getParameterNames(fcsObj);
    fcsObj = fsApply(fcsObj, useTargetNames);
    return(fcsObj);
  }
  stop("You didn't give me a flowFrame or flowSet to swap names on :-(");
}

