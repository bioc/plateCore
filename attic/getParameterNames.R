`getParameterNames` <-
function(fcsObj) {

  if (as.vector(class(fcsObj)) == "flowSet") {
    fcsObj = fcsObj[[1]];
  }

  if (as.vector(class(fcsObj)) == "flowFrame") {
    channelNames = as.vector((fcsObj)@parameters$desc);
    for (i in 1:length(channelNames)) {
      if (is.na(channelNames[i])) {
        channelNames[i] = as.vector((fcsObj)@parameters$name)[i];
      }
    }
    return(channelNames);
  }
  
  stop(paste("\n --> You tried to get the parameter names off of something",
             "that wasn't a 'flowFrame' or 'flowSet'.\n\n"));
}

