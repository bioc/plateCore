`flowSetFromFlowJo` <-
function(xmlFile, autoCompensate=TRUE, ...) {

  flowJoDOM       = readFlowJoFile(xmlFile);
  dataSetPath     = "/Workspace/SampleList/Sample/DataSet";
  uris            = xpathApply(flowJoDOM, path=dataSetPath, xmlGetAttr, "uri");
  uris            = unique(unlist(uris));

  # Normally the URI attribute should hold the full path to the FCS file.  However,
  # sometimes it will hold a relative path.  This tries to catch that problem.
  if (sum(file.exists(uris)) != length(uris)) {
    path = dirname(xmlFile);
    uris = sub("\\.", "", uris);
    uris = paste(path, uris, sep="");
  }

  # Otherwise, just give up.
  if (sum(file.exists(uris)) != length(uris)) {
    stop("I can't find the FCS files referenced in the FlowJo workspace!");
  }

  aFlowSet        = read.flowSet(uris, ...);
  compCheck       = lookForCompensation(flowJoDOM);
  didCompensation = FALSE;
  
  if (compCheck$allHaveCompMats & compCheck$sameCompMats & autoCompensate) {
    spillOverMatrix  = getFlowJoCompMatrix(flowJoDOM, silent=FALSE)[[1]];

    print("I am now applying 1/max(compMatrix) to the data to make the concept");
    print("of 'compensation' in flowCore the same as 'compensation' in FlowJo.");
    spillOverMatrix  = spillOverMatrix / max(spillOverMatrix);  # Funny FlowJo thing
    aFlowSet         = compensate(aFlowSet, spillOverMatrix);
    didCompensation  = TRUE;
  }
  
  if (compCheck$allHaveCompMats & (compCheck$hasCompGates == FALSE)) {
    cat("\n\tWARNING: I see non-trivial compensation matrices in the FlowJo workspace,\n");
    cat("\tbut I don't see any of the gates labelled with 'Comp-' in them as FlowJo\n");
    cat("\tseems to normally do in such cases.\n\n");
  }

  
  if (compCheck$sameCompMats == FALSE) {
    cat("\n\tWARNING: I noticed that different compensation matrices are stored for different FCS\n");
    cat("\tfiles in this FlowJo workspace.  This is unexpected?!? I will NOT compensate right now.\n\n");
  }

  if (didCompensation == FALSE) {
    cat("\n\tNo compensation was performed on this flow data.\n\n");
  } else {
    cat("\n\tThe flow data were compensated on the way in.\n\n");
  }

  gc();
  return(aFlowSet);
}

