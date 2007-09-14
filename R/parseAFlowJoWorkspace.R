`parseAFlowJoWorkspace` <- 
function (xmlFile, keep="all", ...) {

  # keep       = c("all", "leaf", "internal");
  filterVec     = vector();
  filterNameVec = vector();
  fcsPathVec    = vector();
  fcsNameVec    = vector();
  resultList    = list();
  
  
  # Load the FlowJo workspace in as a DOM tree
  flowJoDOM   = readFlowJoFile(xmlFile);
  gatePaths   = recurseFlowJoGates(flowJoDOM=flowJoDOM, keepMiddleGates=FALSE);
  
  # Will want to keep track of where the corresponding FCS files live
  dataSetPath  = "/Workspace/SampleList/Sample/DataSet";
  uris         = unlist(xpathApply(flowJoDOM, path=dataSetPath, xmlGetAttr, "uri"));
  uris         = unique(uris);
#  print(uris);
  
  # Note that the recursion process returns a list that contains the xpath
  # to both the 'Gate' and 'Population' level branches in the FlowJo DOM.
  for (i in 1:length(gatePaths)) {
  
    # Discard the gatePaths that don't end in gates.  Could have done
    # this with a regexpr...
    gatePath      = gatePaths[i];
    originalGate  = getAFlowJoGate(gatePath, flowJoDOM, quiet=TRUE, ...);
    if (is.null(originalGate)) {
      next;  # ignore the xpaths that don't end in gates

    # Now look back up every level of the tree and pull out the sections
    # that correspond to "parent" gates off of this gate.
    } else {
      subGatesVec  = vector();
      subXpathVec  = vector();
      fcsNameIndex = regexpr("'[^']*'", gatePath, perl=TRUE);
      fcsName      = substring(gatePath, first=fcsNameIndex+1,
                           last=(fcsNameIndex + attr(fcsNameIndex, "match.length") -2));
                           
      fcsPath      = uris[grep(fcsName, uris)];
                           
      subPopIndex = unlist(gregexpr("Subpopulations", gatePath, perl=TRUE));
      for (i in 1:length(subPopIndex)) {
        aParentXpath = substring(gatePath, first=1, last=(subPopIndex[i] - 2));

        nextChildren = xpathApply(flowJoDOM, aParentXpath, xmlChildren);
        for (j in 1:length(nextChildren[[1]])) {
          nodeName = xmlName(nextChildren[[1]][[j]]);
          if (regexpr("Gate", nodeName) > 0) {
#            print(paste("See a", nodeName, "at", aParentXpath));
            goodGatePath = paste(aParentXpath, nodeName, sep="/");
            parentGate   = getAFlowJoGate(goodGatePath, flowJoDOM, quiet=TRUE, ...);
            subGatesVec  = append(subGatesVec, parentGate);
#            subXpathVec  = append(subXpathVec, goodGatePath);
          }
        }
      }


      # Assemble these in order from most general to most specific gate
      subGatesVec   = append(subGatesVec, originalGate);
#      subXpathVec   = append(subXpathVec, gatePath);
      filterChain   = subGatesVec[[1]];
      
      for (g in 1:length(subGatesVec)) {
        if (g > 1) { filterChain = filterChain & subGatesVec[[g]]; }
        filterName      = subGatesVec[[g]]@filterId;
        filterChain@filterId = filterName;
        filterNameVec   = append(filterNameVec, filterName);
        filterVec       = append(filterVec, filterChain);
        fcsNameVec      = append(fcsNameVec, fcsName);
        fcsPathVec      = append(fcsPathVec, fcsPath);
      }
    }
    
  } # Now go back and find the next leaf node gate
  
  resultList[[1]] = fcsNameVec;
  resultList[[2]] = filterNameVec;
  resultList[[3]] = filterVec;
  resultList[[4]] = fcsPathVec;
  names(resultList) <- c("fcsName", "filterName", "filter", "fcsPath");
  
  gc();
  return(resultList);
}

