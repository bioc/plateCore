`recurseFlowJoGates` <-
function(nodeVec = vector(),
                               flowJoDOM = newXMLDoc(),
                               keepMiddleGates = TRUE) {

  if (length(nodeVec) == 0) { nodeVec[1] = "/Workspace/SampleList/Sample/SampleNode"; }
  prevPath      = nodeVec[length(nodeVec)];
  childrenNames = xpathApply(flowJoDOM, path=prevPath, xmlGetAttr, "name");
  if (length(childrenNames) == 0) { return(nodeVec); }

  for (i in 1:length(childrenNames)) {
  
    nextSubPopLayer = paste(prevPath, "[@name='",
                            childrenNames[i],
                            "']/Subpopulations/Population", sep="");

    layersAhead = length(xpathApply(flowJoDOM, nextSubPopLayer, xmlChildren));
#    print(paste("Layers ahead:", layersAhead));

    # Now add the xpaths to upcoming gates.  Note in some cases I want to add
    # every gate I run across, other times, only the terminal gates.
    if ((keepMiddleGates == TRUE) |
        ((keepMiddleGates == FALSE) & (layersAhead == 0)) ) {
      nextLayer    = paste(prevPath, "[@name='", childrenNames[i], "']", sep="");
      nextChildren = xpathApply(flowJoDOM, nextLayer, xmlChildren);
      for (j in 1:length(nextChildren[[1]])) {
        nodeName = xmlName(nextChildren[[1]][[j]]);
#        print(paste("See a", nodeName, "at", nextLayer));
        if (regexpr("Gate", nodeName) > 0) {
          nodeVec = append(nodeVec, paste(nextLayer, nodeName, sep="/"));
        }
      }
    }

    # Then keep burrowing down further
    nodeVec = append(nodeVec, nextSubPopLayer);
    nodeVec = recurseFlowJoGates(nodeVec, flowJoDOM, keepMiddleGates);
  }
  
  return(nodeVec);
}

