`lookForCompensation` <-
function (flowJoDOM) {

  compList                = getFlowJoCompMatrix(flowJoDOM);
  answer                  = list();  # default to all FALSE
  answer$allHaveCompMats  = FALSE;
  answer$sameCompMats     = FALSE;
  answer$hasCompGates     = FALSE;

  compMatCount  = 0;
  diffCompCount = 0;
  compGateCount = 0;
  
  for (i in 1:length(compList)) {
    # Track number of non trivial individual compensations
    if (sum(compList[[i]] != 1) > 1) {
      compMatCount = compMatCount + 1;
    }
    
    # Track  if the compensation matrices aren't all the same (as the first one)
    if (! isTRUE(all.equal(compList[[1]], compList[[i]]))) {
      diffCompCount = diffCompCount + 1;
    }
  }
  
  # Count how many gates start with the word "Comp-" for their parameters.
  axisNodes     = getNodeSet(flowJoDOM, "//PolygonGate/Axis");
  axisNames     = unlist(lapply(axisNodes, xmlGetAttr, "name"));
  rangeNodes    = getNodeSet(flowJoDOM, "//RangeGate/Axis");
  rangeNames    = unlist(lapply(rangeNodes, xmlGetAttr, "name"));
  compGateCount = length(grep("\\bComp-", c(axisNames, rangeNames), perl=TRUE))

  # So what is the final answer for this FlowJo workspace?
  if (length(compList) == compMatCount) { answer$allHaveCompMats  = TRUE; }
  if (diffCompCount == 0) {               answer$sameCompMats     = TRUE; }
  if (compGateCount > 0) {                answer$hasCompGates     = TRUE; }

  return(answer);
}

