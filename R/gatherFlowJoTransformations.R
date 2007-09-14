`gatherFlowJoTransformations` <-
function(xmlFile) {

  returnList = list();

  # The DivaSettings node in the FlowJo Workspace contains general transformation
  # parameters for the relevant FCS files.  These can be overwritten if there are
  # specific transformations for specific FCS files and parameters (see below).
  #
  # My transcript of the conversation w/ the TreeStar folks (Aaron Hart ahart@treestar.com):
  # 1. Look under <TransformSettings>.  
  #      If the only transform there is "Global Diva Settings", then *don't* apply it
  #      If there are other transforms available, then apply the "Global Diva Settings"
  #        first, and override them (as appropriate) with the other settings.
  #
  # 2. If there is a single transformation, apply that one.
  #
  # * 8/2/07 In retrospect, these instructions don't seem to make sense.  In any
  # event, the code below does *not* reflect these instructions.
  # 

  flowJoDOM       = readFlowJoFile(xmlFile);
  divaXpath       = "/Workspace/DivaSettings";
  divaNode        = getNodeSet(flowJoDOM, divaXpath);
  
  divaVec    = vector();
  divaAttrs  = c("numDecades", "logSideSc", "logForwSc", "logAllNonSc", "logAllFluors", "baseChannel",
                 "pulseHeight", "pulseArea", "autoTransform", "numDecades", "extraNegs",
                 "widthBasis", "linearParameters");
  
  for (i in 1:length(divaAttrs)) {
    if (is.null(xmlGetAttr(divaNode[[1]], divaAttrs[i]))) {
      divaVec[i] = NA;
    } else {
      divaVec[i] = xmlGetAttr(divaNode[[1]], divaAttrs[i]);
    }
  }
  names(divaVec) <- divaAttrs;
  
  # In addition, there may be one or more Transform elements that describe
  # transformations on specific parameters for specific flow files.
  transformXpath  = "/Workspace/TransformSettings/Transform";
  transformNodes  = getNodeSet(flowJoDOM, transformXpath);
  numTransforms   = length(transformNodes);
  transAttrs      = c("sampleName", "parameterName", "compensated", "extraNegDecades",
                      "posDecades", "widthBasis", "isDiva", "autoCalculateExtraNegative");

  # If the FlowJo file has transforms then loop over all of them and put them
  # into an array.
  if (numTransforms > 0) {
    sampleNames     = unlist(lapply(transformNodes,xmlGetAttr, "sampleName"));
    transformArray  = array(dim=c(numTransforms, length(transAttrs )),  
                            dimnames=list(NULL, transAttrs));

    for (i in 1:numTransforms) {
      aNode         = transformNodes[[i]];
      isCompensated = FALSE;
#      parameterName = xmlGetAttr(aNode, "parameterName");

      for (transAttr in transAttrs) {
        value = xmlGetAttr(aNode, transAttr);
        if (is.null(value)) { value = NA; }
        if (transAttr == "compensated") { value = isCompensated; }
        
        if ((transAttr == "parameterName") &
            (length(grep("Comp-", value)) > 0)) {
            value = gsub("Comp-", "", value);
            isCompensated = TRUE;
        }
        transformArray[i, transAttr]  = value;
      }
    }
    
  } else {
    # No transform nodes
    transformArray = NA;
  }

  returnList[[1]]   = divaVec;
  returnList[[2]]   = transformArray;
  names(returnList) = c("DivaSettings", "TransformSettings");
  gc();
  return(returnList);
}

