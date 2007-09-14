`getAFlowJoGate` <- 
function(gatePath=character(), flowJoDOM=newXMLDoc(),
                           quiet=FALSE, stripComp=TRUE) {

  # Recursive parsing of the FlowJo workspace using other functions above will yield
  # a list of gates.  Here we pick those off one at a time.  For example, if:
  # gatePath = "/Workspace/SampleList/Sample/SampleNode[@name='20070222c_A9_A09.fcs']/Subpopulations/Population[@name='Lymphocytes']/PolygonGate"
  # gateType will get "PolygonGate"
  # gateName will get "Lymphocytes"

  gateType       = gsub(".*/", "", gatePath);

  knownGateTypes = c("PolygonGate", "RectangleGate",
                     "RangeGate", "EllipsoidGate");
  if (!(gateType %in% knownGateTypes)) {
    if (quiet == FALSE) {
      print(paste("Sorry, I don't know what a '", gateType,
                "' gate is...", sep=""));
    }
    return(NULL);
  }

  # A clever trick to splice together the informative parts of the
  # name part of the xpath that give us the gate.
  partList    = unlist(strsplit(gatePath, split="'", perl=TRUE));
  goodParts   = seq(2, length(partList), by=2);
  gateName    = paste(partList[goodParts], collapse=":");


  # Note that FlowJo mysteriously bases all their FSC and SSC gates on values
  # divided by 64.  Hence I slip this in just in case.  Note, I haven't figured
  # out how to gracefully diagnose what type of scatter gate ("FSC-A" and "SSC-A",
  # or "FSC-H" and "SSC-H") are being used, and how to selectively transform on them.
  div64      <- linearTransform(transformationId="Linear-transformation", a=1/64, b=0);
  
  
  # Note below that FlowJo, when creating gates on compensated data, adds the
  # prefix "Comp-" to all of the gate parameter names.  I need to strip this
  # out in order to work cleanly in the flowCore paradigm.  Conversely, I may
  # sometimes need to track the fact of whether the gates are for compensated
  # or uncompensated data.
  
  # To fix log parameters, go to DivaSettings, look for "logAllFluors" and "linearParameters".

  # To find the transformations go to TransformSettings to see what to apply.  Read the
  # attributes in there.  Now NOTE, if one of the Transform elements has sampleName="GlobalDiVa Settings"
  # Then apply the Global DiVa Settings first, then the other transforms overlay them. (If applicable).
  
  # First case, a polygon gate.
  if (gateType == "PolygonGate") {
    axisPath  = paste(gatePath, "/Axis", sep="");
    axesNames = unlist(xpathApply(flowJoDOM, axisPath, xmlGetAttr, "name"));
    if (stripComp == TRUE) { axesNames = gsub("\\bComp-", "", axesNames, perl=TRUE); }

    pointPath = paste(gatePath, "/Point", sep="");
    xPoints   = as.numeric(unlist(xpathApply(flowJoDOM, pointPath, xmlGetAttr, "x")));
    yPoints   = as.numeric(unlist(xpathApply(flowJoDOM, pointPath, xmlGetAttr, "y")));
    
    verticies = matrix(c(xPoints, yPoints), ncol=2,
                         nrow=length(xPoints), byrow=FALSE);
                         
    colnames(verticies) = axesNames;
    
    # Do the gating on transformed data, if necessary.  I've got to figure out how
    # to do this more gracefully
    thisGate = polygonGate(filterId=gateName, boundaries=verticies);
    if (length(grep("FSC-H|SSC-H", axesNames)) > 0) {
      flowJo64  <- transform("FSC-H"=div64, "SSC-H"=div64);
      thisGate  <- polygonGate(filterId=gateName, boundaries=verticies) %on% flowJo64;
    
    } else if (length(grep("FSC-A|SSC-A", axesNames)) > 0) {
      flowJo64  <- transform("FSC-A"=div64, "SSC-A"=div64);
      thisGate  <- polygonGate(filterId=gateName, boundaries=verticies) %on% flowJo64;
    }
    
    # For some reason the filterId keeps getting overwritten when transforms are present
    thisGate@filterId = gateName;
    return(thisGate);
  }
  
  
  # Next case is a bit trickier.  FlowJo uses the term RectangleGate when
  # they could mean a 'range gate' or a 'rectangle gate' (ie made from two
  # subnested range gates.
  if (gateType == "RectangleGate") {
    rangePath   = paste(gatePath, "/RangeGate", sep="");
    axisPath    = paste(rangePath, "/Axis", sep="");
    rangeNames  = unlist(xpathApply(flowJoDOM, axisPath, xmlGetAttr, "name"));
    if (stripComp == TRUE) { rangeNames = gsub("\\bComp-", "", rangeNames, perl=TRUE); }

    rangeMins   = as.numeric(unlist(xpathApply(flowJoDOM, rangePath, xmlGetAttr, "min")));
    rangeMaxes  = as.numeric(unlist(xpathApply(flowJoDOM, rangePath, xmlGetAttr, "max")));

    verticies = matrix(c(rangeMins, rangeMaxes), ncol=2,
                         nrow=length(rangeMins), byrow=TRUE);

    colnames(verticies) = rangeNames;

    # Do the gating on transformed data, if necessary.
    thisGate = rectangleGate(filterId=gateName, verticies);

    if (length(grep("FSC-H|SSC-H", rangeNames)) > 0) {
       flowJo64  <- transform("FSC-H"=div64, "SSC-H"=div64);
       thisGate  <- polygonGate(filterId=gateName, boundaries=verticies) %on% flowJo64;

    } else if (length(grep("FSC-A|SSC-A", rangeNames)) > 0) {
       flowJo64  <- transform("FSC-A"=div64, "SSC-A"=div64);
       thisGate  <- polygonGate(filterId=gateName, boundaries=verticies) %on% flowJo64;
    }

    # For some reason the filterId keeps getting overwritten when transforms are present
    thisGate@filterId = gateName;
    return(thisGate);
  }

  # HAVEN'T FIGURED THESE OUT YET.
  if (gateType == "EllipseGate") {
    if (quiet == FALSE) { print("Sorry, I haven't worked out EllipseGates yet."); }
    return(NULL);
  }
}

