`getFlowJoCompMatrix` <-
function(flowJoDOM, silent=TRUE) {

  sampleXpath  = "/Workspace/SampleList/Sample";
  sampleNodes  = xpathApply(flowJoDOM, path=sampleXpath, xmlChildren);
  compMatrix   = matrix();
  compList     = list();
  if (silent == FALSE) {
    print("Be aware that FlowJo and flowCore have different notions about the the spill over matrix.");
    print("I have not yet set the getFlowJoCompMatrix method to do the 1/max automatically.");
  }
  
  # This loops of the samples in your Workspace.  There can theoretically be
  # as many different compensation matrices as there are samples.  Though in
  # practice this doesn't happen very often.  Under each "Sample" node are
  # 'Keywords', 'DataSet', 'CompensationMatrix' and 'SampleNode'
  for (i in 1:length(sampleNodes)) {
    aSampleNode     = sampleNodes[[i]];
    sampleUriIndex  = grep("DataSet", lapply(aSampleNode, xmlName));
    sampleName      = xmlGetAttr(aSampleNode[[sampleUriIndex]], "uri");
    sampleName      = basename(sampleName);
    compIndex       = grep("CompensationMatrix", lapply(aSampleNode, xmlName));
    assemblyList    = list();
    
    # Burrow into the "CompensationMatrix" node.
    for (j in compIndex) {
      aSubNode     = aSampleNode[[j]];
      channelNodes = xmlChildren(aSubNode);
      channelIndex = grep("Channel", lapply(channelNodes, xmlName));

      # Walk through each of the "Channel" nodes.
      for (k in channelIndex) {
        channelName = xmlGetAttr(channelNodes[[k]], "name")
        overLaps    = xmlChildren(channelNodes[[k]]);
#        print(paste(channelName, length(overLaps)));

        # And collect the overlap values in a growing list.
        for (m in grep("ChannelValue", lapply(overLaps, xmlName))) {
          overlapName  = xmlGetAttr(overLaps[[m]], "name");
          overlapValue = xmlGetAttr(overLaps[[m]], "value");
#          print(paste("   ", overlapName, overlapValue));
          assemblyList[[channelName]][[overlapName]] = overlapValue;
        }
      }
    }

    # Then convert the "list" into a "matrix".  Note, I could have done this
    # earlier but I'm not sure that the order of the channel names will be
    # the same in each "Channel" node.
    numChannels     = length(assemblyList);
    channelNamesOrd = names(assemblyList);
    assemblyMat     = matrix(nrow=numChannels, ncol=numChannels, byrow=FALSE,
                          dimnames=list(channelNamesOrd, channelNamesOrd));
                          
    for (topChannel in channelNamesOrd) {
      for (leftChannel in channelNamesOrd) {
        assemblyMat[leftChannel, topChannel] = as.numeric(assemblyList[[leftChannel]][[topChannel]]);
      }
    }
    
    # Note that what finally gets returned is a list of compensation matricies
    compList[[sampleName]] = assemblyMat;
  }

  return(compList);
}

