`readFlowJoFile` <-
function(xmlFile) {

  ## Initial load of XML data file into in-RAM DOM tree
  flowJoDOM <- xmlTreeParse(xmlFile, useInternalNodes=TRUE)

  ## Get the root of the DOM tree; expect it to be a "Workspace" node
  root <- xmlRoot(flowJoDOM);
  if(xmlName(root) != "Workspace") {
    stop("Top node is not a 'Workspace' node.  Is this really a FlowJo workspace?!?");
  }

  ## This Workspace should refer to at least one or more FCS files
  numSamples = length(xpathApply(flowJoDOM,
                                 path="/Workspace/SampleList/Sample/DataSet",
                                 xmlChildren));
                                 
  if (numSamples < 1) {
    stop("This workspace doesn't refer to any FCS files.  I give up!");
  }
  return(flowJoDOM);
}

