plotWellCounts.FCS <- function(aFlowSet) {

    # Create a plot showing the number of events in each raw FCS
    # file for each well on the plate.
    pd        = phenoData(aFlowSet);
    counts    = unlist(fsApply(aFlowSet, keyword, "$TOT"));
    index     = grep("\\.\\$TOT1$", names(counts), perl=TRUE);
    runCounts = as.numeric(as.vector(counts[index]));
    plateName = as.vector(unique(pd@data$PlateBarcode));
    
    # Put it all in a good order
    plateMap    = setUpPlateMap(aFlowSet);
    pointLabels = plateMap$wellNames[plateMap$plotOrder];
    runCounts   = runCounts[plateMap$plotOrder];
    
    plotPlate(runCounts, ind = plateMap$present, 
        main = paste("Total Cell Events for", plateName), 
        col = c("yellow", "darkblue"), cex.main = 1, cex.desc = 1, 
        desc = c("Lots", "Few"));    
}

