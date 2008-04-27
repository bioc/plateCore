`rank.FCS` <-
function(fcsObj, windoSize=3) {

  rankList   = list();
  
  if (is(fcsObj, "flowFrame") == FALSE) {
    stop("Hey, that isn't a flowFrame you are trying to rank!");
  }

  paramNames  = getParameterNames(fcsObj);
  timeChannel = guessTimeChannel(fcsObj);

  channels  = seq(1:length(paramNames));
  channels  = channels[-timeChannel];
  cellData  = fcsObj@exprs[,channels];
  
  # DERIVE A BASIC RANKING FOR THE CELLS FOR EACH CHANNEL
  ranks     = apply(cellData, 2, rank);       # eg. a matrix of numbers 1 ... 21317
  percRanks = ranks/dim(ranks)[1];            # now just 0...1

  # DERIVE A SCORE FOR A SINGLE CELL BASED ON ALL IT'S MEASUREMENTS
  prodRanks    = as.vector(apply(ranks, 1, prod));              # multiply the numbers across for each cell
  normRanks    = prodRanks / (dim(ranks)[1])^(dim(ranks)[2]);   # scale the values by the size of the total data
  negLogRanks  = -1 * log(normRanks);                           # convert to "better looking" numbers
  
  # DERIVE A SCORE FOR A SINGLE CELL AND IT'S NEIGHBORS.  THE ASSUMPTION HERE
  # IS THAT OCCASIONALLY YOU MAY FIND A DISRUPTION IN THE FLOW THAT CAUSES SEVERAL
  # CELLS IN A ROW TO BE "UNUSUAL".
  windowRanks = prodRanks;
  for (i in 1:length(windowRanks)) {
    leftBound  = round(i - (windowSize - 1)/2);
    rightBound = round(i + (windowSize - 1)/2);
    windowRanks[i] = prod(prodRanks[leftBound:rightBound], na.rm=FALSE);
  }

  # RESCALE BECAUSE OF ALL THE NUMBERS CHAINED TOGETHER
  normWindowRanks = windowRanks /(dim(ranks)[1])^(windowSize * (dim(ranks)[2]));
  logWindowRanks = -1 * log(normWindowRanks);


  # DERIVE A SCORE FOR EACH CELL EVENT BASED ON HOW MUCH TIME HAS ELAPSED
  # SINCE THE PREVIOUS CELL WENT THROUGH THE LASER.  I WOULD EXPECT THAT
  # THIS SHOULD BE A NORMAL DISTRIBUTION.
  deltaTime = deltaTime.FCS(fcsObj);
  
  # PREP THE DATA FOR RETURN
  rankList$ranks           = ranks;
  rankList$percRanks       = percRanks;
  rankList$negLogRanks     = negLogRanks;
  rankList$logWindowRanks  = logWindowRanks;
  rankList$deltaTime       = deltaTime;

  return(rankList);
}

