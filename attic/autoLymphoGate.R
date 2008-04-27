`autoLymphoGate` <-
function (aFlowFrame, fsc="FSC-A", ssc="SSC-A") {

  filterId = "autoGatedLymphocytes";
  nf = norm2Filter(filterId=filterId,
                   x=c(fsc, ssc), method="covMcd",
                   scale.factor=1.0);

  return(Subset(aFlowFrame, nf));
}

