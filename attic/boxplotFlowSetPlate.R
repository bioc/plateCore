`boxplotFlowSetPlate` <-
function (aFlowSet, parameter="FSC-A", do.out=FALSE, ...) {
  pd         = phenoData(aFlowSet);
  cols       = as.vector(pd@data$Column);
  rows       = as.vector(pd@data$RowSymbol);
  plotOrder  = rev(order(rows, cols));
  wellNames  = paste(rows[plotOrder], cols[plotOrder], sep="");
  plateName = as.vector(unique(pd@data$PlateBarcode));

  dat        = list();
  par(mfrow=c(1,1));
  for (i in 1:length(plotOrder)) {
#  for (i in 1:12) {
    # lower whisker, lower hinge, median, upper hinge, upper whisker, outliers
    index  = plotOrder[i];
    dat[i] = list(boxplot.stats(aFlowSet[[index]]@exprs[,parameter], ...)$stats);
  }

  par(mfrow=c(1,1));
  boxplot(dat, names=wellNames, col=palette(rainbow(24))[cols],
          horizontal=TRUE,  main=plateName,
          srt=45, notch=TRUE, ylab="Wells", xlab=parameter);
          
  mat = matrix(unlist(dat), nrow=length(plotOrder), ncol=5, byrow=TRUE);
  abline(v=mean(mat[,3]), lty=2);
  
  invisible(dat);
}

