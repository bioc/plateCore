`facsBox360` <- 
function(aFlowFrame, params=c(1,2,3,4), method="polygons", ...) {

  # Choice of method = "polygons", "dots"
  pNames   = getParameterNames(aFlowFrame);
  fileName = as.character(keyword(aFlowFrame, "$FIL"));

  a     = stripGutters(aFlowFrame, method="ltZero");
  b     = stripGutters(a, method="quantile");
  mat   = b@exprs[, params];
  mat   = log(mat);
  dimnames(mat)[[2]] = pNames[params];
  dat   = data.frame(mat[,1:3]);
  mat   = mat[order(mat[,4]),];
  colorSet = rainbow(dim(mat)[1]);

  if (method == "dots") {
    open3d();
    plot3d(mat, col=colorSet, size=2, axes=FALSE,
          main=fileName,
          sub=paste("Colored by", pNames[params[4]]));

# Get the parameters like this
#  a = par3d();    # Don't think it works after the window is closed.
#  a$userMatrix;   # I think this is the viewing angle.
#  a$zoom;         # How zoomed in are you?


  } else if (method == "polygons") {
    featureSignif(dat, scaleData=TRUE, ...);
  }
}
