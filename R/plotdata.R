
#' @title plotdata
#' @description Plotting predicted and true response values together
#' @param x True values
#' @param y Predicted values
#' @param txt Main label to show
#' @param yr Range of values
#' @param grp A factor vector for grouped points
#' @export

plotdata <- function(x,y,txt="",yr=NULL,grp=NULL) {
  if(is.null(yr)) yr <- range(c(x,y))
  plot(x, y,ylab="Predicted response", xlab="True response",main=txt,xlim=yr,ylim=yr,axes=FALSE,ty="n")
  axis(1)
  axis(2)
  abline(a=0,b=1,col=1,lty=1,lwd=1.3)
  if(!is.null(grp)) {
    unGrp = levels(grp) #get unique grps
    if(is.null(unGrp)) unGrp = unique(grp) #get unique grps
    for(gg in unGrp) {
      ind = grp==gg
      points(x[ind], y[ind], col=which(gg==unGrp))
    }
    legend("topleft",legend=unGrp,col=1:length(unGrp),pch=1,bg="white")
  } else {
    points(x, y)
  }
}
