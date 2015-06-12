# Purpose: MVTB Plots (Revised)
# Author: Patrick Miller
# RRV: June 12, 2015
# Purpose: plots

#' Plots the model implied effect of 1 predictor for one outcome.
#' 
#' @param out mvtb output object
#' @param response.no index of the response variable
#' @param predictor.no index of the predictor variable
#' @param n.trees desired number of trees (default: best trees)
#' @param X predictor. If included, a rug is included with the plot showing the density of the variable.
#' @param xlab label of the x variable
#' @param ... extra arguments are passed to plot. See ?par
#' @return Function is called for it's side effect, a plot.
#' @export
mvtb.plot1 <- function(out,response.no,predictor.no,n.trees=min(unlist(out$best.trees)),X=NULL,xlab=NULL,...){
  ri <- relative.influence(gbm.obj,n.trees=n.trees)/sum(relative.influence(gbm.obj,n.trees=n.trees))*100
  ri <- ri[predictor.no]
  #gbm.obj <- convert.mvtb.gbm(out,k=response.no)
  gbm.obj <- out$models[[k]] 
  if(is.null(xlab)){ xlab <- names(ri)}
  xlab <- paste0(xlab," ", formatC(ri,2), "%")
  grid <- plot.gbm(gbm.obj,i.var = predictor.no,n.trees = bi,perspective=TRUE,return.grid=TRUE)
  plot(y=grid$y,x=grid[,1],bty="n",xlab=xlab,...)
  if(!is.null(X)) { rug(jitter(X[,predictor.no])) }
}


#' Perspective plot for 2 predictors and 1 response. 
#' 
#' This is a plot of the model implied function of 2 predictors averaged over the other predictors
#' included in the model. This is called a partial dependence plot.
#' As an alternative to the perspective (3D) plot, a 2D heat plot can be obtained directly
#' using ?plot.gbm.
#' 
#' @param out mvtb output object
#' @param response.no index of the response variable
#' @param predictor.no vector containing indeces of the predictor variables to plot
#' @param n.trees desired number of trees (default: best trees)
#' @param phi angle of viewing direction. See ?persp.
#' @param theta angle of viewing direction See ?persp.
#' @param r distance from eye to center. See ?persp
#' @param d strength of perspective. See ?persp. 
#' @param ticktype 'detailed' gives axis points. See ?persp for other options.
#' @param ... extra arguments are passed to persp. See ?persp
#' @return Function is called for it's side effect, a plot.
#' @seealso plot.gbm 
#' @export
mvtb.perspec <- function(out,response.no,predictor.no,n.trees=min(unlist(out$best.trees)),
                         phi=15,theta=-55,r=sqrt(10),d=3,ticktype="detailed",...) {
  #gbm.obj <- convert.mvtb.gbm(out,k=response.no)
  gbm.obj <- out$models[[k]]
  grid <- plot.gbm(gbm.obj,i.var = predictor.no,n.trees = n.trees,perspective=TRUE,return.grid=TRUE)
  x <- unique(grid[,1])
  y <- unique(grid[,2])
  z <- matrix(grid[,3],length(unique(x)),length(unique(y)))
  persp(x=as.numeric(x),y=as.numeric(y),z=z,d=d,r=r,phi=phi,theta=theta,ticktype=ticktype,...)
}


# Pairwise plot for 2 predictors and 1 response. 
plot.pw.perspec <- function(out,response.no,predictor.no,npairs=3,nonlin.rank=NULL,p1=NULL,p2=NULL,theta=rep(-55,npairs),...){
  bi <- out$best.iter[[2]]
  pred.names <- out$iter.models[[1]][[1]]$var.names
  if(is.null(nonlin.rank)){
    ris <- sort(out$ri[[2]][,response.no],decreasing=T)
    if(is.null(p1)) p1 <- rep(match(names(ris[predictor.no]),pred.names),npairs)
    if(is.null(p2)) p2 <- match(names(ris[(predictor.no+1):(predictor.no+npairs+1)]),pred.names)
  } else {
    r <- nonlin.rank[[response.no]]
    if(is.null(p1)) p1 <- r$var1.index[predictor.no:(predictor.no+npairs-1)]
    if(is.null(p2)) p2 <- r$var2.index[predictor.no:(predictor.no+npairs-1)]
  }
  for(i in 1:npairs) {
    mvtb.perspec(out,response.no=response.no,predictor.no=c(p1[i],p2[i]),
                 xlab=pred.names[p1[i]],ylab=pred.names[p2[i]],theta=theta[i],...)
  }
}

#' Simple heatmap of the covariance explained matrix.
#' 
#' @param out mvtb output object
#' @param clust.method clustering method for rows and columns. See ?hclust
#' @param dist.method  method for computing the distance between two lower triangluar covariance matrices. See ?dist for alternatives.
#' @param numformat function to format the covex values into strings. Defaults to removing leading 0 and rounding to 2 decimal places.
#' @param ... extra arguments are passed to image, then to plot. See ?image, ?par
#' @return heatmap of the clustered covariance matrix.
#' @export 
heat.covex <- function(out,clust.method="ward.D",dist.method="manhattan",numformat=function(val){sub("^(-?)0.", "\\1.", sprintf("%.2f", val))},...) {
  x <- out$covex
  hcr <- hclust(dist(x,method=dist.method),method=clust.method)
  ddr <- as.dendrogram(hcr)
  rowInd <- order.dendrogram(ddr)
  hcc <- hclust(dist(t(x),method=dist.method),method=clust.method)
  ddc <- as.dendrogram(hcc)
  colInd <- order.dendrogram(ddc)
  x <- x[rowInd,colInd]
  cellnote <- matrix(numformat(x),dim(x))
  cellnote <- cellnote[rowInd,colInd]
  x <- t(x)
  nc <- nrow(x) # final number of columns (usually predictors)
  nr <- ncol(x) # final number of rows    (usually dvs)
  cellnote <- t(cellnote)
  image(x=1:nc,y=1:nr,x,col=wbr(50),xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr),ylab="",xlab="",axes=F,...)
  #axis(1,at=seq(0,1,length=nrow(x)))
  cexRow <- .2+1/log10(nc)
  axis(1, 1:nc, labels = rep("",nc), las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  axis(2, 1:nr, labels = colnames(x), las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  text(x =  c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
       col = "white", cex = 1)
  text(1:nc,rep(0,nc), las=2,cex.axis=cexRow,adj=1,
       labels = rownames(x), xpd = TRUE,srt=45,
       col="black")
}
