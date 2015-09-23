# Purpose: MVTB Plots (Revised)
# Author: Patrick Miller
# RRV: June 12, 2015
# Purpose: plots

#' Plots the model implied effect of 1 predictor for one outcome 
#' @param x mvtb output object
#' @param predictor.no index of the predictor variable
#' @param response.no index of the response variable
#' @param n.trees desired number of trees (default: best trees)
#' @param X predictor. If included, a rug is included with the plot showing the density of the variable.
#' @param xlab label of the x axis
#' @param ylab label of the y axis
#' @param ... extra arguments are passed to plot. See ?par
#' @return Function does not return a value, but produces a plot of the model implied effect along with the relative influence of the predictor.
#' @seealso \code{plot.gbm}, \code{mvtb.perspec}, for other plots, \code{mvtb.heat} to plot the covariance explained by predictors in a heatmap
#' @export
#' @importFrom graphics plot rug
#' @details 
#' This is the classic partial dependence plot, where the model implied effect of the chosen predictor is plotted
#' controlling for the other predictors. In addition to the model-implied effect, the relative influence
#' of the predictor is included in the x-axis label. If this is not desired, a custom xlabel can be provided.
#' 
plot.mvtb <- function(x,predictor.no=1,response.no=1,n.trees=NULL,X=NULL,xlab=NULL,ylab=NULL,...){
  if(any(unlist(lapply(x,function(li){is.raw(li)})))){
    x <- mvtb.uncomp(x)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(x$best.trees)) }
  gbm.obj <- x$models[[response.no]]
  ri <- gbm::relative.influence(gbm.obj,n.trees=n.trees)/sum(gbm::relative.influence(gbm.obj,n.trees=n.trees))*100
  ri <- ri[predictor.no]
  #gbm.obj <- convert.mvtb.gbm(out,k=response.no)
  if(is.null(xlab)){ 
    xlab <- paste0(names(ri)," ", formatC(ri,2), "%")
  }
  if(is.null(ylab)) { ylab <- x$ynames[response.no]}
  grid <- gbm::plot.gbm(gbm.obj,i.var = predictor.no,n.trees = n.trees,perspective=TRUE,return.grid=TRUE)  
  plot(y=grid$y,x=grid[,1],type="l",bty="n",xlab=xlab,ylab=ylab,...)
  if(!is.null(X)) { rug(jitter(X[,predictor.no])) }
}


#' Perspective plot for 2 predictors and 1 response. 
#' 
#' This is a plot of the model implied function of 2 predictors averaged over the other predictors
#' included in the model. This is called a partial dependence plot.
#' As an alternative to the perspective (3D) plot, a 2D heat plot can be obtained directly
#' using ?plot.gbm.
#' 
#' @param object mvtb output object
#' @param response.no index of the response variable
#' @param predictor.no vector containing indices of the predictor variables to plot
#' @param n.trees desired number of trees (default: best trees)
#' @param phi angle of viewing direction. See ?persp.
#' @param theta angle of viewing direction See ?persp.
#' @param r distance from eye to center. See ?persp
#' @param d strength of perspective. See ?persp. 
#' @param ticktype 'detailed' gives axis points. See ?persp for other options.
#' @param xlab, title for x axis, must be character strings. 
#' @param ylab, title for y axis, must be character strings.  
#' @param zlab, title for z axis, must be character strings. 
#' @param ... extra arguments are passed to persp. See ?persp
#' @return Function returns a plot.
#' @seealso \code{plot.gbm}, \code{plot.mvtb}, \code{mvtb.heat}
#' @export
#' @importFrom graphics persp
mvtb.perspec <- function(object,response.no=1,predictor.no=1:2,n.trees=NULL,
                         phi=15,theta=-55,r=sqrt(10),d=3,xlab=NULL,ylab=NULL,zlab=NULL,ticktype="detailed",...) {
  if(any(unlist(lapply(object,function(li){is.raw(li)})))){
    object <- mvtb.uncomp(object)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(object$best.trees)) }
  gbm.obj <- object$models[[response.no]]
  grid <- gbm::plot.gbm(gbm.obj,i.var = predictor.no,n.trees = n.trees,perspective=TRUE,return.grid=TRUE)
  x <- as.numeric(unique(grid[,1]))
  y <- as.numeric(unique(grid[,2]))
  if(is.null(xlab)){
      xlab <- colnames(grid)[1]
  }
  if(is.null(ylab)){
      ylab <- colnames(grid)[2]
  }
  if(is.null(zlab)){
    zlab <- object$ynames[response.no]
  }
  z <- matrix(grid[,3],length(unique(x)),length(unique(y)))
  persp(x=x,y=y,z=z,d=d,r=r,phi=phi,theta=theta,ticktype=ticktype,xlab=xlab,ylab=ylab,zlab=zlab,...)
}


# Pairwise plot for 2 predictors and 1 response. 
plot.pw.perspec <- function(out,response.no,predictor.no,npairs=3,nonlin.rank=NULL,p1=NULL,p2=NULL,theta=rep(-55,npairs),...){
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- mvtb.uncomp(out)
  }
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

#' Simple (clustered) heatmap of tables from mvtb (relative influence, covariance explained)
#' 
#' @param x Any table. For example: the covariance explained (\code{res$covex}), or relative influence \code{mvtb.ri(res)}
#' @param clust.method clustering method for rows and columns. See ?hclust. If NULL, unclustered.
#' @param dist.method  method for computing the distance between columns See ?dist for alternatives.
#' @param numformat function to format the covex values into strings. Defaults to removing leading 0 and rounding to 2 decimal places.
#' @param col A list of colors mapping onto covex explained values. A white to black gradient is default.
#' @param cexRow, See \code{cex.axis} from par. The magnification used for the row axis labels. A useful default is provided.
#' @param cexCol, See \code{cex.axis} from par. The magnification used for the col axis labels. The default is set equal to the row axis labels.
#' @param ... extra arguments are passed to image, then to plot. See ?image, ?par
#' @return heatmap of x, usually a covariance explained matrix or a matrix of (relative) influences.
#' @details You will probably want to modify the default colors.
#' @export 
#' @seealso \code{plot.mvtb}, \code{mvtb.perspec}
#' @importFrom graphics image axis text
mvtb.heat <- function(x,clust.method="ward.D",dist.method="manhattan",numformat=NULL,col=NULL,cexRow=NULL,cexCol=NULL,...) {
  if(class(x) %in% "mvtb"){
    if(any(unlist(lapply(x,function(li){is.raw(li)})))){
      x <- mvtb.uncomp(x)
    }
    x <- x$covex
  }
  if(!is.null(clust.method)){
    x <- mvtb.cluster(x,clust.method=clust.method,dist.method=dist.method)
  }
  if(is.null(numformat)){ numformat <- function(val){sub("^(-?)0.", "\\1.", sprintf("%.2f", val))}}
  cellnote <- matrix(numformat(x),dim(x))
  #cellnote <- cellnote[rowInd,colInd] DONT BE TEMPTED TO DO THIS
  x <- t(x)
  cellnote <- t(cellnote)
  nc <- nrow(x) # final number of columns (usually predictors)
  nr <- ncol(x) # final number of rows    (usually dvs)
  if(is.null(col)) { col <- colorRampPaletteAlpha(RColorBrewer::brewer.pal(9,"Greys"),100)}
  image(x=1:nc,y=1:nr,abs(x),xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr),ylab="",xlab="",axes=F,col=col)
  #axis(1,at=seq(0,1,length=nrow(x)))
  #cexRow <- .2+1/log10(max(nc,nr))
  if(is.null(cexRow)) { cexRow <- log10(max(nc,nr)) }
  if(is.null(cexCol)) { cexCol <- cexRow}
  axis(1, 1:nc, labels = rep("",nc), las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  axis(2, 1:nr, labels = colnames(x), las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  text(x =  c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
       col = "white", cex = 1)
  text(1:nc,rep(0,nc), las=2,cex.axis=cexRow,adj=1,
       labels = rownames(x), xpd = TRUE,srt=45,
       col="black")
}
# note that this works for a single predictor, but isn't pretty

#' Add alpha
#'
#' Internal to colorRampPaletteALpha but exported in case a user wants to modify
#'
#' @param colors list of colors
#' @param alpha value of alpha
#' @importFrom grDevices col2rgb rgb 
#' @export
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


#' Add alpha and ramps between colors
#' 
#' Internal to mvtb.heat, but exported for easy modification
#' 
#' @param colors original colors
#' @param n number of colors
#' @param interpolate linear, otherwise spline interpolation is used
#' @importFrom grDevices colorRampPalette col2rgb
#' @importFrom stats approx spline
#' @export
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
  if (interpolate=='linear') {
    l <- approx(a, n=n)
  } else {
    l <- spline(a, n=n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y/255.0)
  return(cr)
}

