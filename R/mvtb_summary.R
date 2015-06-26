
#gbm.ri <- function(out,n.trees=out$best.iter[[2]]){
#  k <- length(out05$finaltree)
#  ri <- matrix(0,nrow=nrow(out05$ri[[1]]),ncol=ncol(out05$ri[[1]]),dimnames=dimnames(out05$ri[[1]]))
#  for(i in 1:k) {
#    gbm.obj <- convert.mvtb.gbm(out,k=i)
#    ri[,i] <- relative.influence(gbm.obj,n.trees=n.trees)
#  }
#  return(ri)  
#}

#' Computes the relative influence of each predictor for each outcome.
#' 
#' The relative influence of a predictor is the reduction in sums of squares attributable to splits on individual predictors.
#' It is often expressed as a percent (sums to 100).
#' @param out mvtb output object
#' @param n.trees number of trees to use
#' @param weighted T/F. Reductions in SSE are weighted according the covariance explained by each predictor.
#' @param relative If 'col', each column sums to 100. If 'tot', the whole matrix sums to 100 (a percent). If 'n', the raw reductions in SSE are returned.
#' @return Matrix of (relative) influences.
#' @export 
mvtb.ri <- function(out,n.trees=NULL,weighted=F,relative="col"){
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- uncomp.mvtb(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  if(weighted) {
    ri <- apply(out$w.rel.infl[,,1:n.trees,drop=FALSE],1:2,sum)
  } else {
    ri <- apply(out$rel.infl[,,1:n.trees,drop=FALSE],1:2,sum)
  }
  if(relative == "col"){
    ri <- apply(ri,2,function(col){col/sum(col)})*100
  } else if (relative=="tot") {
    ri <- ri/sum(ri)*100
  }
  colnames(ri) <- out$ynames
  rownames(ri) <- out$xnames
  return(ri)
}

r2 <- function(out,Y,X,n.trees=NULL){
  if(is.null(n.trees)) { n.trees <- out$best.iter[[2]] }
  p <- predict.mvtb(out,n.trees,newdata=X)
  1-apply(Y - p,2,var)/apply(Y,2,var)
}

#' Simple default printing of the mvtb output object
#' @param x mvtb output object
#' @param ... unused
#' @export
print.mvtb <- function(x,...) {
  if(any(unlist(lapply(x,function(li){is.raw(li)})))){
    x <- uncomp.mvtb(x)
  }
  str(x,1)
}

#' Computes a summary of the multivariate tree boosting model
#' 
#' @param object mvtb output object
#' @param print result (default is TRUE)
#' @param n.trees number of trees used to compute relative influence. Defaults to the minimum number of trees by CV, test, or training error
#' @param ... unused
#' @return Returns the best number of trees, the univariate relative influence of each predictor for each outcome, and covariance explained in pairs of outcomes by each predictor
#' @seealso \code{mvtb.ri}, \code{gbm.ri}, \code{cluster.covex}
#' @export
summary.mvtb <- function(object,print=TRUE,n.trees=NULL,...) {
  out <- object
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- uncomp.mvtb(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  ri <- mvtb.ri(out,n.trees=n.trees)
  cc <- cluster.covex(out)
  sum <- list(best.trees=n.trees,relative.influence=ri,cluster.covex=cc)
  if(print){ print(lapply(sum,function(o){round(o,2)})) }
  invisible(sum)
}

#' Computing a clustered covariance explained matrix
#' 
#' For each pair of predictors, computes the distance between the correlation matrices of the outcomes explained by those predictors.
#'  
#' @param out mvtb output
#' @param clust.method clustering method for rows and columns. See \code{?hclust}
#' @param dist.method  method for computing the distance between two lower triangluar covariance matrices. See \code{?dist} for alternatives.
#' @return clustered covariance matrix, with rows and columns.
#' @seealso \code{heat.covex}
#' @export
cluster.covex <- function(out,clust.method="ward.D",dist.method="manhattan") {
    if(any(unlist(lapply(out,function(li){is.raw(li)})))){
      out <- uncomp.mvtb(out)
    }
    x <- out$covex
    hcr <- hclust(dist(x,method=dist.method),method=clust.method)
    ddr <- as.dendrogram(hcr)
    rowInd <- order.dendrogram(ddr)
    hcc <- hclust(dist(t(x),method=dist.method),method=clust.method)
    ddc <- as.dendrogram(hcc)
    colInd <- order.dendrogram(ddc)
    x <- x[rowInd,colInd]
    return(x)
}

#' Uncompress a compressed mvtb output object
#' @param out an object of class \code{mvtb}
#' @export
uncomp.mvtb <- function(out) { 
  o <- lapply(out,function(li){unserialize(memDecompress(li,type="bzip2"))})
  class(o) <- "mvtb"
  return(o)
}
