#' Computes the relative influence of each predictor for each outcome from gbm.
#' 
#' The relative influence of a predictor is the reduction in sums of squares attributable to splits on individual predictors.
#' It is often expressed as a percent (sums to 100).
#' @param object mvtb output object
#' @param n.trees number of trees to use
#' @param relative How to scale the multivariate influences. If 'col', each column sums to 100. If 'tot', the whole matrix sums to 100 (a percent). If 'n', the raw reductions in SSE are returned.
#' @param ... Additional arguments passed to \code{gbm::relative.influence}
#' @return Matrix of (relative) influences.
#' @export 
gbm.ri <- function(object,n.trees=NULL,relative="col",...){
  out <- object
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- mvtb.uncomp(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  k <- length(out$models)
  ri <- matrix(0,nrow=length(out$xnames),ncol=k)
  for(i in 1:k) {
    gbm.obj <- out$models[[i]]
    ri[,i] <- gbm::relative.influence(gbm.obj,n.trees=n.trees,...)
  }
  if(relative == "col"){
    ri <- matrix(apply(ri,2,function(col){col/sum(col)})*100,nrow=nrow(ri),ncol=ncol(ri))
  } else if (relative=="tot") {
    ri <- ri/sum(ri)*100
  }
  colnames(ri) <- out$ynames
  rownames(ri) <- out$xnames
  return(ri)  
}

#' Computes the relative influence of each predictor for each outcome.
#' 
#' The relative influence of a predictor is the reduction in sums of squares attributable to splits on individual predictors.
#' It is often expressed as a percent (sums to 100).
#' @param object mvtb output object
#' @param n.trees number of trees to use
#' @param weighted T/F. Reductions in SSE are weighted according the covariance explained by each predictor. See \code{?mvtb} for details.
#' @param relative If 'col', each column sums to 100. If 'tot', the whole matrix sums to 100 (a percent). Otherwise, the raw reductions in SSE are returned.
#' @return Matrix of (relative) influences.
#' @export 
mvtb.ri <- function(object,n.trees=NULL,weighted=F,relative="col"){
  out <- object
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- mvtb.uncomp(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  if(weighted) {
    ri <- apply(out$w.rel.infl[,,1:n.trees,drop=FALSE],1:2,sum)
  } else {
    ri <- apply(out$rel.infl[,,1:n.trees,drop=FALSE],1:2,sum)
  }
  if(relative == "col"){
    ri <- matrix(apply(ri,2,function(col){col/sum(col)})*100,nrow=nrow(ri),ncol=ncol(ri))
  } else if (relative=="tot") {
    ri <- ri/sum(ri)*100
  }
  colnames(ri) <- out$ynames
  rownames(ri) <- out$xnames
  return(ri)
}

#' @importFrom stats var
r2 <- function(object,Y,X,n.trees=NULL){
  out <- object
  if(is.null(n.trees)) { n.trees <- out$best.iter[[2]] }
  p <- predict.mvtb(out,n.trees,newdata=X)
  1-apply(Y - p,2,var)/apply(Y,2,var)
}

#' Simple default printing of the mvtb output object
#' @param x mvtb output object
#' @param ... unused
#' @export
#' @importFrom utils str
print.mvtb <- function(x,...) {
  if(any(unlist(lapply(x,function(li){is.raw(li)})))){
    x <- mvtb.uncomp(x)
  }
  str(x,1)
}

#' Computes a summary of the multivariate tree boosting model
#' 
#' @param object mvtb output object
#' @param print result (default is TRUE)
#' @param n.trees number of trees used to compute relative influence. Defaults to the minimum number of trees by CV, test, or training error
#' @param relative relative If 'col', each column sums to 100. If 'tot', the whole matrix sums to 100 (a percent). If 'n', the raw reductions in SSE are returned.
#' @param ... unused
#' @return Returns the best number of trees, the univariate relative influence of each predictor for each outcome, and covariance explained in pairs of outcomes by each predictor
#' @seealso \code{mvtb.ri}, \code{gbm.ri}, \code{mvtb.cluster}
#' @export
summary.mvtb <- function(object,print=TRUE,n.trees=NULL,relative="tot",...) {
  out <- object
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- mvtb.uncomp(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  ri <- mvtb.ri(out,n.trees=n.trees,relative=relative)
  cc <- mvtb.cluster(out)
  sum <- list(best.trees=n.trees,relative.influence=ri,mvtb.cluster=cc)
  if(print){ print(lapply(sum,function(o){round(o,2)})) }
  invisible(sum)
}

#' Clustering the covariance explained or relative influence matrix
#' 
#' The 'covariance explained' by each predictor is the reduction in covariance between each pair of outcomes due to splitting on each predictor over all trees (\code{$covex}).
#' To aid in the interpretability of the covariance explained matrix, this function clusters the rows (pairs of outcomes) and the columns (predictors) of \code{object$covex}
#' so that groups of predictors that explain similar pairs of covariances are closer together.
#' 
#' This function can also be used to cluster the relative influence matrix. In this case, the rows (usually outcomes) and columns (usually predictors) with similar values will
#' be clustered together.
#'  
#' @param Any table (e.g. \code{object$covex}), \code{mvtb.ri(object)}. If an mvtb object, defaults to \code{object$covex}
#' @param clust.method clustering method for rows and columns. See \code{?hclust}
#' @param dist.method  method for computing the distance between two lower triangluar covariance matrices. See \code{?dist} for alternatives.
#' @param plot Produces a heatmap of the covariance explained matrix. see \code{?mvtb.heat}
#' @param ... Arguments passed to \code{mvtb.heat} 
#' @return clustered covariance matrix, with re-ordered rows and columns.
#' @seealso \code{mvtb.heat}
#' @export
#' @details The covariance explained by each predictor is only unambiguous if the predictors are uncorrelated and interaction.depth = 1. 
#' If predictors are not independent, the decomposition of covariance explained is only approximate (like the decomposition of R^2 by each predictor in a linear model). 
#' If interaction.depth > 1, the following heuristic is used: the covariance explained by the tree is assigned to the predictor with the largest influence in each tree.
#'
#' Note that different distances measures (e.g. 'manhattan', 'euclidean') provide different ways to measure distances between 
#' the covariance explained patterns for each predictor. See \code{?dist} for further details.
#' After the distances have been computed, \code{hclust} is used to form clusters. 
#' Different clustering methods (e.g. "ward.D", "complete") generally group rows and columns differently.
#' It is suggested to try different distance measures and clustering methods to obtain the most interpretable solution. 
#' The defaults are for 'manhattan' distances and 'ward.D' clustering, which seem to provide reasonable solutions in many cases.
#' Transposing the rows and columns may also lead to different results.
#'
#' A simple heatmap of the clustered matrix can be obtained by setting \code{plot=TRUE}. Details of the plotting procedure are available via \code{mvtb.heat}.
#' 
#' \code{covex} values smaller than \code{getOption("digits")} are truncated to 0. Note that it is not impossible to obtain negative variance explained
#' due to sampling fluctuaion. These can be truncated or ignored. 
#'
#' @importFrom stats hclust dist as.dendrogram order.dendrogram
mvtb.cluster <- function(x,clust.method="ward.D",dist.method="manhattan",plot=FALSE,...) {
    out <- x
    if(class(out) %in% "mvtb"){
      if(any(unlist(lapply(out,function(li){is.raw(li)})))){
        out <- mvtb.uncomp(out)
      }
      x <- out$covex
    }
    if(nrow(x) > 1) { 
      hcr <- hclust(dist(x,method=dist.method),method=clust.method)
      ddr <- as.dendrogram(hcr)
      rowInd <- order.dendrogram(ddr)
    } else {
        rowInd <- 1
    }
    if(nrow(t(x)) > 1) {
      hcc <- hclust(dist(t(x),method=dist.method),method=clust.method)
      ddc <- as.dendrogram(hcc)
      colInd <- order.dendrogram(ddc)
    } else {
      colInd <- 1
    }
    x <- x[rowInd,colInd,drop=FALSE]
    x <- zapsmall(x)
    if(plot){
      mvtb.heat(out,clust.method=clust.method,dist.method=dist.method,...)
    }
    return(x)
}

#' Uncompress a compressed mvtb output object
#' 
#' This function uncompresses a compressed mvtb output object. All elements are uncompressed.
#' @param object an object of class \code{mvtb}
#' @export
mvtb.uncomp <- function(object) { 
  o <- lapply(object,function(li){unserialize(memDecompress(li,type="bzip2"))})
  class(o) <- "mvtb"
  return(o)
}


