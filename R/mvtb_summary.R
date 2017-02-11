


#' @importFrom stats var
mvtb.r2 <- function(object,Y,X,n.trees=NULL){
  
  if(is.null(n.trees)) { n.trees <- object$best.iter[[2]] }
  p <- predict.mvtb(object,n.trees,newdata=X)
  1-apply(Y - p,2,var)/apply(Y,2,var)
}

#' Simple default printing of the mvtb output object
#' @param x mvtb output object
#' @param ... unused
#' @export
#' @importFrom utils str
print.mvtb <- function(x,...) {
  #if(any(unlist(lapply(x,function(li){is.raw(li)})))){
  #  x <- mvtb.uncomp(x)
  #  cat("COMPRESSED",fill=TRUE)
  #}
  str(x,1)
}

#' Computes a summary of the multivariate tree boosting model
#' 
#' @param object mvtb output object
#' @param print result (default is TRUE)
#' @param n.trees number of trees to use. Defaults to the minimum number of trees by CV, test, or training error for each outcome.
#' @param relative relative If 'col', each column sums to 100. If 'tot', the whole matrix sums to 100 (a percent). If 'n', the raw reductions in SSE are returned.
#' @param ... additional arguments affecting the summary produced.
#' @return Returns the best number of trees, the univariate relative influence of each predictor for each outcome, and covariance explained in pairs of outcomes by each predictor
#' @seealso \code{mvtb.ri}, \code{gbm.ri}, \code{mvtb.cluster}
#' @export
summary.mvtb <- function(object, print = TRUE, n.trees = NULL, relative = "col", ...) {
  
  if(any(unlist(lapply(object,function(li){is.raw(li)})))){
    object <- mvtb.uncomp(object)
  }
  k <- length(object$models)
  if(is.null(n.trees)) { n.trees <- apply(object$best.trees, 1, min, na.rm=T) }
  
  if(length(n.trees) == 1){ n.trees <- rep(n.trees, k)}
  
  ri <- influence.mvtb(object, n.trees = n.trees, relative = relative)
  
  sum <- list(best.trees = n.trees, relative.influence = ri)
  
  if(print){ print(lapply(sum,function(o){round(o,2)})) }
  invisible(sum)
}

#' Clustering the covariance explained or relative influence matrix
#' 
#' The 'covariance explained' by each predictor is the reduction in covariance between each pair of outcomes due to splitting on each predictor over all trees (\code{$covex}).
#' To aid in the interpretability of the covariance explained matrix, this function clusters the rows (pairs of outcomes) and the columns (predictors) of \code{object$covex}
#' so that groups of predictors that explain similar pairs of covariances are closer together.
#' This function can also be used to cluster the relative influence matrix. In this case, the rows (usually outcomes) and columns (usually predictors) with similar values will
#' be clustered together.
#'  
#' @param x Any matrix, such as \code{mvtb.covex(object)}, or \code{mvtb.ri(object)}. 
#' @param clust.method clustering method for rows and columns. This should be (an unambiguous abbreviation of) one of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"} (= UPGMA), \code{"mcquitty"} (= WPGMA), \code{"median"} (= WPGMC) or \code{"centroid"} (= UPGMC).
#' @param dist.method  method for computing the distance between two lower triangular covariance matrices. This must be one of \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} or \code{"minkowski"}. Any unambiguous substring can be given.
#' @param plot Produces a heatmap of the covariance explained matrix. see \code{?mvtb.heat}
#' @param ... Arguments passed to \code{mvtb.heat} 
#' @return clustered covariance matrix, with re-ordered rows and columns.
#' @seealso \code{mvtb.heat}
#' @export
#' @details The covariance explained by each predictor is only unambiguous if the predictors are uncorrelated and interaction.depth = 1. 
#' If predictors are not independent, the decomposition of covariance explained is only approximate (like the decomposition of R^2 by each predictor in a linear model). 
#' If interaction.depth > 1, the following heuristic is used: the covariance explained by the tree is assigned to the predictor with the largest influence in each tree.
#'
#' Note that different distances measures (e.g. \code{"manhattan"}, \code{"euclidean"}) provide different ways to measure (dis)similarities between 
#' the covariance explained patterns for each predictor. See \code{?dist} for further details.
#' After the distances have been computed, \code{hclust} is used to form clusters. 
#' Different clustering methods (e.g. \code{"ward.D"}, \code{"complete"}) generally group rows and columns differently (see \code{?hclust} for further details).
#' It is suggested to try different distance measures and clustering methods to obtain the most interpretable solution. 
#' The defaults are for \code{"euclidean"} distances and \code{"complete"} clustering.
#' Transposing the rows and columns may also lead to different results.
#'
#' A simple heatmap of the clustered matrix can be obtained by setting \code{plot=TRUE}. Details of the plotting procedure are available via \code{mvtb.heat}.
#' 
#' \code{covex} values smaller than \code{getOption("digits")} are truncated to 0. Note that it is possible to obtain negative variance explained
#' due to sampling fluctuation. These can be truncated or ignored. 
#'
#' @importFrom stats hclust dist as.dendrogram order.dendrogram
mvtb.cluster <- function(x, clust.method = "complete", dist.method = "euclidean", plot=FALSE, ...) {
    if(nrow(x) > 1) { 
      hcr <- hclust(dist(x, method=dist.method), method=clust.method)
      ddr <- as.dendrogram(hcr)
      rowInd <- order.dendrogram(ddr)
    } else {
        rowInd <- 1
    }
    if(nrow(t(x)) > 1) {
      hcc <- hclust(dist(t(x), method=dist.method), method=clust.method)
      ddc <- as.dendrogram(hcc)
      colInd <- order.dendrogram(ddc)
    } else {
      colInd <- 1
    }
    x <- x[rowInd,colInd,drop=FALSE]
    x <- zapsmall(x)
    if(plot){
      mvtb.heat(x, clust.method=clust.method, dist.method=dist.method,...)
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




