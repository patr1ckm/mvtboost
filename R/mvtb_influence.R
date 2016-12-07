#' Computes the relative influence of each predictor for each outcome
#' 
#' The relative influence of a predictor is the reduction in sums of squares attributable to splits on individual predictors.
#' It is often expressed as a percent (sums to 100).
#' @param model \code{mvtb} output model
#' @param n.trees number of trees to use. Defaults to the minimum number of trees by CV, test, or training error for each outcome.
#' @param relative How to scale the multivariate influences. If \code{"col"}, each column sums to 100. If \code{"tot"}, the whole matrix sums to 100 (a percent). Otherwise, the raw reductions in SSE are returned.
#' @param ... Additional arguments passed to \code{gbm::relative.influence}
#' @return Matrix of (relative) influences.
#' @export 
mvtb.ri <- function(model, n.trees=NULL, relative="col",...){
  
  if(any(unlist(lapply(model,function(li){is.raw(li)})))){
    model <- mvtb.uncomp(model)
  }
  k <- length(model$models)
  if(is.null(n.trees)) { n.trees <- apply(model$best.trees, 1, min, na.rm=T) }
  if(length(n.trees) == 1){ n.trees <- rep(n.trees, k)}
  
  ri <- matrix(0,nrow=length(model$xnames),ncol=k)
  for(i in 1:k) {
    gbm.obj <- model$models[[i]]
    ri[,i] <- gbm::relative.influence(gbm.obj,n.trees=n.trees[i],...)
  }
  if(relative == "col"){
    ri <- matrix(apply(ri,2,function(col){col/sum(col)})*100,nrow=nrow(ri),ncol=ncol(ri))
  } else if (relative=="tot") {
    ri <- ri/sum(ri)*100
  }
  colnames(ri) <- model$ynames
  rownames(ri) <- model$xnames
  return(ri)  
}


#' @importFrom stats influence
NULL


#' Compute influence scores from mvtb
#' @inheritParams mvtb.ri
#' @export
influence.mvtb <- mvtb.ri

#' @export
influence.twostage <- function(model, n.trees = NULL, relative = TRUE, sort = FALSE, ...){
  if(is.null(n.trees)){ n.trees = model$tr}
  inf <- gbm::relative.influence(model$o.gbm, n.trees = n.trees, scale. = FALSE, sort. = sort)
  if(relative){
    inf <- inf / sum(inf) * 100
  }
  return(inf)
}

#' @export
influence.lmerboost <- function(model, n.trees = NULL, relative = TRUE, sort = FALSE, ...){
  if(is.null(n.trees)){ 
    n.trees <- min(model$best.trees, na.rm = TRUE)
  }
  inf <- influence_from_tree_list(model$trees, n.trees = n.trees, 
                                             var.names = model$xnames)
  inf <- inf[-model$id]
  if(relative) { inf <- (inf / sum(inf)) * 100 }
  if(sort) { inf <- sort(inf, decreasing = TRUE) }
  inf
}


#' Compute the influnce from pcb
#' @export
influence.pcb <- function(model, n.trees = NULL, relative = "col", ...){
  ri <- influence.mvtb(model = model, n.trees=n.trees, relative = FALSE) 
  ri <- data.frame(ri %*% t(model$ev$vectors))
  
  if(relative == "col"){
    ri <- matrix(apply(ri,2,function(col){col/sum(col)})*100, nrow=nrow(ri), ncol=ncol(ri))
  } else if (relative=="tot") {
    ri <- ri/sum(ri)*100
  } # else do nothing
  
  rownames(ri) <- model$xnames
  colnames(ri) <- model$ynames
  return(ri)
}

influence_from_tree_list <- function(model, n.trees=1, var.names) {
  get.rel.inf <- function(obj) {
    lapply(split(obj[[6]], obj[[1]]), sum)
  }
  temp <- unlist(lapply(model[1:n.trees], get.rel.inf))
  rel.inf.compact <- unlist(lapply(split(temp, names(temp)), 
                                   sum))
  rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) != 
                                       "-1"]
  rel.inf <- rep(0, length(var.names))
  i <- as.numeric(names(rel.inf.compact)) + 1
  rel.inf[i] <- rel.inf.compact
  names(rel.inf) <- var.names 
  
  return(rel.inf = rel.inf)
}
