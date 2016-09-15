#' boost principal components of outcomes
#' @inheritParams  mvtb 
#' @export
pcb <- function(Y,X,n.trees=100,
                    shrinkage=.01,
                    interaction.depth=1,
                    distribution="gaussian",
                    train.fraction=1,
                    bag.fraction=1,
                    cv.folds=1,
                    keep.data=FALSE,
                    s=NULL,
                    seednum=NULL,
                    compress=FALSE,
                    save.cv=FALSE,
                    iter.details=TRUE,
                    verbose=FALSE, 
                    mc.cores=1, ...){
  
  ev <- eigen(cov(Y))$values
  Ystar <- Y %*% ev$vectors
  out <- mvtb.sep(Y=Ystar, X=X, n.trees=n.trees,
           shrinkage=shrinkage,
           interaction.depth = interaction.depth,
           train.fraction = train.fraction,
           bag.fraction = bag.fraction,
           cv.folds = cv.folds,
           keep.data = keep.data,
           s = s,
           seednum = seednum,
           compress = compress,
           save.cv = save.cv,
           iter.details = iter.details,
           verbose = verbose,
           mc.cores = mc.cores)
  out$ev <- ev
  class(out) <- c(class(out), "pcb")
  return(out)
}

#' Predicted values from principal components boosting
#' @inheritParams  predict.mvtb 
#' @export
predict.pcb <- function(object, n.trees = NULL, newdata, drop=TRUE, ...){
  Yhat <- predict(object, n.trees = n.trees, newdata = newdata, drop = drop) 
  Pred  <- Pred %*% t(object$ev$vectors)
  return(Pred)
}


#' Compute the influnce from pcb
#' @export
influence.pcb <- function(object, n.trees = NULL, relative = "col", ...){
  ri <- influence.mvtb(object = object, n.trees=n.trees, relative = FALSE) 
  ri %*% t(object$ev$vectors)
  if(relative == "col"){
    ri <- matrix(apply(ri,2,function(col){col/sum(col)})*100,nrow=nrow(ri),ncol=ncol(ri))
  } else if (relative=="tot") {
    ri <- ri/sum(ri)*100
  } # else do nothing
  return(ri)
}


