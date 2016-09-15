#' boost principal components of outcomes
#' @inheritParams  mvtb 
#' @export
pcboost <- function(Y,X,n.trees=100,
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
  Ystar <- Y %*% ev
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
  class(out) <- c(class(out), "pcboost")
  return(out)
}

#' Predicted values from principal components boosting
#' @inheritParams  predict.mvtb 
#' @export
predict.pcboost <- function(object, n.trees=NULL, newdata, drop=TRUE, ...){
  Yhat <- predict(object, n.trees = n.trees, newdata = newdata, drop = drop) 
  Pred  <- Pred %*% t(object$ev$vectors)
  return(Pred)
}