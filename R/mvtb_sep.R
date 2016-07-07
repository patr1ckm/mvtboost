#' Separate tuning of multivariate tree boosting models. 
#' @export
mvtb.sep <- function(Y,X,n.trees=100,
                 shrinkage=.01,
                 interaction.depth=1,
                 distribution="gaussian",
                 train.fraction=1,
                 bag.fraction=1,
                 cv.folds=1,
                 s=NULL,
                 seednum=NULL,
                 compress=FALSE,
                 mc.cores=1,
                 save.cv=FALSE,
                 iter.details=TRUE,
                 verbose=FALSE, ...) {
  
  Y <- as.data.frame(Y)
  X <- as.data.frame(X)
  n <- nTrain <- nrow(X)
  
  if(!is.null(seednum)) set.seed(seednum)
  if(train.fraction < 1) nTrain <- ceiling(train.fraction*n) 
  
  if(is.null(s)) { 
    s <- 1:n 
    snull <- TRUE
  } else {
    snull <- FALSE
    nTrain <- NULL
  }

  do.one <- function(y, x, s, ...){ gbm::gbm.fit(y=y[s], x=x[s,], ...) }
  
  get.pred.err <- function(o, y, x, n.trees){
    yhat.iter <- predict(o, newdata=x, n.trees=1:n.trees) 
    apply(yhat.iter, 2, function(yhat, y){
      var(y - yhat)
    }, y=y)
  }
  get.test.err <- function(oL, Y, x, s, n.trees){
    test.err <- vector(mode = "list", length=ncol(Y))
    names(test.err) <- colnames(Y)
    for(i in 1:ncol(Y)){
      test.err[[i]] <- get.pred.err(o=oL[[i]], y=Y[-s, i], x=x[-s, ], n.trees=n.trees)
    }
    return(test.err)
  }
  
  cv.mods <- cv.mod.err <- vector(mode="list", length(ncol(Y)))
  
  if(cv.folds > 1){
    folds <- sample(1:cv.folds, size=length(s), replace=T)
    for(i in 1:cv.folds){
      train <- s[folds != i]
      test <- s[folds == i]
      cv.mods[[i]] <- parallel::mclapply(Y, FUN=do.one, x=X, n.trees=n.trees, shrinkage=shrinkage, 
                             interaction.depth=interaction.depth, distribution=distribution,
                             nTrain=nTrain, bag.fraction=bag.fraction, s=s, 
                             keep.data=FALSE, verbose=verbose, mc.cores=mc.cores,...)
      cv.mod.err[[i]] <- get.test.err(mods, Y=Y, x=X, n.trees=n.trees, s=test)
    }
  }
  
  mods <- parallel::mclapply(Y, FUN=do.one, x=X, n.trees=n.trees, shrinkage=shrinkage, 
                 interaction.depth=interaction.depth, distribution=distribution,
                 nTrain=nTrain, bag.fraction=bag.fraction, s=s, 
                 keep.data=FALSE, verbose=verbose, mc.cores=mc.cores,...)
  
  train.err <- lapply(mods, function(o){ o$train.error})
  oob.err <- lapply(mods, function(o){ -cumsum(o$oobag.improve)})
  if(cv.folds > 1) {
    cv.err <- lapply(1:ncol(Y), function(i) { rowMeans(sapply(cv.mod.err, "[[", i)) })
  } else {
    cv.err <- lapply(1:ncol(Y), function(i){ rep(NaN, n.trees) })
  }
  names(cv.err) <- colnames(Y)
  if(snull){
    test.err <- lapply(mods, function(o){ o$valid.error})
  } else {
    test.err <- get.test.err(mods, Y=Y, x=X, s=s, n.trees=n.trees)
  }
  
  # which.min can return length 0, return NA instead
  which.min.na <- function(x){
    res <- which.min(x)
    if(length(res) == 0){
      res <- NA
    }
    res
  }
  
  best.trees <- data.frame(train = sapply(train.err, which.min.na), 
                     test = sapply(test.err, which.min.na),
                     oob = sapply(oob.err, which.min.na), 
                     cv = sapply(cv.err, which.min.na))
  rownames(best.trees) <- colnames(Y)
  fl <- list(mods=mods, best.trees=best.trees, 
              train.err=train.err, oob.err=oob.err, cv.err=cv.err, test.err=test.err,
              s=s,n=nrow(X),xnames=colnames(X),ynames=colnames(Y))
  class(fl) <- "mvtb"
  return(fl)
   
}

