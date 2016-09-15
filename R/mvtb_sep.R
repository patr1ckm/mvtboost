#' mvtb tuning each outcome separately.
#' @inheritParams  mvtb 
#' @export
mvtb_sep <- function(Y,X,n.trees=100,
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
                 mc.cores=1, ...) {
  # Gives Y and X different automatic names. Y and X need to be data frames b/c loops over columns.
  if(!is.data.frame(Y)) Y <- data.frame(Y)  
  if(!is.data.frame(X)) X <- as.data.frame(X)
  n <- nrow(X)
  
  params <- c(as.list(environment()),list(...)) # this won't copy y and x
  ## Checks
  if(shrinkage > 1 | shrinkage <= 0){ stop("shrinkage should be > 0, < 1")}
  if(train.fraction > 1 | train.fraction <= 0){ stop("train.fraction should be > 0, < 1")}
  if(bag.fraction > 1 | bag.fraction <= 0){ stop("bag.fraction should be > 0, < 1")}

  
  if(!is.null(seednum)) set.seed(seednum)
  if(is.null(s)) { 
    s <- sample(1:n, floor(n*train.fraction), replace=F) #force round down if odd
  }

  do.one <- function(y, x, s, ...){ gbm::gbm.fit(y=y[s], x=as.data.frame(x[s,]), ...) }
  
  get.pred.err <- function(o, y, x, n.trees){
    yhat.iter <- as.matrix(gbm::predict.gbm(o, newdata=x, n.trees=1:n.trees))
    apply(yhat.iter, 2, function(yhat, y){
      var(y - yhat)
    }, y=y)
  }
  get.test.err <- function(oL, Y, x, s, n.trees){
    test.err <- vector(mode = "list", length=ncol(Y))
    names(test.err) <- colnames(Y)
    for(i in 1:ncol(Y)){
      test.err[[i]] <- get.pred.err(o=oL[[i]], y=Y[-s, i], x=x[-s,], n.trees=n.trees)
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
                             bag.fraction=bag.fraction, s=s, 
                             keep.data=keep.data, verbose=verbose, mc.cores=mc.cores,...)
      cv.mod.err[[i]] <- get.test.err(cv.mods[[i]], Y=Y, x=X, n.trees=n.trees, s=test)
    }
  }
  
  models <- parallel::mclapply(Y, FUN=do.one, x=X, n.trees=n.trees, shrinkage=shrinkage, 
                 interaction.depth=interaction.depth, distribution=distribution,
                 bag.fraction=bag.fraction, s=s, 
                 keep.data=keep.data, verbose=verbose, mc.cores=mc.cores,...)
  
  train.err <- lapply(models, function(o){ o$train.error})
  oob.err <- lapply(models, function(o){ -cumsum(o$oobag.improve)})
  if(cv.folds > 1) {
    cv.err <- lapply(1:ncol(Y), function(i) { rowMeans(sapply(cv.mod.err, "[[", i)) })
  } else {
    cv.err <- lapply(1:ncol(Y), function(i){ rep(NaN, n.trees) })
  }
  names(cv.err) <- colnames(Y)
  if(all(1:n %in% s)){
    test.err <- lapply(models, function(m){ m$valid.error})
  } else {
    test.err <- get.test.err(models, Y=Y, x=X, s=s, n.trees=n.trees)
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
  
  if(!save.cv){cv.mods <- NULL}
  if(!iter.details){train.err <- NULL; test.err <- NULL; cv.err = NULL}

  fl <- list(models=models, best.trees=best.trees, params=params,
             train.err=train.err, test.err=test.err, cv.err=cv.err,
             cv.mods=cv.mods,
             s=s,n=nrow(X), xnames=colnames(X), ynames=colnames(Y))

  if(compress) {
    # compress each element using bzip2
    fl <- lapply(fl,comp)
  }

  class(fl) <- "mvtb"
  return(fl)
}

