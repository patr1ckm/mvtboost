#' CV Tune gbm fits over a grid of parameters in parallel
#' 
#' This function cross validates a grid of meta-parameter arguments implied by \code{expand.grid(...)}. 
#' The final model is returned with the best fitting parameters.
#' 
#' @param y vector of outcomes
#' @param x matrix of predictors
#' @param cv.folds number of cross-validation folds
#' @param mc.cores number of cores
#' @param subset subset of the data to use for training; cross-validation will be performed within subset
#' @param ... args to gbm.fit. Arguments can be passed vectors. The function tunes across all rows of \code{expand.grid(...)}
#' @return A list with \code{$gbm} (fitted model with best meta parameters), \code{$best_args} (best meta parameters), and \code{$args} (all meta parameters with cv err)
#' @details Non-grid arguments (such as \code{distribution}, \code{weights}, etc) can be passed as normal.
#' @export
gbm_grid <- function(y, x, cv.folds, mc.cores=1, subset=NULL, ...){
  dots <- list(...)
#  x <- as.matrix(x)
  if(is.null(colnames(x))) colnames(x) <- paste0("X", 1:ncol(x))
  
  n <- length(y)
  if(is.null(subset)){
    train <- 1:n
  } else if(is.logical(subset)){
    train <- which(subset)
  } else {
    train <- subset
  }
  if(all(train == 1:n) & cv.folds==1){ stop("must set subset for training or cv.folds > 1")}
  
  args <- expand.grid(dots, stringsAsFactors = FALSE)
  args_cv <- expand.grid(append(list(k=1:cv.folds), dots), stringsAsFactors = FALSE)
  args.ls <- split(args_cv, 1:nrow(args_cv))
  
  folds <- sample(1:cv.folds, size=length(train), replace=TRUE)
  unique_args <- rep(1:nrow(args), each=cv.folds)
  
  if(mc.cores > 1){
    cl <- parallel::makePSOCKcluster(mc.cores)
    #parallel::clusterExport(cl, gbm::gbm.fit envir = environment())
    parallel::clusterExport(cl, "do_one_fold", envir = environment())
    
    #ocv <- lapply(seq_along(args.ls), FUN=do_one_row, args=args.ls, train=train, folds=folds, Xm=x, y=y)
    ocv <- parallel::clusterApply(cl=cl, x=seq_along(args.ls), 
                        fun=do_one_row, args=args.ls, train=train, folds=folds, Xm=x, y=y)
    
    parallel::stopCluster(cl)
  } else {
    ocv <- lapply(seq_along(args.ls), FUN=do_one_row, args=args.ls, train=train, 
                  folds=folds, Xm=x, y=y)
  }
  
  cv_err <- aggregate_cv_err(ocv, unique_args)
  min_cv_err <- sapply(cv_err, function(e){ min(e) })
  best_arg_idx <- which.min(min_cv_err)
  best_args <- args[best_arg_idx, ]
  args$err <- min_cv_err
  
  d <- data.frame()
  
  out <- do.call(gbm::gbm.fit, append(list(y=y[train], x=x[train, ,drop=F]), best_args))
  out$cv.error <- cv_err[[best_arg_idx]]
  # todo: out$valid.error
  return(list(gbm=out, best_args=best_args, args=args))
}

aggregate_cv_err <- function(ocv, unique_args){
  tapply(ocv, unique_args, FUN=function(l){
    rowMeans(data.frame(l))
  }, simplify=FALSE)
}

do_one_fold <- function(k, folds, train, y, x, ...){
  s <- train[folds != k]
  if(length(s) == 0){ s <- train}
  o <- gbm::gbm.fit(y=y[s], x=x[s, ,drop=F], ...)
  yhat <- as.matrix(predict(o, newdata=data.frame(x[-s,,drop=F]), n.trees=1:o$n.trees))
  test_err <- apply(yhat, 2, function(yh, y){mean((yh-y)^2)}, y=y[-s])
  return(test_err)
}

do_one_row <- function(i, args, train, folds, Xm, y){
  do.call(do_one_fold, append(args[[i]], list(x=Xm, y=y, train=train, folds=folds))) 
}
