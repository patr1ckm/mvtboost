
mvtb_grid <- function(Y,X,n.trees=100,
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
  
  n <- nrow(X)
  if (is.null(s)) {                                                        
    train <- sample(1:n, floor(n * train.fraction), replace = F)             
  } else {
    train <- s
  }
  stopifnot(distribution=="guassian", "Distribution must be gaussian for now")
  grid <- expand.grid(n.trees=n.trees, 
                      shrinkage=shrinkage, 
                      interaction.depth=interaction.depth)
  grid_cv <- expand.grid(k=1:cv.folds, 
                         n.trees=n.trees, 
                         shrinkage=shrinkage, 
                         interaction.depth=interaction.depth)
                      
                      
  folds <- sample(1:cv.folds, size=length(train), replace=T)
  grid_ls <- split(grid_cv, 1:nrow(grid))
  
  do_one <- function(args, folds, train, x, ...){
    new_args <- append(args, list(s=train[folds != args$k], X=x, ...))
    new_args$k <- NULL
    do.call(mvtb_sep, new_args)
  }
  mods <- parallel::mclapply(grid_ls, FUN=do_one, folds=folds, train=train, 
                             Y=Y, x=X, distribution=distribution,
                             bag.fraction=bag.fraction, 
                             keep.data=keep.data,
                             seednum=seednum,
                             compress=compress, save.cv=save.cv,
                             iter.details=iter.details,
                             verbose=verbose,
                             mc.cores=mc.cores)
  unique_args <- factor(rep(1:nrow(grid), each=cv.folds))
  # Choose best args by averaging over outcomes
  cv_err <- lapply(mods, function(x){rowMeans(data.frame(x$test.err))})
  # then by meta parameter
  arg_cv_err <- tapply(cv_err, unique_args, function(x){min(rowMeans(data.frame(x)))})
  
  best <- which.min(arg_cv_err)
  best_args <- grid[best, ]
  best_args$k <- cv.folds + 1
  
  # Fit a model using best args on full training sample
  best_mod <- do_one(best_args, folds, train=train, Y=Y, x=X, 
                     distribution=distribution,
                     bag.fraction=bag.fraction, 
                     keep.data=keep.data,
                     seednum=seednum,
                     compress=compress, save.cv=save.cv,
                     iter.details=iter.details,
                     verbose=verbose,
                     cv.folds=cv.folds,
                     mc.cores=mc.cores)
  
  best_args <- data.frame(grid[best, ], err=arg_cv_err[best])
  args <- data.frame(grid, err=arg_cv_err)
  fl <- append(best_mod, list(best_args=best_args, args=args))
  class(fl) <- c("mvtb")
  return(fl)
}

which.min.na <- function(x){
  res <- which.min(x)
  if(length(res) == 0){
    res <- NA
  }
  res
}