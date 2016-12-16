# Testing gbm_grid
context("test_gbm_grid")

n <- 1000
x <- matrix(rnorm(n))
colnames(x) <- "x"
y <- x*.5
n <- length(y)
cv.folds <- 3
train <- seq_along(y)
folds <- sample(1:cv.folds, size=length(train), replace=TRUE)
tol = 1E-4

args_cv <- expand.grid(
  k=1:3,
  interaction.depth=1:3, 
  shrinkage=c(.5,1),
  n.minobsinnode=c(5,10), 
  bag.fraction=c(.5, 1),
  distribution="gaussian", 
  verbose=FALSE, 
  stringsAsFactors = FALSE)

args <- expand.grid(
  interaction.depth=1:3, 
  shrinkage=c(.5,1),
  n.minobsinnode=c(5,10), 
  bag.fraction=c(.5, 1),
  distribution="gaussian", 
  verbose=FALSE, 
  stringsAsFactors = FALSE)

args.ls <- split(args_cv, 1:nrow(args_cv))

dots <- list(interaction.depth=1:3, n.trees=500, distribution="gaussian", verbose=FALSE)

test_that("gbm_grid parallel", {
  system.time(o <- gbm_grid(y=y, x=x, cv.folds=3, 
                            interaction.depth=1:3, shrinkage=c(.1,.5,.8),
                            n.trees=500, distribution="gaussian",
                            verbose=FALSE, mc.cores=6))
  expect_named(o)
  system.time(o <- gbm_grid(y=y, x=x, cv.folds=3, 
                            interaction.depth=1:3, shrinkage=c(.1, .5, .8),
                            n.trees=500, distribution="gaussian",
                            verbose=FALSE, mc.cores=1))  
  expect_named(o)
})

test_that("gbm_grid grid", {
  o <- gbm_grid(y=y, x=x, cv.folds=3, 
                interaction.depth=1:3,
                shrinkage=c(.5,1),
                n.minobsinnode=c(5,10), 
                bag.fraction=c(.5, 1),
                distribution="gaussian", 
                verbose=FALSE, 
                mc.cores=1)

  expect_equivalent(o$args[,-ncol(o$args)], args)
})

test_that("gbm_grid cv.folds=1", {
  expect_error(o <- gbm_grid(y=y, x=x[,1], cv.folds=1, mc.cores=1, distribution="gaussian"))
  o <- gbm_grid(y=y, x=x[,1,drop=F], cv.folds=1, mc.cores=1, subset=1:250, distribution="gaussian", verbose=FALSE)
  expect_named(o)
  o <- gbm_grid(y=y, x=x[,1,drop=F], cv.folds=1, mc.cores=3, subset=1:250, distribution="gaussian", verbose=FALSE)
  expect_named(o)
})

test_that("gbm_grid do_one_fold", {
  for(k in 1:3){
    s <- train[folds != k]
    set.seed(102)
    o <- gbm::gbm.fit(y=y[s], x=data.frame(x[s, ,drop=F]), 
                      distribution="gaussian", 
                      verbose=FALSE, n.trees=5)
    yhat <- predict(o, newdata=data.frame(x), n.trees=1:5)
    test_err <- rep(0, 5)
    for(i in 1:5){
      test_err[i] <- mean((y[-s] - yhat[-s,i])^2)
    }
    set.seed(102)
    one_fold <- mvtboost:::do_one_fold(k=k, folds=folds, train=1:n, y=y, x=x, 
                                       distribution="gaussian", verbose=FALSE, 
                                       n.trees=5)
    expect_equal(unname(test_err), unname(one_fold), info=paste0("k = ",k), 
                 tolerance=tol, check.attributes=F)
  }
})

test_that("gbm_grid do_one_row", {
  set.seed(101)
  one_fold <- mvtboost:::do_one_fold(y=y, x=x, train=train, k=1, folds=folds, 
                                     interaction.depth=1,
                                     shrinkage=.5,
                                     n.minobsinnode=5, 
                                     bag.fraction=.5,
                                     distribution="gaussian", verbose=FALSE)
                                      
  set.seed(101)
  one_row <- mvtboost:::do_one_row(1, args=args.ls, train=train, 
                                   folds=folds, Xm=x, y=y)
  expect_equal(one_fold, one_row)
  
})

test_that("gbm_grid aggregate_cv_err", {
  set.seed(102)
  folds <- sample(1:cv.folds, size=length(train), replace=TRUE)
  unique_ids <- rep(1:nrow(args), each=cv.folds)
  ocv <- lapply(seq_along(args.ls), mvtboost:::do_one_row, 
                args=args.ls, train=train, folds=folds, Xm=x, y=y)
  ag1 <- mvtboost:::aggregate_cv_err(ocv, unique_ids)
  ag2 <- list()
  for(i in 1:nrow(args)){
    cond_cv_err <- ocv[unique_ids==i]
    ag2[[i]] <- rowMeans(data.frame(cond_cv_err))
  }
  
  # tapply adds a couple of attributes, so we just add them here for 
  # the most stringent test
  
  attr(ag2, "dim") <- nrow(args)
  attr(ag2, "dimnames") <- list(paste(1:nrow(args)))

  expect_equal(ag1, ag2)
  
  cv_err <- sapply(ag2, function(e){min(e)})
  new_args <- args
  new_args$err <- cv_err
  bargs <- args[which.min(sapply(ag1, min)), ]
  
  set.seed(102)
  o <- gbm_grid(y=y, x=x, cv.folds=3, 
                interaction.depth=1:3,
                shrinkage=c(.5,1),
                n.minobsinnode=c(5,10), 
                bag.fraction=c(.5, 1),
                distribution="gaussian", 
                verbose=FALSE, 
                mc.cores=1)

  expect_equal(bargs, o$best_args)
  expect_equal(new_args, o$args)
  
})

test_that("gbm_grid replicability", {
  set.seed(104)
  o <- gbm_grid(y=y, x=x, cv.folds=3, mc.cores=1, subset=1:250, 
                distribution="gaussian", verbose=FALSE, bag.fraction=.5)
  set.seed(104)
  o2 <- gbm_grid(y=y, x=x, cv.folds=3, mc.cores=1, subset=1:250, 
                distribution="gaussian", verbose=FALSE, bag.fraction=.5)
  expect_equal(o, o2)
})