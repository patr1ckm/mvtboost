context("test_cv")
#fp <- paste0(getwd(),"/test_cv.R")

set.seed(123)
n <- 1000
B <- matrix(0,nrow=3,ncol=4)
B[3,1:2] <- 2
B[2,2:3] <- 1
B[1,1] <- .5
X <- matrix(rbinom(n*nrow(B),size=1,prob=.5),n,nrow(B))
E <- matrix(rnorm(n*ncol(B)),nrow=n,ncol=ncol(B))
Y <- X %*% B + E

n.trees <- 25

## check cv, s, seednum

s <- 1:500
cv.folds <- 3
save.cv=TRUE

test_that("mvtb - CV param", {

  out1 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
  out2 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
  
  expect_equal(out1,out2)
  
  out1 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
  out2 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
  
  expect_equal(out1,out2)
  
  out1 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,bag.frac=.5,seednum=1)
  out2 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,bag.frac=.5,seednum=1)
  
  expect_equal(out1,out2)
})

## 6. Final model in cv.folds=3 should be the same as cv.folds =1
## This should be true for all combinations of train.fraction, bag.fraction, and s
test_that("final_model", {
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,s=1:1000,seednum=1)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,s=1:1000,seednum=1)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2) 
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,train.fraction=.5,cv.folds=3,compress=F,seednum=1)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,train.fraction=.5,cv.folds=1,compress=F,seednum=1)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1,s=1:500)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1,s=1:500)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1,bag.fraction=.5)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1,bag.fraction=.5)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1,bag.frac=.5,s=1:500)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1,bag.frac=.5,s=1:500)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1,bag.frac=.5,train.fraction=.5)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1,bag.frac=.5,train.fraction=.5)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
})




