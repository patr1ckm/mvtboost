context("test_mvtb_cv")
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

  set.seed(1)
  out1 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3)
  set.seed(1)
  out2 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3)
  
  expect_equal(out1,out2)
  
  set.seed(1)
  out1 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3)
  set.seed(1)
  out2 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3)
  
  expect_equal(out1,out2)
  
  set.seed(1)
  out1 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,bag.frac=.5)
  set.seed(1)
  out2 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,bag.frac=.5)
  
  expect_equal(out1,out2)
})

## This should be true for all combinations of train.fraction, and s

test_that("final_model", {
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,s=1:1000)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,s=1:1000)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2) 
  
  set.seed(104)
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F)
  set.seed(104)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  set.seed(104)
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,train.fraction=.5,cv.folds=3,compress=F)
  set.seed(104)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,train.fraction=.5,cv.folds=1,compress=F)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
  out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,s=1:500)
  out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,s=1:500)
  out$params <- out2$params
  out$best.trees <- out2$best.trees
  out$cv.err <- out2$cv.err <- NULL
  expect_equal(out,out2)
  
})




