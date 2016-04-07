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

check.samp <- function(ocv,s=1:500,folds=cv.folds,n=1000) {
  ## get the training sample used in every fold
  fold.obs <- lapply(ocv$models.k[1:folds],function(out){unique(out$s)})
  ## check which obs in s (1:500) are in each f. The row sums of the resulting matrix should be k-1
  expect_true(all(rowSums(t(plyr::laply(fold.obs,function(f){matrix(s %in% f)})))==(folds-1)))
  expect_equal(sum(sapply(fold.obs,length)),n)
  for(k in 1:cv.folds) {
    expect_true(!any(intersect(which(ocv$cv.groups==k),fold.obs[[k]])))
  }
}

test_that("mvtbCV", {
  ocv <- mvtboost:::mvtbCV(Y=Y, X=X, distribution="gaussian", n.trees=n.trees, cv.folds=3, s=s, save.cv=TRUE, mc.cores=1, verbose=F)
  # 1. check that each observation is left out once.
  check.samp(ocv)
  
  # 2. check that observations in k are not in training set

  ocv <- mvtboost:::mvtbCV(Y=Y, X=X, distribution="gaussian", n.trees=n.trees, cv.folds=3, s=s, save.cv=TRUE, mc.cores=1, verbose=F)
  k.obs <- lapply(1:3,function(k){which(ocv$cv.groups==k)})
  fold.obs <- lapply(ocv$models.k[1:3],function(out){unique(out$s)})
  expect_true(length(unlist(mapply(intersect,k.obs,fold.obs))) < 1)
  
  # 3. no intersection between the observations supposed to be in the kth fold, and in the training set
  # observations in a training fold are only from s at each iteration
  expect_true(all(unlist(lapply(fold.obs,function(f){
      all(unlist(lapply(f,function(col){col %in% s})))
  }))))
})

test_that("mvtb - CV param", {
# 4. check an initial split of training and test, with CV only in the training set
out <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,save.cv=TRUE)
out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,train.fraction=.5,shrinkage=.5,cv.folds=3,s=NULL,save.cv=TRUE)
check.samp(out$ocv,s=out$params$s,folds=3)
s <- out2$s
fold.obs <- lapply(out2$ocv$models.k[1:3],function(out){unique(out$s)})
expect_true(all(unique(unlist(fold.obs)) %in% s))
expect_true(all(s %in% unique(unlist(fold.obs))))

expect_true(all(!unlist(lapply(out$ocv$models.k,function(m){any(m$s > 500)})))) # expect that not any observation numbers > 500 (s was 500)

# 5. check that setting the seed obtains the same observations in each fold 
out1 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
out2 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
out2 <- mvtb(X=X,Y=Y,train.fraction=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)

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




