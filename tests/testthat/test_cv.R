
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

context("mvtbCV")
## check cv, s, samp.iter, seednum

params <- formals(mvtb)[-c(length(formals(mvtb)))]
params$s <- s <- 1:500
params$cv.folds <- 3
params$X <- X
params$Y <- Y
params$save.cv=TRUE
plist <- params
plist$alpha <- NULL
plist$cov.discrep <- NULL
plist$weight.type <- NULL


# 0. check that each observation is left out once. (samp.iter off)

ocv <- mvtboost:::mvtbCV(params=plist)

context("check samples")
check.samp <- function(ocv,s=1:500,folds=params$cv.folds,n=1000) {
## get the training sample used in every fold
    fold.obs <- lapply(ocv$models.k[1:folds],function(out){unique(out$s)})
## check which obs in s (1:500) are in each f. The row sums of the resulting matrix should be k-1
    expect_true(all(rowSums(t(plyr::laply(fold.obs,function(f){matrix(s %in% f)})))==(folds-1)))
    expect_equal(sum(sapply(fold.obs,length)),n)
    for(k in 1:params$cv.folds) {
        expect_true(!any(intersect(which(ocv$cv.groups==k),fold.obs[[k]])))
    }
}
check.samp(ocv)

# check that observations in k are not in s[-k] (are not trained on, no intersect)

# 2. check that samp.iter produces new draws of only -k samples, equivalently: observations in k are not in training set
params$samp.iter <- TRUE
ocv <- mvtboost:::mvtbCV(params=plist)
k.obs <- lapply(1:3,function(k){which(ocv$cv.groups==k)})
fold.obs <- lapply(ocv$models.k[1:3],function(out){unique(out$s)})
expect_true(length(unlist(mapply(intersect,k.obs,fold.obs))) < 1)
# no intersection between the observations supposed to be in the kth fold, and in the training set
# observations in a training fold are only from s at each iteration
expect_true(all(unlist(lapply(fold.obs,function(f){
    all(unlist(lapply(f,function(col){col %in% s})))
}))))

context("mvtb - CV param")
n.trees <- 25

# 4. check an initial split of training and test, with CV only in the training set
out <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,save.cv=TRUE)
out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,trainfrac=.5,shrinkage=.5,cv.folds=3,samp.iter=FALSE,s=NULL,save.cv=TRUE)
check.samp(out$ocv,s=out$params$s,folds=3)
s <- out2$s
fold.obs <- lapply(out2$ocv$models.k[1:3],function(out){unique(out$s)})
expect_true(all(unique(unlist(fold.obs)) %in% s))
expect_true(all(s %in% unique(unlist(fold.obs))))

expect_true(all(!unlist(lapply(out$ocv$models.k,function(m){any(m$s > 500)})))) # expect that not any observation numbers > 500 (s was 500)

# 5. check that setting the seed obtains the same observations in each fold (with & w/o samp.iter)

out1 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
out2 <- mvtb(X=X,Y=Y,s=1:500,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,trainfrac=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)
out2 <- mvtb(X=X,Y=Y,trainfrac=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,trainfrac=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,samp.iter=TRUE,seednum=1)
out2 <- mvtb(X=X,Y=Y,trainfrac=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,samp.iter=TRUE,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,trainfrac=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,samp.iter=TRUE,bag.frac=.5,seednum=1)
out2 <- mvtb(X=X,Y=Y,trainfrac=.5,n.trees=n.trees,shrinkage=.5,cv.folds=3,samp.iter=TRUE,bag.frac=.5,seednum=1)

expect_equal(out1,out2)

## okay, these are just some rough tests, but it seems like it works

## 6. test with comp=TRUE

out <- mvtb(X=X,Y=Y,s=s,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=T)

## 7. Compare cv.folds=3 to cv.folds =1
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

out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,trainfrac=.5,cv.folds=3,compress=F,seednum=1)
out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,trainfrac=.5,cv.folds=1,compress=F,seednum=1)
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

out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1,bag.frac=.5)
out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1,bag.frac=.5)
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

out <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=3,compress=F,seednum=1,bag.frac=.5,trainfrac=.5)
out2 <- mvtb(X=X,Y=Y,n.trees=n.trees,shrinkage=.5,cv.folds=1,compress=F,seednum=1,bag.frac=.5,trainfrac=.5)
out$params <- out2$params
out$best.trees <- out2$best.trees
out$cv.err <- out2$cv.err <- NULL
expect_equal(out,out2)



