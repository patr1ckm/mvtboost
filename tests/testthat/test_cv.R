
#library(testthat)
#setwd( "/Users/pmille13/Documents/Projects/mvtboost/Tests/")
#source("../13/mvtboost_v14.R")
#source("wt_helper_functions.R")
#library(MASS)
library(plyr)

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

ocv <- mvtbCV(params=plist)

context("check samples")
check.samp <- function(ocv,s=1:500,folds=params$cv.folds,n=1000) {
## get the training sample used in every fold
    fold.obs <- lapply(ocv$models.k[1:folds],function(out){unique(out$ss[,1])})
## check which obs in s (1:500) are in each f. The row sums of the resulting matrix should be k-1
    expect_true(all(rowSums(t(laply(fold.obs,function(f){matrix(s %in% f)})))==(folds-1)))
    expect_equal(sum(sapply(fold.obs,length)),n)
    for(k in 1:params$cv.folds) {
        expect_true(!any(intersect(which(ocv$cv.groups==k),fold.obs[[k]])))
    }
}
check.samp(ocv)

# 1. check dimensions of s (should be n*(k-1)/k x niter)

expect_true(all(data.frame(lapply(ocv$models.k,function(out){dim(out$ss)}))[2,]==100))

# check that observations in k are not in s[-k] (are not trained on, no intersect)

# 2. check that samp.iter produces new draws of only -k samples, equivalently: observations in k are not in training set
params$samp.iter <- TRUE
ocv <- mvtbCV(params=plist)
k.obs <- lapply(1:3,function(k){which(ocv$cv.groups==k)})
fold.obs <- lapply(ocv$models.k[1:3],function(out){apply(out$ss,2,unique)})
expect_true(length(unlist(mapply(intersect,k.obs,fold.obs))) < 1)
# no intersection between the observations supposed to be in the kth fold, and in the training set
# observations in a training fold are only from s at each iteration
expect_true(all(unlist(lapply(fold.obs,function(f){
    all(unlist(lapply(f,function(col){col %in% s})))
}))))

context("mvtb - CV param")
niter <- 25

# 4. check an initial split of training and test, with CV only in the training set
out <- mvtb(X=X,Y=Y,s=1:500,niter=niter,shrinkage=.5,cv.folds=3,save.cv=TRUE)
out2 <- mvtb(X=X,Y=Y,niter=niter,trainfrac=.5,shrinkage=.5,cv.folds=3,samp.iter=FALSE,s=NULL,save.cv=TRUE)
check.samp(out$ocv,s=out$params$s,folds=3)
check.samp(out2$ocv,s=out$params$s,folds=3)


expect_true(all(!unlist(lapply(out$ocv$models.k,function(m){any(m$s > 500)})))) # expect that not any observation numbers > 500 (s was 500)
expect_true(all(!unlist(lapply(out2$ocv$models.k,function(m){any(m$s > 500)})))) # expect that not any observation numbers > 500 (n*trainfrac = 500)

# this is only a minor check of the train fraction. In CV, we only use X[s,],Y[s,], which are then subsequently indexed by 1:length(s) not (1:n)[s]
## To really check this you need to check that the rownames of X (within each cv group) are always in s, not just row indices s.cv
check.samp(out2$ocv,s=1:500,folds=3) 

# 5. check that setting the seed obtains the same observations in each fold (with & w/o samp.iter)

out1 <- mvtb(X=X,Y=Y,s=1:500,niter=niter,shrinkage=.5,cv.folds=3,seednum=1)
out2 <- mvtb(X=X,Y=Y,s=1:500,niter=niter,shrinkage=.5,cv.folds=3,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,trainfrac=.5,niter=niter,shrinkage=.5,cv.folds=3,seednum=1)
out2 <- mvtb(X=X,Y=Y,trainfrac=.5,niter=niter,shrinkage=.5,cv.folds=3,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,trainfrac=.5,niter=niter,shrinkage=.5,cv.folds=3,samp.iter=TRUE,seednum=1)
out2 <- mvtb(X=X,Y=Y,trainfrac=.5,niter=niter,shrinkage=.5,cv.folds=3,samp.iter=TRUE,seednum=1)

expect_equal(out1,out2)

out1 <- mvtb(X=X,Y=Y,trainfrac=.5,niter=niter,shrinkage=.5,cv.folds=3,samp.iter=TRUE,bag.frac=.5,seednum=1)
out2 <- mvtb(X=X,Y=Y,trainfrac=.5,niter=niter,shrinkage=.5,cv.folds=3,samp.iter=TRUE,bag.frac=.5,seednum=1)

expect_equal(out1,out2)

## okay, these are just some rough tests, but it seems like it works

## 6. test with comp=TRUE

out <- mvtb(X=X,Y=Y,s=s,niter=niter,shrinkage=.5,cv.folds=3,compress=T)


## 7. test cv error vs training error (trainfrac = 1). Expect that training errors at every iteration are < cv err
## context("cv error vs training error")
## niter <- 1000
## out <- mvtb(X=X,Y=Y,niter=niter,trainfrac=1,shrinkage=.5,cv.folds=3,compress=FALSE,stop.crit=TRUE)

## ## par(mfcol=c(2,2))
## ## plot(1:out$maxiter,out$testerr,type="l",sub="mse")
## ## points(1:out$maxiter,rowMeans(out$ocv$cv.err$testerr,na.rm=TRUE)[1:out$maxiter],type="l",col="red")
## ## plot(1:out$maxiter,out$srcm,type="l",sub="srcm")
## ## points(1:out$maxiter,rowMeans(out$ocv$cv.err$srcm,na.rm=TRUE)[1:out$maxiter],type="l",col="red")
## ## plot(1:out$maxiter,out$wsrcm,type="l",sub="wsrcm")
## ## points(1:out$maxiter,rowMeans(out$ocv$cv.err$wsrcm,na.rm=TRUE)[1:out$maxiter],type="l",col="red")

## expect_true(all(out$testerr < rowMeans(out$ocv$cv.err$testerr,na.rm=TRUE)[1:out$maxiter]))
## expect_true(all(out$srcm < rowMeans(out$ocv$cv.err$srcm,na.rm=TRUE)[1:out$maxiter]))
## expect_true(all(out$wsrcm < rowMeans(out$ocv$cv.err$wsrcm,na.rm=TRUE)[1:out$maxiter]))

## #dev.off()


## ## 8. test save.cv=FALSE
## context("save.cv")
## niter <- 25
## out <- mvtb(X=X,Y=Y,niter=niter,trainfrac=1,shrinkage=.5,cv.folds=3,compress=TRUE,stop.crit=TRUE,save.cv=TRUE)
## out1 <- mvtb(X=X,Y=Y,niter=niter,trainfrac=1,shrinkage=.5,cv.folds=3,compress=TRUE,stop.crit=TRUE,save.cv=FALSE)
## expect_true(object.size(out) > object.size(out1))
## expect_true(is.null(uncomp(out1$ocv)$models.k))

