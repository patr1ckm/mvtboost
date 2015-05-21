#library(testthat)
#setwd( "/Users/pmille13/Documents/Projects/mvtboost/Tests/")
#source("../13/mvtboost_v14.R")
#source("wt_helper_functions.R")
library(MASS)

context("test_params")
fp <- paste0(getwd(),"/test_params.R")

## purpose: make sure all the buttons work
set.seed(123)
n <- 1000
B <- matrix(0,nrow=3,ncol=4)
B[3,1:2] <- 2
B[2,2:3] <- 1
B[1,1] <- .5
X <- matrix(rbinom(n*nrow(B),size=1,prob=.5),n,nrow(B))
E <- mvrnorm(n,rep(0,4),Sigma=diag(4))
Y <- X %*% B + E

## test n.trees
context("n.trees")
r <- mvtb(X=X,Y=Y,n.trees=50)
expect_equal(r$maxiter,50)

## test alpha
context("alpha")
for(i in c(.25,.5,.75)) {
    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=1,weight.type=1,s=1:nrow(X))
    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=1)),r$wm.raw[1])
    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=2,weight.type=1,s=1:nrow(X))
    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=2)),r$wm.raw[1])
    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=3,weight.type=1,s=1:nrow(X))
    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=3)),r$wm.raw[1])    
    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=4,weight.type=1,s=1:nrow(X))
    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=4)),r$wm.raw[1])        
}

## test trainfrac
context("trainfrac")
for(i in seq(.1,.9,by=.1)) {
    r <- mvtb(X=X,Y=Y,n.trees=1,alpha=i, trainfrac=i,cov.discrep=1,weight.type=2)
    expect_equal(r$models[[1]]$nTrain,floor(n*i))
}

## test samp.iter
context("samp.iter")
r <- mvtb(X=X,Y=Y,n.trees=10,alpha=.5, trainfrac=1,samp.iter=TRUE,cov.discrep=1,weight.type=2,save.cv=TRUE)
#expect_equal(apply(r$s,2,function(col) {length(unique(col))}),rep(n,ncol(r$s))) ## all unique

ord.s <- apply(r$s,2,function(col){col[order(col)]})
expect_true(!any(apply(combn(5,2),2,function(idx){all(ord.s[,idx[1]]==ord.s[,idx[2]])})),info="one or more iterations have identical observations") 

r <- mvtb(X=X,Y=Y,n.trees=10,alpha=.5, trainfrac=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,save.cv=TRUE)
a_ply(r$s,1,function(row){expect_equal(row,rep(row[sample(1:length(row),1)],length(row)))}) # grab a random element, check to see if the whole row is equal to it


## test bag.frac
context("bag.frac")
r <- mvtb(X=X,Y=Y,n.trees=10,alpha=.5, trainfrac=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5)
for(i in 1:4) { expect_equal(r$models[[i]]$bag.fraction,.5) }

## stop crit - no longer necessary. all weights should be positive
context("all wm positive")
r <- mvtb(X=X,Y=Y,n.trees=103,alpha=.5, trainfrac=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5)
expect_equal(r$maxiter,103)
Y <- E
for(i in 1:4) { 
    r <- mvtb(X=X,Y=Y,n.trees=25,alpha=.5, shrinkage=1,trainfrac=1,samp.iter=FALSE,cov.discrep=i,weight.type=2,bag.frac=.5)
    expect_true(all(r$wm >= 0),info=paste0("not all weights > 0, loss ",i))
}

## s
context("s")
r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,save.cv=TRUE)
a_ply(r$s,1,function(row){expect_equal(row,rep(row[sample(1:length(row),1)],length(row)))}) # grab a random element, check to see if the whole row is equal to it
#expect_error(r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:10))
expect_true(all(r$s %in% 1:500))

## seednum
context("seednum")
r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8)
expect_equal(r$params$seednum,8)
#expect_equal(r$params$seed,55543*8)
r1 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=TRUE,cov.discrep=1,weight.type=2,bag.frac=.5,seednum=8)
r2 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=TRUE,cov.discrep=1,weight.type=2,bag.frac=.5,seednum=8)
expect_equal(r1,r2)
r3 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=TRUE,cov.discrep=1,weight.type=2,bag.frac=.5,seednum=9)
expect_true(!identical(r1,r3))

context("compress")
r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE,cv.folds=3,save.cv=TRUE)
r1 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, trainfrac=.5,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=TRUE,cv.folds=3,save.cv=TRUE)

uncomp.l <- function(l) { lapply(l,function(li){unserialize(memDecompress(li,type="bzip2"))}) }
r2 <- uncomp.l(r1)
r$params$compress <- TRUE # set to TRUE so that the comparison is legitimate
expect_equal(r,r2)

r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE)
r1 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=TRUE)
r$params$compress <- TRUE
r2 <- uncomp.l(r1)
expect_equal(r,r2)

context("parallel")

r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, cv.folds=3,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE,mc.cores=3)



context("trainfrac, samp.iter, and s")
## this tests s
for(i in seq(.1,.9,by=.1)) {    
    r <- mvtb(X=X,Y=Y,n.trees=10,alpha=.5, trainfrac=i,samp.iter=TRUE,cov.discrep=1,weight.type=2,s=1:floor(n*i),save.cv=TRUE)
    #expect_equal(apply(r$s,2,function(col) {length(unique(col))}),rep(floor(n*i),ncol(r$s))) ## 500 unique
    expect_true(all(apply(r$s,2,function(col){all(col %in% 1:floor(n*i))})))
    #for(j in 1:10) {
    #    expect_equal(sum(1:n %in% r$s[,j]),floor(n*i))
    #}
}


context("interaction depth")
r <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE,interaction.depth=1)
r2 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE,interaction.depth=2)
r3 <- mvtb(X=X,Y=Y,n.trees=5,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE,interaction.depth=3)
expect_equal(r2$params$interaction.depth,2)

# totalnumber of nodes for a given "interaction depth" (which is the number of splits) is 3*n + 1:
# = {l,r,NA}*n + root
# interaction depth = number of splits.
n <- 1:3
nodes <- 3*n+1
expect_true(all(unlist(lapply(r$finaltree[[1]],function(t){length(t[[1]])})) == nodes[1]))
expect_true(all(unlist(lapply(r2$finaltree[[1]],function(t){length(t[[1]])})) == nodes[2]))
expect_true(all(unlist(lapply(r3$finaltree[[1]],function(t){length(t[[1]])})) == nodes[3]))


#r <- mvtb(X=X,Y=Y,n.trees=100,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE)
#r1 <- mvtb(X=X,Y=Y,n.trees=100,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=FALSE,save.cv=TRUE)
#r <- mvtb(X=X,Y=Y,n.trees=100,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=TRUE)
#r1 <- mvtb(X=X,Y=Y,n.trees=100,alpha=.5, cv.folds=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5,s=1:500,seednum=8,compress=TRUE,save.cv=TRUE)
