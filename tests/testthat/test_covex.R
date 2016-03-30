context("test_covex")
k <- 3
#set.seed(100)

# no error, linear
q <- 3
p <- 1
n <- 1000
B <- matrix(0,p,q)
B[1,c(1,3)] <- 1
X <- matrix(rnorm(n),n,1)

Y <- X %*% B 

out <- mvtb(Y=Y,X=X)
covex <- mvtb.covex(out, Y=Y,X=X)
expect_output(print(covex),"")
expect_equal(dim(covex),c((q*(q+1))/2,p))

a <- stats::cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)

test_that("covex-noE", { 
for(i in 1:length(shrink)){
  out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=1000)
  covex <- mvtb.covex(out,Y=Y,X=X)
  expect_equal(covex[c(2,4,5)],rep(0,3))
  expect_true(all(abs(covex[covex > 0] - a) < .01))
}
})


E <- matrix(rnorm(n*q,0,.0001),n,q)
Y <- X %*% B + E
a <- stats::cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)

test_that("covex-E", {
  for(i in seq_along(shrink)){
    out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=1000,cv.folds=3)
    covex <- mvtb.covex(out, Y=Y, X=X)
    expect_equal(covex[c(2,4,5)],rep(0,3), tolerance=.001, info=shrink[i])
    expect_equal(covex[covex > .01],rep(a,3), info=shrink[i], tolerance=.01)
  }
})

## No approximation error
Xb <- ifelse(X < .5,0,1)
Y <- Xb %*% B
a <- stats::cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)
test_that("covex-exact", {
for(i in seq_along(shrink)){
  out <- mvtb(Y=Y,X=Xb,shrinkage=shrink[i],n.trees=100)
  covex <- mvtb.covex(Y=Y,X=Xb,out)
  # expect_equal(out$covex[c(2,4,5)],rep(0,3))
  expect_equal(covex[covex > .01],rep(a,3), tolerance=1E-10)
}
})
## test alpha

for(i in c(.25,.5,.75)) {
  #    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=1,weight.type=1,s=1:nrow(X))
#    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=1)),r$wm.raw[1])
  #    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=2,weight.type=1,s=1:nrow(X))
  #    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=2)),r$wm.raw[1])
  #    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=3,weight.type=1,s=1:nrow(X))
  #    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=3)),r$wm.raw[1])    
  #    r <- mvtb(X=X,Y=Y[,1,drop=F],n.trees=1,alpha=i, shrinkage=1,trainfrac=1,bag.frac=1,cov.discrep=4,weight.type=1,s=1:nrow(X))
  #    expect_equal(c(eval.loss(r$resid,Y[,1,drop=F],alpha=i,type=4)),r$wm.raw[1])        
}


## stop crit - no longer necessary. all weights should be positive

#r <- mvtb(X=X,Y=Y,n.trees=103,alpha=.5, trainfrac=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5)
#expect_equal(r$maxiter,103)
Y <- E
for(i in 1:4) { 
  #r <- mvtb(X=X,Y=Y,n.trees=25,alpha=.5, shrinkage=1,trainfrac=1,samp.iter=FALSE,cov.discrep=i,weight.type=2,bag.frac=.5)
  #    expect_true(all(r$wm >= 0),info=paste0("not all weights > 0, loss ",i))
}
