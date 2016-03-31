context("test_params")
set.seed(123)
n <- 1000
B <- matrix(0,nrow=3,ncol=4)
B[3,1:2] <- 2
B[2,2:3] <- 1
B[1,1] <- .5
X <- matrix(rbinom(n*nrow(B),size=1,prob=.5),n,nrow(B))
E <- MASS::mvrnorm(n,rep(0,4),Sigma=diag(4))
Y <- X %*% B + E

test_that("n.trees", {
r <- mvtb(X=X,Y=Y,n.trees=50)
expect_equal(r$best.trees$last,50)
})

test_that("train.fraction", {
for(i in seq(.1,.9,by=.1)) {
    r <- mvtb(X=X,Y=Y,n.trees=1,train.fraction=i)
    expect_equal(r$models[[1]]$nTrain,floor(n*i))
}
})

test_that("bag.fraction", {
r <- mvtb(X=X,Y=Y,n.trees=10,alpha=.5, trainfrac=1,samp.iter=FALSE,cov.discrep=1,weight.type=2,bag.frac=.5)
for(i in 1:4) { expect_equal(r$models[[i]]$bag.fraction,.5) }
})

test_that("subsetting", {
r <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5,bag.fraction=.5,s=1:500,save.cv=TRUE)
expect_equal(r$s, 1:500) # should be in this order correct order
})

test_that("seednum", {
r <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8)
expect_equal(r$params$seednum,8)
r2 <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8)
expect_equal(r,r2)
r3 <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=9)
expect_true(!identical(r,r3))
})

test_that("compress", {
  r <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=FALSE)
  expect_true(all(!(unlist(lapply(r,class)) %in% "raw"))) # none of the objects are raw
  r <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=TRUE)
  expect_equal_to_reference(unlist(lapply(r,class)),"raw") # all of the objects are raw
  r <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=FALSE, cv.folds=3, save.cv=T)
  expect_true(all(!(unlist(lapply(r,class)) %in% "raw"))) # none of the objects are raw
  expect_true(all(!(unlist(lapply(r$ocv,class)) %in% "raw"))) # none of the cv objects are raw
  r1 <- mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=TRUE, cv.folds=3, save.cv=T)
  expect_equal_to_reference(unlist(lapply(r1,class)),"raw") # all of the objects are raw
  expect_equal(class(r1$ocv),"raw") # the ocv object is raw
  expect_lt(object.size(r1), object.size(r)) # just verify that the object sizes are as they should be
})

test_that("mvtb.uncomp", {
  rc <- mvtb(X=X,Y=Y,n.trees=5, seednum=8, compress=TRUE)
  r  <- mvtb(X=X,Y=Y,n.trees=5, seednum=8, compress=FALSE)
  r2 <- mvtb.uncomp(rc)
  r$params$compress <- TRUE # set to TRUE so that the comparison is legitimate
  expect_equal(r,r2)
})

test_that("iter.details", {
r <- mvtb(X=X,Y=Y,n.trees=5, seednum=8, compress=FALSE, cv.folds=3, save.cv=T, iter.details = T)
expect_true(all(c("trainerr", "testerr", "cv.err", "ocv") %in% names(r)))
expect_length(r$testerr, r$best.trees$last)
expect_length(r$trainerr, r$best.trees$last)
expect_length(r$cv.err, r$best.trees$last)

r <- mvtb(X=X,Y=Y,n.trees=5, seednum=8, compress=FALSE, cv.folds=3, save.cv=F, iter.details = F)
expect_named(r,c("models", "best.trees", "params", "s", "ocv", "n", "xnames", "ynames"))
expect_null(r$ocv)
r <- mvtb(X=X,Y=Y,n.trees=5, seednum=8, compress=FALSE, cv.folds=3, save.cv=T, iter.details = F)
})

test_that("verbose", {
expect_silent(mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=FALSE, cv.folds=3, save.cv=T, verbose = F))
expect_output(mvtb(X=X,Y=Y,n.trees=5, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=FALSE, cv.folds=3, save.cv=T, verbose = T))
})

#context("parallel")

# r <- mvtb(X=X,Y=Y,n.trees=500, train.fraction=.5, bag.fraction=.5,s=1:500,seednum=8, compress=FALSE, cv.folds=3, save.cv=T, mc.cores=3)

test_that("keep.data", {
  r <- mvtb(X=X,Y=Y,n.trees=50, keep.data=FALSE)
  expect_null(r$models[[1]]$data)
  r <- mvtb(X=X,Y=Y,n.trees=50, keep.data=TRUE)
  expect_is(r$models[[1]]$data,"list")
})

test_that("distribution", {
r <- mvtb(X=X,Y=Y,n.trees=50)
expect_equal(r$models[[1]]$distribution$name,"gaussian") # default is gaussian
r <- mvtb(X=X,Y=Y,n.trees=50, distribution="gaussian")   
expect_equal(r$models[[1]]$distribution$name,"gaussian") 
expect_error(mvtb(X=X,Y=Y,n.trees=50, distribution="bernoulli")) # correctly passes bernoulli to gbm 
})

test_that("train.fraction and s", {
## Test 
for(i in seq(.1,.9,by=.1)) {    
    r <- mvtb(X=X,Y=Y,n.trees=10,train.fraction=i,s=1:floor(n*i),save.cv=TRUE)
    expect_true(all(r$s %in% 1:floor(n*i)))
}
})

test_that("interaction depth", {
r <- mvtb(X=X,Y=Y,n.trees=50, interaction.depth=1)
r2 <- mvtb(X=X,Y=Y,n.trees=50, interaction.depth=2)
r3 <- mvtb(X=X,Y=Y,n.trees=50, interaction.depth=2)

expect_equal(r2$params$interaction.depth,2)

# totalnumber of nodes for a given "interaction depth" (which is the number of splits) is 3*n + 1:
# = {l,r,NA}*n + root
# interaction depth = number of splits.
n <- 1:3
nodes <- 3*n+1
expect_true(all(unlist(lapply(r$finaltree[[1]],function(t){length(t[[1]])})) == nodes[1]))
expect_true(all(unlist(lapply(r2$finaltree[[1]],function(t){length(t[[1]])})) == nodes[2]))
expect_true(all(unlist(lapply(r3$finaltree[[1]],function(t){length(t[[1]])})) == nodes[3]))
})

test_that("checks", {
# test errors
expect_error(mvtb(X=X,Y=Y,n.trees=50, shrinkage=2))
expect_error(mvtb(X=X,Y=Y,n.trees=50, shrinkage=0))
expect_error(mvtb(X=X,Y=Y,n.trees=50, shrinkage=-1))
expect_error(mvtb(X=X,Y=Y,n.trees=50, train.fraction=2))
expect_error(mvtb(X=X,Y=Y,n.trees=50, train.fraction=0))
expect_error(mvtb(X=X,Y=Y,n.trees=50, train.fraction=-1))
expect_error(mvtb(X=X,Y=Y,n.trees=50, bag.fraction=-1))
expect_error(mvtb(X=X,Y=Y,n.trees=50, bag.fraction=0))
expect_error(mvtb(X=X,Y=Y,n.trees=50, bag.fraction=2))
Y[1,1] <- NA
expect_error(mvtb(X=X,Y=Y,n.trees=50, bag.fraction=2))
})

test_that("input", {
# test data frame
Y <- X %*% B + E
Xf <- as.data.frame(X)
Yf <- as.data.frame(Y)
out <- try(mvtb(Y=Yf,X=Xf))
expect_is(out,"mvtb")

# test single predictor case
set.seed(123)
n <- 1000
B <- matrix(0,nrow=1,ncol=4)
B[1,1:2] <- 1
X <- matrix(rbinom(n,size=1,prob=.5),n,nrow(B))
E <- matrix(rnorm(n*4),nrow=n,ncol=4)
Y <- X %*% B + E
expect_is(mvtb(Y=Y,X=X),"mvtb")
expect_is(mvtb(Y=Y,X=as.data.frame(X)), "mvtb")
Xf <- as.factor(X)
expect_is(mvtb(Y=Y,X=as.data.frame(X)), "mvtb")

# test single outcome, single predictor
set.seed(123)
n <- 1000
B <- matrix(0,nrow=1,ncol=1)
B[1,1] <- 1
X <- matrix(rbinom(n,size=1,prob=.5),n,nrow(B))
E <- matrix(rnorm(n*nrow(B)),nrow=n,ncol=nrow(B))
Y <- X %*% B + E
expect_is(mvtb(Y=Y,X=X), "mvtb")

# test single outcome, single predictor
set.seed(123)
n <- 1000
B <- matrix(0,nrow=1,ncol=1)
B[1,1] <- 1
X <- matrix(rbinom(n,size=1,prob=.5),n,nrow(B))
E <- matrix(rnorm(n*nrow(B)),nrow=n,ncol=nrow(B))
Y <- X %*% B + E
yf <- as.numeric(Y < 0)
expect_is(mvtb(Y=yf,X=X,distribution="bernoulli"), "mvtb")

# test vectors
expect_is(mvtb(Y=Y[,,drop=TRUE],X=X[,,drop=TRUE]), "mvtb")

x <- rnorm(1000)
y <- x*5 + rnorm(1000)
expect_is(mvtb(Y=y,X=x), "mvtb")
})


