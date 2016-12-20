context("test_covex")

k <- 3
q <- 3
p <- 1
n <- 1000
B <- matrix(0,p,q)
B[1,c(1,3)] <- 1
X <- matrix(rnorm(n),n,1)

Y <- X %*% B 

test_that("mvtb.covex", {
  
  out <- mvtb(Y=Y, X=X)
  covex <- mvtb.covex(out, Y=Y, X=X)
  
  expect_output(print(covex),"")
  expect_equal(dim(covex), c((q*(q+1))/2, p))
  
  expect_is(mvtb.covex(out,Y=Y,X=X),"matrix")
  expect_is(mvtb.covex(o, Y=y, X=x),"matrix")
  expect_is(mvtb.covex(o, Y=y, X=x, iter.details=T),"list")
})


test_that("covex-exact", {
  Xb <- ifelse(X < .5,0,1)
  Y <- Xb %*% B
  a <- stats::cov(Y)[1,3]
  
  out <- mvtb(Y=Y, X=Xb, shrinkage=1, n.trees=100)
  covex <- mvtb.covex(Y=Y,X=Xb,out)
  # expect_equal(out$covex[c(2,4,5)],rep(0,3))
  expect_equal(covex[covex > .01],rep(a,3), tolerance=1E-10)
})

