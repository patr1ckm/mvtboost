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

a <- stats::cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)

for(i in 1:length(shrink)){
  out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=1000)
  covex <- mvtb.covex(Y=Y,X=X,object=out)
  expect_equal(covex[c(2,4,5)],rep(0,3))
  expect_true(all(abs(covex[covex > 0] - a) < .01))
}

E <- matrix(rnorm(n*q,0,.0001),n,q)
Y <- X %*% B + E
a <- stats::cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)
for(i in seq_along(shrink)){
  out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=5000)
 expect_equal(out$covex[c(2,4,5)],rep(0,3))
  expect_true(all(abs(out$covex[out$covex > .01] - a) < .01))
}

## No approximation error
Xb <- ifelse(X < .5,0,1)
Y <- Xb %*% B
a <- stats::cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)
for(i in seq_along(shrink)){
  out <- mvtb(Y=Y,X=Xb,shrinkage=shrink[i],n.trees=100)
  covex <- mvtb.covex(Y=Y,X=Xb,out)
  # expect_equal(out$covex[c(2,4,5)],rep(0,3))
  expect_true(all(abs(covex[covex > .01] - a) < 1E-10))
}