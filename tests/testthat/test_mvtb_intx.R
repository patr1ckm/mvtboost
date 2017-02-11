context("test_mvtb_intx")
set.seed(123)
n <- 1000
B <- matrix(0, nrow=5, ncol=4)
B[3,1:2] <- 0
B[2,2:3] <- 0
B[4,1] <- 1
B[5,3:4] <- 1

X <- matrix(rbinom(n*(nrow(B)-2), size=1, prob=.5), n, nrow(B)-2)
X2 <- cbind(x1x2=X[,1]*X[,2], x2x3=X[,2]*X[,3])
X <- cbind(X, X2)
E <- matrix(rnorm(n*4),nrow=n,ncol=4)
Y <- X %*% B + E

out <- mvtb(Y=Y, X=X[,1:3], n.trees=50,interaction.depth = 5, shrinkage = .5)

test_that("intx runs", {
  expect_is(mvtb.nonlin(out, X=X[,1:3], Y=Y, n.trees=50, detect = "grid"), "list")
  expect_is(mvtb.nonlin(out, X=X[,1:3], Y=Y, n.trees=50, detect = "influence"), "list")
  expect_is(mvtb.nonlin(out, X=X[,1:3], Y=Y, n.trees=50, detect = "lm"), "list")
})
