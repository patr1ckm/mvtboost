context("test_summary")

set.seed(123)
n <- 1000
B <- matrix(0,nrow=5,ncol=4)
B[3,1:2] <- 0
B[2,2:3] <- 0
B[4,1] <- 1
B[5,3:4] <- 1

X <- matrix(rbinom(n*(nrow(B)-2),size=1,prob=.5),n,nrow(B)-2)
X2 <- cbind(x1x2=X[,1]*X[,2],x2x3=X[,2]*X[,3])
Xf <- cbind(X,X2)
E <- matrix(rnorm(n*4),nrow=n,ncol=4)
Y <- Xf %*% B + E
p <- 3
ncovs <- 10
q <- 4
out <- mvtb(Y=Y,X=X,shrinkage=.1,n.trees=100)

test_that("summary",{ 
  ## Again, tests just to make sure that they run
  expect_output(print(summary(out)),"trees")
  expect_output(summary(out),"influence")
  expect_equal(sum(summary(out,print=FALSE,relative="tot")$relative.influence),100)
})

test_that("mvtb.cluster",{
  covex <- mvtb.covex(out, Y=Y,X=X)
  expect_output(print(mvtb.cluster(covex)))
  expect_equal(dim(mvtb.cluster(covex)),c(ncovs,p))
  mvtb.cluster(covex,plot=TRUE)
  expect_output(print(mvtb.cluster(covex,dist.method="manhattan",clust.method="complete")))
})

test_that("mvtb.ri",{
  expect_output(print(mvtb.ri(out)))
  expect_equal(dim(mvtb.ri(out)),c(p,q))
  expect_output(print(mvtb.ri(out)))
  expect_equivalent(sum(mvtb.ri(out,relative = "tot")),100)
  expect_equal(sum(colSums(mvtb.ri(out,relative = "col")))-q*100,0,tolerance=1E-12)
  expect_output(print(mvtb.ri(out,relative = "n")))
  expect_output(print(out), "List of ") # verifies that print.mvtb is being called
})
  
test_that("test_summary_compression",{
  out <- mvtb(Y=Y,X=X,compress = T)
  expect_output(summary(out))
  expect_output(print(out), "List of")
  covex <- mvtb.covex(out, Y=Y, X=X)
  expect_output(print(mvtb.cluster(covex)))
  expect_output(print(mvtb.ri(out)))
})
