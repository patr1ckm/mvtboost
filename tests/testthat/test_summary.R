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

out <- mvtb(Y=Y, X=X, shrinkage=.1, n.trees=100)
outpc <- pcb(Y=Y, X=X, shrinkage=.1, n.trees=100)
out.comp <- mvtb(Y=Y, X=X, shrinkage=.1, n.trees=100, compress=TRUE)
mods <- list(out, outpc, out.comp)

test_that("summary",{ 
  ## Tests just to make sure that they run
  for(i in seq_along(mods)){
    expect_output(print(summary(mods[[i]])),"trees")
    expect_output(summary(mods[[i]]),"importance")
    suminf <- sum(summary(mods[[i]], print=FALSE, relative="tot")$relative.influence)
    expect_equal(suminf, 100)
  }
})

test_that("mvtb.cluster",{
  covex <- mvtb.covex(out, Y=Y, X=X)
  expect_output(print(mvtb.cluster(covex)))
  
  # test dimensions
  expect_equal(dim(mvtb.cluster(covex)), c(ncovs,p))
  
  # run plot
  cluster.covex <- mvtb.cluster(covex, plot=TRUE)
  
  cluster.covex <- mvtb.cluster(covex, dist.method="manhattan", 
                                clust.method="complete")
  expect_output(print(cluster.covex))
  
  # test clustering importances
  for(i in seq_along(mods)){
    expect_output(print(mvtb.cluster(importance(mods[[i]]))))
  }
})

test_that("importance",{
  for(i in seq_along(mods)){
    expect_output(print(importance(mods[[i]])))
    
    # dimensions
    expect_equal(dim(importance(mods[[i]])),c(p,q))
    
    # testthat sums of importance correctly equal 100 
    expect_equivalent(sum(importance(mods[[i]],relative = "tot")),100)
    
    # testthat sum of each column is 100
    expect_equal(sum(colSums(importance(mods[[i]],relative = "col")))-q*100,
                 0,tolerance=1E-12)
    
    # should produce raw
    expect_output(print(importance(mods[[i]],relative = "n")))
    
    # verifies that print.mvtb is being called
    expect_output(print(mods[[i]]), "List of ")
  }
})
  
