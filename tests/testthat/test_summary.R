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
out2 <- mvtb_sep(Y=Y, X=X, shrinkage=.1, n.trees=100)
outpc <- pcb(Y=Y, X=X, shrinkage=.1, n.trees=100)
out.comp <- mvtb(Y=Y, X=X, shrinkage=.1, n.trees=100, compress=TRUE)
out2.comp <- mvtb_sep(Y=Y, X=X, shrinkage=.1, n.trees=100, compress=TRUE)
outpc.comp <- mvtb_sep(Y=Y, X=X, shrinkage=.1, n.trees=100, compress=TRUE)
mods <- list(out, out2, outpc, out.comp, out2.comp, outpc.comp)

test_that("summary",{ 
  ## Tests just to make sure that they run
  for(i in seq_along(mods)){
    expect_output(print(summary(mods[[i]])),"trees")
    expect_output(summary(mods[[i]]),"influence")
    suminf <- sum(summary(mods[[i]], print=FALSE, relative="tot")$relative.influence)
    expect_equal(suminf, 100)
  }
})

test_that("mvtb.cluster",{
  covex <- mvtb.covex(out, Y=Y,X=X)
  expect_output(print(mvtb.cluster(covex)))
  
  # test dimensions
  expect_equal(dim(mvtb.cluster(covex)),c(ncovs,p))
  
  # run plot
  expect_output(mvtb.cluster(covex,plot=TRUE))
  
  cluster.covex <- mvtb.cluster(covex,dist.method="manhattan",clust.method="complete")
  expect_output(print(cluster.covex))
  
  # test clustering influences
  for(i in seq_along(mods)){
    expect_output(print(mvtb.cluster(influence(mods[[i]]))))
  }
})

test_that("influence",{
  for(i in seq_along(mods)){
    expect_output(print(influence(mods[[i]])))
    
    # dimensions
    expect_equal(dim(influence(mods[[i]])),c(p,q))
    
    # testthat sums of influence correctly equal 100 
    expect_equivalent(sum(influence(mods[[i]],relative = "tot")),100)
    
    # testthat sum of each column is 100
    expect_equal(sum(colSums(influence(mods[[i]],relative = "col")))-q*100,0,tolerance=1E-12)
    
    # should produce raw
    expect_output(print(influence(mods[[i]],relative = "n")))
    
    # verifies that print.mvtb is being called
    expect_output(print(mods[[i]]), "List of ")
  }
})
  
