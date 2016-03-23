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

## Again, tests just to make sure that they run
out <- mvtb(Y=Y,X=X,shrinkage=.1,n.trees=1000)
expect_output(summary(out),"trees")
expect_output(summary(out),"influence")
expect_equal(sum(summary(out,print=FALSE,relative="tot")$relative.influence),100)

expect_output(mvtb.covex(Y=Y,X=X,out=out),"")
expect_equal(dim(mvtb.covex(out)),c(ncovs,p))
expect_identical(mvtb.covex(out),out$covex)
expect_output(mvtb.cluster(out),"")
expect_equal(dim(mvtb.cluster(out)),c(ncovs,p))
mvtb.cluster(out,plot=TRUE)
expect_output(mvtb.cluster(out,dist.method="manhattan",clust.method="complete"),"")
#expect_identical(mvtb.cluster(out,clust.method=NULL),out$covex)

expect_output(mvtb.ri(out),"")
expect_equal(dim(mvtb.ri(out)),c(p,q))
expect_output(mvtb.ri(out,weighted = TRUE),"")
expect_equivalent(sum(mvtb.ri(out,relative = "tot")),100)
expect_equal(sum(colSums(mvtb.ri(out,relative = "col")))-q*100,0,tolerance=1E-12)
expect_output(mvtb.ri(out,relative = "n"),"")
expect_output(gbm.ri(out),"")
expect_output(print.mvtb(out),"")

context("test_summary_compression")

out <- mvtb(Y=Y,X=X,compress = T)
expect_output(summary(out),"")
expect_output(print(out),"")
expect_output(mvtb.cluster(out),"")
expect_output(gbm.ri(out),"")

