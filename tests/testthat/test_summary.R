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

## Again, tests just to make sure that they run
out <- mvtb(Y=Y,X=X)
summary(out)
summary(out,covex=FALSE)
summary(out,relative="tot",covex=FALSE)
mvtb.cluster(out)
mvtb.cluster(out,plot=TRUE)
mvtb.cluster(out,dist.method="manhattan",clust.method="complete")
#expect_identical(mvtb.cluster(out,clust.method=NULL),out$covex)

mvtb.ri(out)
mvtb.ri(out,weighted = TRUE)
mvtb.ri(out,relative = "tot")
mvtb.ri(out,relative = "col")
mvtb.ri(out,relative = "n")
gbm.ri(out)
print.mvtb(out)

out <- mvtb(Y=Y,X=X,compress = T)
summary(out)
print(out)
mvtb.cluster(out)
gbm.ri(out)

