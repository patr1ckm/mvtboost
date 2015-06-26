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
cluster.covex(out)
mvtb.ri(out)
print.mvtb(out)
