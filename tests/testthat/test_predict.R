context("test_predict")

k <- 1
B <- matrix(c(1,1,rep(0,(k-1)*2)),nrow=2,ncol=k)
varx <- 2
vare <- diag(k)
diag(vare) <- 1:k
vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- .5

n <- 1000
X <- matrix(c(rnorm(n,1,.7),rbinom(n,2,.5)),ncol=2)
E <- MASS::mvrnorm(n,mu=rep(0,k),Sigma=vare)
Y <- X %*% B + E
colnames(X) <- paste0("X",1:varx)
X2 <- factor(X[,2])
Xf <- data.frame(X1=X[,1],X2)
d1 <- data.frame(y=Y[,1],X1=X[,1],X2=X2)

m <- mvtb(X=Xf,Y=Y,n.trees=50,interaction.depth=3,shrinkage=.5,bag.fraction=1,train.fraction=1,compress=FALSE,cv.folds=1,s=1:1000)
m1 <- gbm::gbm(y~.,distribution="gaussian",n.trees=50,interaction.depth=3,data=d1,bag.fraction=1,train.fraction=1,shrinkage=.5)
m2 <- mvtb_sep(Y=Y, X=Xf, n.trees=50, interaction.depth=3, bag.fraction=1, compress=F, cv.folds=1, s=1:1000, shrinkage=.5)

test_that("predictions - mixed continuous categorical", {
  expect_equal(m$models[[1]]$c.splits,m1$c.splits)
  expect_equal(m2$models[[1]]$c.splits,m1$c.splits)
  expect_equal(m1$initF,m$models[[1]]$initF)
  expect_equal(m1$initF,m2$models[[1]]$initF)
  
  for(i in c(1:10,20,50)) {
       p <- predict(m,newdata=Xf,n.trees=i) 
       p1 <- predict(m1,newdata=Xf,n.trees=i)
       p2 <- predict(m2,newdata=Xf,n.trees=i)
       expect_equal(c(p),p1,info=paste0("ntrees = ",i))
       expect_equal(c(p2),p1,info=paste0("ntrees = ",i))
  }
})

n <- 1000
k <- 1
X <- matrix(rbinom(n*2,3,.5),ncol=2)
vare <- diag(k)
diag(vare) <- 1:k
vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- .5
E <- MASS::mvrnorm(n,mu=rep(0,k),Sigma=vare)
B <- matrix(c(1,1,rep(0,(k-1)*2)),nrow=2,ncol=k)
Y <- X %*% B + E
colnames(X) <- paste0("X",1:2)
X <- data.frame(X1=factor(X[,1]),X2=factor(X[,2]))
newdata <- X
d1 <- data.frame(y=Y[,1],X)

m <- mvtb(X=X,Y=Y,n.trees=50,shrinkage=.5,bag.fraction=1,train.fraction=1,compress=FALSE,cv.folds=1,s=1:1000)
m1 <- gbm::gbm(y~.,distribution="gaussian",n.trees=50,data=d1,bag.fraction=1,train.fraction=1,shrinkage=.5)
m2 <- mvtb_sep(X=X,Y=Y,n.trees=50,shrinkage=.5,bag.fraction=1,train.fraction=1,compress=FALSE,cv.folds=1,s=1:1000)

test_that("predictions - categorical", {
  for(i in c(1,2,3,10,50)) {
      p <- predict(m,newdata=newdata,n.trees=i) 
      p1 <- predict(m1,newdata=newdata,n.trees=i)
      p2 <- predict(m2,newdata=newdata,n.trees=i) 
      expect_equal(p,p1,info=paste0("ntrees = ",i))
      expect_equal(p2,p1,info=paste0("ntrees = ",i))
  }
})

n <- 1000
k <- 1
X <- matrix(rnorm(n*2,0,1),ncol=2)
vare <- diag(k)
diag(vare) <- 1:k
vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- .5
E <- MASS::mvrnorm(n,mu=rep(0,k),Sigma=vare)
B <- matrix(c(1,1,rep(0,(k-1)*2)),nrow=2,ncol=k)
Y <- X %*% B + E
colnames(X) <- paste0("X",1:2)
Xf <- data.frame(X1=X[,1],X2=X[,2])
newdata <- Xf
d2 <- data.frame(y=Y[,1],Xf)

set.seed(100)
m <- mvtb(X=Xf,Y=Y,n.trees=50,shrinkage=.5,interaction.depth=3,bag.fraction=1,train.fraction=1,compress=FALSE,cv.folds=1,s=1:1000)
m1 <- gbm::gbm(y~.,distribution="gaussian",n.trees=50,data=d2,interaction.depth=3,bag.fraction=1,train.fraction=1,shrinkage=.5)
m2 <- mvtb_sep(X=Xf,Y=Y,n.trees=50,shrinkage=.5,interaction.depth=3,bag.fraction=1,train.fraction=1,compress=FALSE,cv.folds=1,s=1:1000)

test_that("two continuous", {
  for(i in c(1:10,20,50)) {
      p <- predict(m,newdata=newdata,n.trees=i) 
      p1 <- predict(m1,newdata=newdata,n.trees=i)
      p2 <- predict(m2,newdata=newdata,n.trees=i)
      expect_equal(c(p),p1,info=paste0("ntrees = ",i))
      expect_equal(c(p2),p1,info=paste0("ntrees = ",i))
  }
})



