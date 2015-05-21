#library(testthat)
#library(MASS)
#library(gbm)
#source("~/Documents/Projects/mvtboost/13/mvtboost_v14.R")
context("test_predict")

k <- 1
B <- matrix(c(1,1,rep(0,(k-1)*2)),nrow=2,ncol=k)
varx <- 2
vare <- diag(k)
diag(vare) <- 1:k
vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- .5

n <- 1000
X <- matrix(c(rnorm(n,1,.7),rbinom(n,2,.5)),ncol=2)
E <- mvrnorm(n,mu=rep(0,k),Sigma=vare)
Y <- X %*% B + E
colnames(X) <- paste0("X",1:varx)

context("predictions - mixed continuous categorical")

X2 <- factor(X[,2])
Xf <- data.frame(X1=X[,1],X2)
d1 <- data.frame(y=Y[,1],X1=X[,1],X2=X2)

m <- mvtb(X=Xf,Y=Y,n.trees=50,interaction.depth=3,shrinkage=.5,alpha=.5,weight.type=3,bag.frac=1,trainfrac=1,samp.iter=FALSE,compress=FALSE,cv.folds=1,s=1:1000)
m1 <- gbm(y~.,distribution="gaussian",n.trees=50,interaction.depth=3,data=d1,bag.fraction=1,train.fraction=1,shrinkage=.5)

expect_equal(m$models[[1]]$c.splits,m1$c.splits)
expect_equal(m1$initF,m$models[[1]]$initF)

for(i in c(1:10,20,50)) {
     p <- predict.mvtb(m,newdata=Xf,n.trees=i) 
     p1 <- predict(m1,newdata=Xf,n.trees=i)
     expect_equal(c(p),p1,info=paste0("ntrees = ",i))
 }


context("predictions - categorical")

n <- 1000
k <- 1
X <- matrix(rbinom(n*2,3,.5),ncol=2)
vare <- diag(k)
diag(vare) <- 1:k
vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- .5
E <- mvrnorm(n,mu=rep(0,k),Sigma=vare)
B <- matrix(c(1,1,rep(0,(k-1)*2)),nrow=2,ncol=k)
Y <- X %*% B + E
colnames(X) <- paste0("X",1:2)
X <- data.frame(X1=factor(X[,1]),X2=factor(X[,2]))
newdata <- X
d1 <- data.frame(y=Y[,1],X)

m <- mvtb(X=X,Y=Y,n.trees=50,shrinkage=.5,alpha=.05,weight.type=3,bag.frac=1,trainfrac=1,samp.iter=FALSE,compress=FALSE,cv.folds=1,s=1:1000)
m1 <- gbm(y~.,distribution="gaussian",n.trees=50,data=d1,bag.fraction=1,train.fraction=1,shrinkage=.5)


for(i in c(1,2,3,10,50)) {
    p <- predict.mvtb(m,newdata=newdata,n.trees=i) 
    p1 <- predict(m1,newdata=newdata,n.trees=i)
    expect_equal(p[,1,1],p1,info=paste0("ntrees = ",i))
}


context("two continuous")

n <- 1000
k <- 1
X <- matrix(rnorm(n*2,0,1),ncol=2)
vare <- diag(k)
diag(vare) <- 1:k
vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- .5
E <- mvrnorm(n,mu=rep(0,k),Sigma=vare)
B <- matrix(c(1,1,rep(0,(k-1)*2)),nrow=2,ncol=k)
Y <- X %*% B + E
colnames(X) <- paste0("X",1:2)
Xf <- data.frame(X1=X[,1],X2=X[,2])
newdata <- Xf
d2 <- data.frame(y=Y[,1],Xf)



set.seed(100)
m <- mvtb(X=Xf,Y=Y,n.trees=50,shrinkage=.5,interaction.depth=3,alpha=.05,weight.type=3,bag.frac=1,trainfrac=1,samp.iter=FALSE,compress=FALSE,cv.folds=1,s=1:1000)
m1 <- gbm(y~.,distribution="gaussian",n.trees=50,data=d2,interaction.depth=3,bag.fraction=1,train.fraction=1,shrinkage=.5)

for(i in c(1:10,20,50)) {
    p <- predict.mvtb(m,newdata=newdata,n.trees=i) 
    p1 <- predict(m1,newdata=newdata,n.trees=i)
    expect_equal(c(p),p1,info=paste0("ntrees = ",i))
}

context("test error")

# Test error
E2 <-  mvrnorm(n,mu=rep(0,k),Sigma=vare)
Y2 <- X %*% B + E2
dt <- data.frame(y=Y2[,1],Xf)

pe.wt <- function(mod,dt,bi=NULL,method=1){
    yhat <- predict.mvtb(mod,newdata=Xf,n.trees=50)
    mse <- mean((dt$y-yhat)^2)
    return(mse)
}
pe.gbm <- function(mod,dt){
    yhat <- predict(mod,newdata=Xf,n.trees=50)
    mse <- mean((dt$y-yhat)^2)
    return(mse)
}
    
expect_equal(pe.wt(m,dt=dt,method=3),pe.gbm(m1,dt=dt))


#source("../Simulations/4_Sim_nonlin/crc/dg_nonlin.R")
#d <- sparse.nl.dgen(k=1,lin=TRUE)
#df <- data.frame(y=d$Y,d$X)
#o1 <- gbm(y~.,data=df,distribution="gaussian",n.trees=100,shrinkage=.5,bag.frac=1)
#o2 <- mvtb(Y=d$Y,X=d$X,n.trees=100,shrinkage=.5,bag.frac=1)

#p <- cbind(predict(o1,newdata=data.frame(d$X),n.trees=100),c(predict.mvtb(o2,n.trees=100,newdata=data.frame(d$X))))
#expect_equal(p[,1],p[,2])

