## Author: Patrick Miller
## RRV: 6-18
## Purpose: verify that weighted tree boosting with one dichotomized predictor and no shrinkage is equivalent to simple linear regression

context("test_wt_lr")


for(i in 1:4) {
    
context(paste0("Univariate loss",i))
test_that("single tree from weighted.trees is equivalent to linear regression", {
    
x <- rbinom(1000,1,.5)
xf <- factor(x)
b <- .4
e <- rnorm(1000)
y <- matrix(b*x+e)
X <- matrix(c(rep(1,1000),x),nrow=1000,ncol=2)
est.coef <- coef(summary(lm.out <- lm(y~x)))[1:2,1]
r <- resid(lm.out)
#X %*% est.coef


#out <- weighted.trees(X=matrix(x),Y=matrix(y),n.trees=1,mv.shrink=1,shrinkage=1,trainfrac=1, weight.type=1,cov.discrep=i,s=1:1000,bag.frac=1)
out <- mvtb(X=as.data.frame(x=xf),Y=y,n.trees=1,interaction.depth=1,shrinkage=1,trainfrac=1, weight.type=1,cov.discrep=i,s=1:1000,bag.frac=1)
out.gbm <- gbm::gbm(y~x,distribution="gaussian",data=data.frame(y=y,x=xf),n.trees=1,shrinkage=1,bag.frac=1,interaction.depth=1)

tol <- 1E-9

## test 1: mean of the squared residuals should equal the trainner
expect_true(abs(mean(r^2)-out$trainerr) < tol)
## test 2: test of loss function 1, should be the reduction in variance*alpha
if(i==1) {
    expect_true((var(y)-var(r))*out$params$alpha-out$wm.raw < tol)
}


#(var(y)-var(r))/var(y)

## test 3: fit and predicted values from tree
#m1 <- out$iter.models[[1]][[1]]
#p1 <- m1$fit+out$init
p2 <- predict(lm.out)
p3 <- predict(out.gbm,newdata=data.frame(x=xf),n.trees=1)
p4 <- predict.mvtb(out,newdata=data.frame(x=xf),n.trees=1)

mse <- function(y,yhat){sum((y-yhat)^2)}
P <- data.frame(lm=p2,gbm=p3,mvt=p4)
apply(P,2,mse,y=y)


expect_true(all(abs(p4)-abs(predict(lm.out)) < tol))
})
}

