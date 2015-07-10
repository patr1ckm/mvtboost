## purpose: testing equivalence of stagewise

context("test_wt_stagewise")

fp <- paste0(getwd(),"/test_wt_stagewise.R")
#library(lars)
tol <- 1E-7
n <- 1000
mse <- function(x){mean((x-mean(x))^2)}
for(i in 1:4) {
    for(wt in 1:4) {
        context(paste0("Univariate Stagewise loss",i," wt = ", wt))
        test_that(paste0("1 var: loss ",i," wt = ", wt), {

            np <- 1
            set.seed(123)
            x1 <- rbinom(n,1,.5)
            b <- matrix(c(1,.4),nrow=np+1,ncol=1)
            e <- rnorm(n)
            X <- matrix(c(rep(1,n),x1),ncol=np+1,nrow=n)
            y <- X %*% b + e
            X <- matrix(c(x1),ncol=np)

            lfs <- lars::lars(x=X,y=y,type="forward.stagewise")
            ##lsw <- lars(x=X,y=y,type="stepwise")
            out <- mvtb(X=data.frame(X),Y=matrix(y),n.trees=4,shrinkage=1,trainfrac=1, weight.type=wt,cov.discrep=i,s=1:nrow(X))

            p1 <- c(predict.mvtb(out,newdata=data.frame(X),n.trees=out$maxiter))
            p2 <- predict(lfs,newx=X)$fit[,np+1]
            ##p3 <- predict(lsw,newx=X)$fit[,3] all(p2==p3)

            expect_true(mse(p1-p2) < 1E-6)
            #expect_true(all(abs(p1-p2) < tol))        
            expect_true(cor(p1,p2,method=) > .9999)
            expect_true(var(p1)/var(y)-var(p2)/var(y) < 1E-4) #r2
            expect_true(mse(y-p1)-mse(y-p2) < 1E-6)

            ## test 2 - are the predictors selected in the correct order?
            ord1 <- unlist(lfs$actions)
            ord2 <- out$bestxs[1:length(ord1)]
            expect_equal(ord2,1:np)
            expect_equal(ord1,ord2)
        })
    }
}

for(i in 1:4) {
    for(wt in 1:2) {

        context(paste0("Multivariate Stagewise loss",i," wt = ", wt))    
        test_that(paste0("3 vars: loss",i," wt = ", wt), {
            
            n <- 1000
            np <- 3
            set.seed(123)
            x1 <- rbinom(n,1,.5)
            x2 <- rbinom(n,1,.5)
            x3 <- rbinom(n,1,.5)
            b <- matrix(c(1,.4,.2,0),nrow=np+1,ncol=1)
            e <- rnorm(n)
            X <- matrix(c(rep(1,n),x1,x2,x3),ncol=np+1,nrow=n)
            y <- X %*% b + e
            X <- matrix(c(x1,x2,x3),ncol=np)
            y <- matrix(y)

            lfs <- lars::lars(x=X,y=y,type="forward.stagewise")
            out <- mvtb(X=X,Y=y,n.trees=4,shrinkage=1,trainfrac=1, weight.type=wt,cov.discrep=i,s=1:nrow(X))
            
            p1 <- c(predict.mvtb(out,newdata=data.frame(X),n.trees=out$maxiter))
            p2 <- predict(lfs,newx=X)$fit[,4]
          
            expect_true(mse(p1-p2) < 1E-6)
            expect_true(all(abs(p1-p2) < .0005))
            expect_true(cor(p1,p2,method="spearman") > .99)
            expect_true(cor(p1,p2,method=) > .9999)           

            ## test 2 - are the predictors selected in the correct order?
            ord1 <- unlist(lfs$actions)
            ord2 <- out$bestxs[1:length(ord1)]
            expect_equal(ord2,1:np)
            expect_equal(ord1,ord2)
            
        })
    }
}


## x <- matrix(rnorm(n),n,1)
## y <- x * .4 + e
## xc <- cut(x,breaks=c(-Inf,.5,Inf))
## summary(lm(y~x))
## summary(lm(y~xc))

## out$bestxs


## summ <- function(y,xc){
##     f.list <- list(length=length,var=var,mean=mean,ssq=ssq,mse=mse)
##     lapply(f.list,function(f){aggregate(x=y,by=list(xc),FUN=f)})
## }
## summ(y,xc)
## var(y)

## o <- gbm(y~xc,data=d,shrinkage=1,n.trees=1,distribution="gaussian",bag.fraction=1)
## pretty.gbm.tree(o)
## pretty.gbm.tree(o)$Prediction+o$initF
## ## The improvement is the reduction in sums of squares within each node.
## pretty.gbm.tree(o)$ErrorReduction[1] - (ssq(y)-sum(summ(y,xc)$ssq$V1)) < tol

