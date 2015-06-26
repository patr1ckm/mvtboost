## purpose: test shrinkage. 

context("test_shrinkage")

tol <- 1E-12
n <- 1000
np <- 3
set.seed(123)
e <- rnorm(n)
X <- matrix(rbinom(n*np,1,.5),n,np)
b <- c(.8,.1,0)
Y <- X %*% b + e
X <- data.frame(X)
shrink <- .5

g <- gbm::gbm(Y~.,distribution="gaussian",data=data.frame(X=X,Y=Y),n.trees=5,shrinkage=shrink,interaction.depth=1,bag.fraction=1)

## approach 1 - shrunken effects -> choose predictors ->
## basically, with one outcome, this should reproduce a gbm results directly, and results using LM
#source("~/Documents/play/boost_ls_tree.R")
for(j in 1:4) {
    context(paste("uv.shrinkage-order, loss ",j))
    shrink <- .5
    out <- mvtb(X=X,Y=matrix(Y),n.trees=5,shrinkage=shrink,trainfrac=1, weight.type=2,cov.discrep=1,s=1:nrow(X))
    
    test_that("variables selected in order", {
        expect_true(all(out$bestxs[1:3] == c(1,1,1)))
        expect_true(all(out$bestxs[4:5] %in% c(1,2)))
        expect_true(!any(out$bestxs == 3))
        expect_true(all(out$bestxs == unlist(lapply(g$trees,function(t){t[[1]][[1]][[1]]}))+1))
    })

    context(paste("uv.shrinkage-size, loss ",j))
    test_that("correct effect sizes are estimated", {
        init <- mean(Y)
        r <- Y-init
        p2 <- p1 <- 0
        for(i in 1) {
            p <- 0
            if(i == 5) {
                o <- lm(r~X[,2])
            } else {
                o <- lm(r~X[,1])
            }
            p <- predict(out$models[[1]],n.trees=i,newdata=X)
            if(i == 1) {
                p2 <- p2+predict(o)*shrink+init
            #    p1 <- p1+p+init
            } else {
                p2 <- p2+predict(o)*shrink
            #    p1 <- p1+p
            }
            p3 <- predict(g,n.trees=i)
            expect_true(all(abs(p-p2) < tol),info=paste("iter",i),label=" same lr prediction")
            expect_true(all(abs(p-p3) < tol),info=paste("iter",i),label=" same gbm prediction")
                                        
            bhat <- diff(unique(p))
            expect_true(bhat-coef(o)[2]*shrink < tol,info=paste("iter",i),label="same lr effect size")                                        
            r <- r - predict(o)*shrink
        }
    })

}

## for weight.type == 2, mv.shrink and shrinkage are the same thing but occur at different points in the algorithm.
## for shrinkage, the shrinkage occurs after the first q fits, before the predictor is selected.
## for mv.shrink, the shrinkage occurs after the predictor, after the predictor is selected.
for(wp in c(TRUE,FALSE)) {
    for(j in 1:4) {
        context(paste("shrink-compare, loss ",j," wp: ",wp))        
        out2 <- mvtb(X=X,Y=matrix(Y),n.trees=5,shrinkage=.5,trainfrac=1, weight.type=2,cov.discrep=j,s=1:nrow(X))
        p2 <- predict.mvtb(out2,newdata=X,n.trees=out2$maxiter)
        p3 <- predict(g,n.trees=5)
        
        expect_true(all(abs(p2-p3) < tol),info="uv.shrink == gbm")
    }
}
