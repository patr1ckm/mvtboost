
## goal here is to test that the selected Xs and Ys are reasonable given the loss functions. Mostly as a sanity check...
context("test_mvloss_shrink")

set.seed(123)
## Test 1
## x1 causes 2 outcomes to covary, x2 effects 1 outcome. expect that x1 should always be selected in the first iteration for Y1 or Y2
n <- 1000
B <- matrix(0,nrow=2,ncol=4)
B[1,1:2] <- 1
B[2,3] <- .3
E <- MASS::mvrnorm(n,rep(0,4),diag(4))
X <- matrix(rbinom(n*2,size=1,prob=.5),n,2) # use binomial variables so that the decision tree can capture the full effect.
#X <- mvrnorm(n,rep(0,nrow(B)),diag(nrow(B)))
Y <- X %*% B + E

## Test 2: 2 bivariate effects
B2 <- matrix(0,nrow=3,ncol=4)
B2[3,1:2] <- 2
B2[2,2:3] <- 1
B2[1,1] <- .5
X2 <- matrix(rbinom(n*nrow(B2),size=1,prob=.5),n,nrow(B2))
#X2 <- mvrnorm(n,rep(0,nrow(B2)),diag(nrow(B2)))
Y2 <- X2 %*% B2 + E

shrinkage <- c(1,.5)

## these are just some reasonable order tests about which Ys and Xs are selected
for(j in 1:2) {
  for(i in 1:4) {
    for(wt in 1:2) {
      shrinking <- ((wt==1) || (shrinkage[j] < 1))                               
      context(paste0("loss = ",i, " wt = ",wt, " shrinkage = ", shrinkage[j]))
      test_that("expectations for 1 iteration, scenario 1", {
        r <- mvtb(X=X,Y=Y,cov.discrep=i,n.trees=1,shrinkage=shrinkage[j],alpha=.5,trainfrac=1,weight.type=wt)
        expect_equal(r$bestxs,1,info=paste(r$bestxs, collapse=" " )) # predictor one should be selected
        expect_true(r$bestys %in% 1:2,info=paste(r$bestys, collapse=" " )) # residuals one or two should be replaced
        expect_true(all(r$wm[1:2] > r$wm[3:4]),info=paste(r$wm, collapse=" " ))
        expect_false(4 %in% r$bestys,info=paste(r$bestys, collapse=" " )) # outcome 4 should never be selected
      })
      if(j == 1) {
        test_that("expectations for 3 iterations, scenario 2, no shrinkage", {
          r <- mvtb(X=X2,Y=Y2,cov.discrep=i,n.trees=3,shrinkage=shrinkage[j],alpha=.5,trainfrac=1,weight.type=wt)
          expect_true(all(r$bestmod[1:2] %in% 1:2))
          expect_equal(r$bestxs[1:2],c(3,2))
        })
        
        break; # don't run the other tests here
      }
      
      if(shrinking) {
        test_that("expectations for 3 iterations, scenario 1, shrinkage", {
          r <- mvtb(X=X,Y=Y,cov.discrep=i,n.trees=3,shrinkage=shrinkage[j],alpha=.5,trainfrac=1,weight.type=wt)
          expect_true(all(r$bestxs[1:2] %in% 1),info=paste(r$bestxs, collapse=" " ))
          ## if wp=F, then on the first iteration, the effects of X1 will be removed from Y1 and Y2. So we can only say anything about the first iteration
          expect_true(all(r$bestys[1] %in% 1:2),info=paste(r$bestys, collapse=" " ))                        
          expect_false(4 %in% r$bestys,info=paste(r$bestys, collapse=" " )) # outcome 4 should never be selected 
        })
        
        test_that("expectations for 5 iterations, scenario 2, shrinkage", {
          r <- mvtb(X=X2,Y=Y2,cov.discrep=i,n.trees=5,shrinkage=shrinkage[j],alpha=.5,trainfrac=1,
                    weight.type=wt)
          expect_equal(r$bestxs[1],3,info=paste(r$bestxs, collapse=" " ))
          #expect_true(all(r$bestxs[1:4] %in% 2:3),info=paste(r$bestxs, collapse=" " )) # doesn't make sense when wp=FALSE
          expect_true(all(r$bestys[1] %in% 1:2),info=paste(r$bestys, collapse=" " ))
          expect_true(all(r$bestys[1:5] %in% 1:3),info=paste(r$bestys, collapse=" " ))                   
          expect_false(4 %in% r$bestys)
        })
        
      } else {
        ## test_that("expectations for 3 iterations, scenario 1, no shrinkage", {
        ##     r <- mvtb(X=X,Y=Y,cov.discrep=i,n.trees=3,shrinkage=shrinkage[j],alpha=.5,trainfrac=1,
        ##                         weight.type=wt,)
        ##     expect_equal(r$bestxs, c(1,1,2),info=paste(r$bestxs, collapse=" " ))
        ##     expect_equal(r$bestys[3],3,info=paste(r$bestys, collapse=" " ))
        ##     expect_true(all(r$bestys[1:2] %in% 1:2),info=paste(r$bestys, collapse=" " ))
        ##     expect_false(4 %in% r$bestys,info=paste(r$bestys, collapse=" " ))
        ## })
        
        ## test_that("expectations for 5 iterations, scenario 2, no shrinkage", {
        ##     r <- mvtb(X=X2,Y=Y2,cov.discrep=i,n.trees=5,shrinkage=shrinkage[j],alpha=.5,trainfrac=1,
        ##                         weight.type=wt,)
        ##     #expect_equal(r$bestxs,c(3,3,2,2,1),info=paste(r$bestxs, collapse=" " ))
        ##     expect_equal(r$bestxs[1],3,info=paste(r$bestxs, collapse=" " ))
        ##     #expect_equal(r$bestxs[3],2,info=paste(r$bestxs, collapse=" " )) # doesn't make sense when wp=FALSE
        ##     expect_true(all(r$bestxs[1:2] %in% 2:3),info=paste(r$bestxs, collapse=" " ))
        ##     expect_true(all(r$bestys[1] %in% 1:2),info=paste(r$bestys, collapse=" " ))
        ##     expect_true(all(r$bestys[1:5] %in% 1:3),info=paste(r$bestys, collapse=" " ))
        ##     expect_equal(r$bestxs[5],1,info=paste(r$bestys, collapse=" " )) 
        ##     expect_false(4 %in% r$bestys)
        ## })
      }
    }
  }
}



#fp <- paste0(getwd(),"/test_mvloss_shrink.R")

#cov(Y)
#cov(X %*% B)

## r <- mvtb(X=X,Y=Y,n.trees=1,mv.shrink=1,shrinkage=1,alpha=.5,cov.discrep=1,trainfrac=1,weight.type=2)
## r1 <- mvtb(X=X,Y=Y,n.trees=1,mv.shrink=.5,shrinkage=1,alpha=.5,cov.discrep=1,trainfrac=1,weight.type=2)
## r2 <- mvtb(X=X,Y=Y,n.trees=1,shrinkage=.5,alpha=.5,cov.discrep=1,trainfrac=1,weight.type=2)
## r3 <- mvtb(X=X,Y=Y,n.trees=2,mv.shrink=.5,shrinkage=1,alpha=.5,cov.discrep=1,trainfrac=1,weight.type=2
#                     )
