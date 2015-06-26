
## ps <- c(10,50)
## niters <- c(10,100)
## lps <- list()
## res <- cbind(expand.grid(niters=niters,ps=ps),time=rep(0,length(ps)*length(niters)))

## context("speed")
## ## set.seed(123)
## ## for(i in 1:nrow(res)) {
## ##     d <- sparse.dgen(p=res$ps[i],numpred=3)
## ##     res$time[i] <- system.time(mvtb(X=d$X,Y=d$Y,niter=res$niters[i],alpha=.5,mv.shrink=.01,trainfrac=1,samp.iter=FALSE,loss.function=1,weight.type=2))[1]
## ## }


## #print(res)

