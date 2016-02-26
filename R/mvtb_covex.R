#' Return the covariance explained matrix
#' 
#' @param object an object of class \code{mvtb}
#' @export
mvtb.covex <- function(Y,X,object,n.trees=NULL,iter.details=F) {
  
  if(any(unlist(lapply(object,function(li){is.raw(li)})))){
    object <- mvtb.uncomp(object)
  }
  if(class(Y) != "matrix") { Y <- as.matrix(Y) }
  if(is.null(ncol(X))){ X <- as.matrix(X)}
  stopifnot(nrow(X) == nrow(Y))
  
  n <- nrow(X); k <- ncol(Y); p <- ncol(X);
  
  # Get models
  if(is.null(n.trees)){n.trees <- min(object$best.trees)}
  finaltree <- list()
  for(m in 1:k) { 
    finaltree[[m]] <- object$models[[m]]$trees
  }
  
  ## Covex
  D <- Rm <- matrix(0,n,k)   
  Res.cov <- array(0,dim=c(k,k,k,n.trees))
  covex <- array(0,dim=c(p,k*(k+1)/2))
  rownames(covex) <- colnames(X)
  names <- outer(1:k,1:k,function(x,y){paste0(colnames(Y)[x],"-",colnames(Y)[y])})
  colnames(covex) <- names[lower.tri(names,diag=TRUE)]
  
  ## Influence
  wm.raw <- wm.rel <- matrix(0,nrow=n.trees,ncol=k)     #raw, relative
  rel.infl <- w.rel.infl <- array(0,dim=c(p,k,n.trees)) # influences at every iteration
  
  init <- unlist(lapply(object$models,function(m){m$initF}))
  D <- Rm <- matrix(0,n,k) 
  for(m in 1:k) { D[,m] <- Y[,m]-init[m] } # current residuals at iteration 1
  yhat <- array(c(rep(init,each=n),yhat),dim=c(n,k,n.trees+1))
  Dpred <- yhat[,,2:(n.trees+1),drop=F]-yhat[,,1:(n.trees),drop=F] # Dpred is the unique contribution of each tree
  bestxs   <- matrix(0,nrow=n.trees,ncol=k) 
  
  ## ss <- matrix(0,nrow=length(s),ncol=n.trees)
  
  ## 1. bootstrap covex if desired.
  ## for(i in 1:n.trees) {
  ##  if(samp.iter) {
  ##    ss[,i] <-  sample(s,length(s),replace=TRUE) # if replace = FALSE, this just permutes the rows.
  ##  } else {
  ##    ss[,i] <- s
  ##  }
  #}  
  
  
  # 2. Compute covariance discrepancy
  for(i in 1:(n.trees)) {        
    
    ## 2.1 From each model get the stuff we need at the current iteration
    for(m in 1:k) {            
      ## 1.2 For each model compute predicted values and influence
      tree.i <- finaltree[[m]][i]
      
      ## 1.3 Replace mth outcome with its residual, compute covariance           
      Rm <- D
      Rm[,m] <- D[,m]-Dpred[,m,i]
      Res.cov[,,m,i] <- cov(as.matrix(Rm),use="pairwise.complete") 
      #Res.cov[,,m,i] <- cov(D) - cov(D[,m]-Dpred[,m,i])
      ## 2. Evaluate loss criteria on training sample. Covariance reduced, correlation reduced, uls, or gls
      wm.raw[i,m] <- eval.loss(Rm=as.matrix(Rm[s,]),D=as.matrix(D[s,]),alpha=alpha,type=cov.discrep)
    }              
    
    wm.raw[i,wm.raw[i,] < 0] <- 0 # truncate negative weights to 0
    
    ## 3. Check early stopping criteria.
    if(all(wm.raw[i,]==0)) {
      wm.rel[i,] <- 0
      final.iter <- TRUE; # break at the end of the loop, so test and train error can still be computed.
    }                     
    
    ## 4. Compute weight types (relative, 0-1, and none)
    if(!final.iter) {
      if(weight.type==1) {
        ##print("relative weights")
        wm.rel[i,] <- wm.raw[i,]/sum(wm.raw[i,]) # relative weights.
      } else if(weight.type==2) {
        wm.rel[i,which.max(wm.raw[i,])] <- 1 # 0-1 weights for covariance explained loss functions (want to maximize)       
      } else {
        ##print("equal weight")
        wm.rel[i,] <- rep(1,k) # equal weight
      }
    }
    # compute the best mod, and which xs were selected
    #besty <- bestys[i] <- which.max(wm.raw[i,])           
    ## as an approximation, choose the best x by relative influence in each col at iteration i
    bestxs[i,] <- bestx <- apply(rel.infl[,,i,drop=F],2,function(col){which.max(col)})    
    
    # compute the covariance reduced (explained) by the best predictor for each outcome
    correc <- 1+(1-shrinkage)
    for(k in 1:m) {
      Sd <- cov(D)-Res.cov[,,k,i]
      Sd[lower.tri(Sd)] <- Sd[lower.tri(Sd)]*correc
      if(k > 1) {
        if(bestx[k] == bestx[k-1]) {
          # the covariance elements will be counted twice. Don't want this to happen.
          # since the covariance elements at column (row) k are already counted, set them to zero for k+1
          # at k = m, this will only update the variance
          Sd[1:(k-1),] <- 0
          Sd[,1:(k-1)] <- 0
        }
      }
      covex[bestx[k],] <- covex[bestx[k],] + Sd[lower.tri(Sd,diag=TRUE)]     
      
    }
    
    #Rm <- D - Dpred[,,i]
    D <- Y - yhat[,,i+1] # i=1 is init (colMeans Y)
    
  }
  covex <- t(covex)
  if(iter.details){
    
    ## Todo to make tests pass probably
    #fl <- list(covex=covex, weights=weights,wm.raw=matrix(wm.raw[1:i,,drop=FALSE],nrow=i,ncol=k),wm.rel=wm.rel[1:i,,drop=FALSE])
  } else {
    return(covex)
  }
}

#' @importFrom stats cov
eval.loss <- function(Rm,D,alpha,type) {
  ## Rm = n x k matrix of residuals
  ## D = n x k matrix of Ytilde at iteration i
  ## S = k x k covariance matrix of Ytilde at iteration i
  ## Res.cov = k x k covariance matrix of Rm (residuals)
  ## Sd = S - Res.cov, k x k the discrepancy between S and Res.cov. The larger the discrepancy, the more covariance is explained.
  Res.cov <- cov(as.matrix(Rm),use="pairwise.complete")
  S <- cov(as.matrix(D),use="pairwise.complete")
  Sd <- S-Res.cov
  ## type is an integer 1:5, maps to loss function 1, 2, ... etc.
  switch(type,
         cov.red(S=S,Res.cov=Res.cov,alpha=alpha),
         cor.red(Rm=Rm,D=D),
         uls(Sd=Sd),
         #gls(Sd=Sd,Sinv=Sinv),
         detcov(S,Res.cov))
}

detcov <- function(S,Res.cov) {
  det(S)-det(Res.cov)
}

## Covariance explained, weighted by alpha. alpha*diagonal (1-alpha)*off-diagonal
cov.red <- function(S,Res.cov,alpha) {
  Sd <- abs(S-Res.cov)
  alpha*sum(diag(Sd),na.rm=TRUE) + (1-alpha)*sum(Sd[lower.tri(Sd)],na.rm=TRUE)
}

#' @importFrom stats cor
cor.red <- function(Rm,D) {    
  if(ncol(Rm) == 1) {
    wm <- 1-cor(Rm,D,use="pairwise.complete")
  } else {
    Sd <- abs(cor(D,use="pairwise.complete") - cor(Rm,use="pairwise.complete"))
    wm <- sum(Sd[upper.tri(Sd)],na.rm=TRUE)
  }
  return(wm)
}

## Some covariance discrepancy functions, motivated from SEM.
## 1/2 the ||R||_F, the frobenius norm of the  covariance discrepancy matrix
uls <- function(Sd) { .5 * sum(diag(Sd %*% t(Sd))) }
