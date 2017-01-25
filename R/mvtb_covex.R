#' Estimate the covariance explained matrix
#' 
#' @param object an object of class \code{mvtb}
#' @param Y vector, matrix, or data.frame for outcome variables with no missing values. To easily compare influences across outcomes and for numerical stability, outcome variables should be scaled to have unit variance.
#' @param X vector, matrix, or data.frame of predictors. For best performance, continuous predictors should be scaled to have unit variance. Categorical variables should converted to factors.
#' @param n.trees number of trees to use. Defaults to the minimum number of trees by CV, test, or training error, and the maximum of this for all outcomes.
#' @param iter.details \code{TRUE/FALSE}. Return the loss, relative loss, and selected predictors at each iteration as a list
# #' @param cov.discrep Norm of the covariance discrepancy (see details)
# #' @param alpha The weight given to covariance explained relative to variance explained in the loss function (see details)
# #' @param weight.type see details
#' @export
#' @return Covariance explained matrix, or a list if \code{iter.details} is \code{TRUE}.
#' @importFrom stats cov
mvtb.covex <- function(object,Y,X,n.trees=NULL,iter.details=FALSE) {
  cov.discrep=1; alpha=.5; weight.type=1; # these are no longer arguments for simplicity.
  
  if(any(unlist(lapply(object,function(li){is.raw(li)})))){
    object <- mvtb.uncomp(object)
  }
  if(class(Y) != "matrix") { Y <- as.matrix(Y) }
  if(is.null(ncol(X))){ X <- as.matrix(X)}
  
  stopifnot(nrow(X) == nrow(Y))
  if(alpha > 1 | alpha < 0){ stop("alpha should be > 0, < 1")}
  if(any(is.na(Y))){ stop("NAs not allowed in outcome variables.")}
  
  if(is.null(n.trees)) { n.trees <- max(apply(object$best.trees, 1, min, na.rm=T))}

  n <- nrow(X); k <- ncol(Y); p <- ncol(X);
  
  # Get models

  finaltree <- list()
  for(m in 1:k) { 
    finaltree[[m]] <- object$models[[m]]$trees
  }
  shrinkage <- object$params$shrinkage
  
  ## Covex
  D <- Rm <- matrix(0,n,k)   
  Res.cov <- array(0,dim=c(k,k,k, n.trees))
  covex <- array(0,dim=c(p,k*(k+1)/2))
  rownames(covex) <- colnames(X)
  names <- outer(1:k,1:k,function(x,y){paste0(object$ynames[x], "-", object$ynames[y])})
  colnames(covex) <- names[lower.tri(names,diag=TRUE)]
  
  ## Loss function evaluations
  ## raw loss, relative loss, overall loss compared to cov(Y)
  wm.raw <- wm.rel <- matrix(NA,nrow=max(n.trees),ncol=k)
  
  ## Influence
  rel.infl <- w.rel.infl <- array(0,dim=c(p,k,n.trees)) # influences at every iteration
  bestxs   <- matrix(0,nrow=max(n.trees),ncol=k) 
  
  
  init <- unlist(lapply(object$models,function(m){m$initF}))
  D <- Rm <- yhatm <- matrix(0,n,k) 
  for(m in 1:k) { D[,m] <- Y[,m]-init[m] } # current residuals at iteration 1
  yhat <- predict.mvtb.array(list(models=object$models),n.trees=1:n.trees,newdata=X,drop=FALSE)
  yhat <- array(c(rep(init,each=n),yhat),dim=c(n,k,n.trees+1))
  Dpred <- yhat[,,2:(n.trees+1),drop=F]-yhat[,,1:(n.trees),drop=F] # Dpred is the unique contribution of each tree
  s <- object$s
  final.iter <- FALSE
  
  # 1. Compute covariance discrepancy
  for(i in 1:(n.trees)) {        
    
    ## From each model get the stuff we need at the current iteration
    for(m in 1:k) {            
      ## For each model compute predicted values and influence
      tree.i <- finaltree[[m]][i]
      rel.infl[,m,i] <- importance_from_tree_list(tree.i, n.trees=1, object$xnames)
      
      ## Replace mth outcome with its residual, compute covariance           
      Rm <- D
      Rm[,m] <- D[,m]-Dpred[,m,i]
      Res.cov[,,m,i] <- cov(as.matrix(Rm),use="pairwise.complete") 
      
      ## Evaluate loss criteria on training sample. Covariance reduced, correlation reduced, uls, or gls
      wm.raw[i,m] <- eval.loss(Rm=as.matrix(Rm[s,]),D=as.matrix(D[s,]),alpha=alpha,type=cov.discrep)
      
    }              
    
    wm.raw[i,wm.raw[i,] < 0] <- 0 # truncate negative weights to 0

    
    ## 4. Compute weight types (relative, 0-1, and none)
    if(!final.iter) {
      if(weight.type==1) {
        wm.rel[i,] <- wm.raw[i,]/sum(wm.raw[i,]) # relative weights.
      } else if(weight.type==2) {
        wm.rel[i,which.max(wm.raw[i,])] <- 1 # 0-1 weights for covariance explained loss functions (want to maximize)       
      } else {
        wm.rel[i,] <- rep(1,k) # equal weight
      }
    }
       
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
    
    D <- Y - yhat[,,i+1] # i=1 is init (colMeans Y)
    
  }
  covex <- t(covex)
  if(iter.details){
    fl <- list(covex=covex,loss=matrix(wm.raw[1:i,,drop=FALSE],nrow=i,ncol=k),rel.loss=wm.rel[1:i,,drop=FALSE],bestxs=bestxs)
  } else {
    fl <- covex
  }
  return(fl)
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

predict.mvtb.array <- function(object, newdata, n.trees, drop=TRUE, ...) {
  
  if(any(unlist(lapply(object,function(li){is.raw(li)})))){
    object <- mvtb.uncomp(object)
  }
  
  K <- length(object$models)
  treedim <- ifelse(length(n.trees) > 1,max(n.trees),1)
  Pred <- array(0,dim=c(nrow(newdata),K,treedim))  
  for(k in 1:K) {                                     
    p <- rep(0,nrow(newdata))        
    p <- predict(object$models[[k]],n.trees=n.trees, newdata=data.frame(newdata))
    Pred[,k,] <- p
  }
  
  if(drop){
    Pred <- drop(Pred)
  }
  return(Pred)
}
