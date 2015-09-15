## Author: Patrick Miller
## Purpose: Multiple response (or multivariate) boosting with decision trees (stumps). Also multi-task, or multi-objective. Only for continuous Y.
## RRV: 9/3/2015
  

#' Fitting a Multivariate Tree Boosting Model
#'
#' Builds on \code{gbm} (Ridgeway 2013; Friedman, 2001) to fit a univariate tree model for each outcome, selecting predictors at each iteration that explain covariance between the outcomes. The number of trees included in the model can be chosen by minimizing the multivariate mean squared error using cross validation or a test set.
#'
#' @param X matrix or data.frame of predictors. For best performance, continuous predictors should be scaled to have unit variance. Categorical variables should converted to factors.
#' @param Y vector, matrix, or data.frame for outcome variables. Missing values must be imputed. To easily compare influences across outcomes and for numerical stability, outcome variables should be scaled to have unit variance.
#' @param n.trees maximum number of trees to be included in the model. Trees are grown until a minimum number observations in each node is reached. This minimum can be modified using additional arguments (below).
#' @param shrinkage a contant multiplier for the predictions from each tree to ensure a slow learning rate. Default is .01. Small shrinkage values may require a large number of trees to provide adequate fit.
#' @param interaction.depth fixed depth of trees to be included in the model. A tree depth of 1 corresponds to fitting stumps (main effects only), higher tree depths capture higher order interactions.
#' @param bag.frac   proportion of the training sample used to fit univariate trees for each response at each iteration. Default: 1
#' @param cv.folds   number of cross validation folds. Default: 1. Runs k + 1 models, where the k models are run in parallel and the final model is run on the entire sample. If larger than 1, the number of trees that minimize the multivariate MSE averaged over k-folds is reported in object$best.trees.
#' @param trainfrac  proportion of the sample used for training the multivariate additive model. The number of trees that minimize the multivariate MSE in the test set are given in out$best.trees. If both cv.folds and trainfrac are specified, the CV is done within the training set only.
#' @param s vector of indices denoting observations to be used for the training sample. If s is specified, trainfrac is ignored.
#' @param seednum integer to ensure that results are reproducible
#' @param compress T/F. Compress output results list using bzip2 (approx 10\% of original size). Default FALSE.
#' @param save.cv  T/F. Save all k-fold cross-validation models.
#' @param mc.cores Number of cores for cross validation.
#' @param alpha optional argument to select predictors that explain more variance or covariance in outcomes. Variance reductions are weighted by alpha, and covariance reductions are weighted by 1-alpha.
#' @param weight.type Experimental.
#' @param cov.discrep Experimental. Choose the type of covariance discrepancy.
#' @param samp.iter   Experimental.
#' @param ... additional arguments passed to \code{gbm}. These include \code{distribution}, \code{weights}, \code{var.monotone}, \code{n.minobsinnode}, \code{keep.data}, \code{verbose}, \code{class.stratify.cv}.  Note that other \code{distribution} arguments have not been tested.
#' @return Fitted model. This is a list containing the following elements:
#' 
#' \itemize{
#'   \item models - list of gbm models for each outcome. Functions from the gbm package (e.g. to compute relative influence, print trees, obtain predictions, etc) can be directly applied to each of these models 
#'   \item covex - covariance explained in each pair of outcomes by each predictor. The covariance explained is only unambiguous if predictors are independent, otherwise it is an approximation. If the interaction.depth is larger than 1, the covariance explained is attributed to the predictor with the largest effect.
#'   \item maxiter - maximum number of iterations run (the number of trees fit)
#'   \item best.trees - list of the best number of trees given by minimizing the multivariate MSE error in the test set, by CV, or just the last tree fit. Many of the functions in the package default to using the minimum value of the three. 
#'   \item rel.infl - n x q x n.trees matrix of relative influences
#'   \item w.rel.infl - n x q x n.trees matrix of weighted relative influences (see details)
#'   \item params - arguments to mvtb
#'   \item trainerr - multivariate training error at each tree.
#'   \item testerr  - multivariate test error at each tree (if trainfrac < 1)
#'   \item cverr    - multivariate cv error at each tree (if cv.folds > 1)
#'   \item bestxs - matrix of predictors selected at each tree
#'   \item resid - n x q matrix of residuals after fitting all trees
#'   \item ocv - if save.cv=TRUE, returns the CV models.
#'   \item wm.raw - raw decreases in covariance attributable to a given tree
#'   \item wm.rel - relative decreases in covariance attributable to a given tree
#'   \item s - indices of training sample
#'   \item n - number of observations
#'   \item xnames
#'   \item ynames
#' }
#' 
#' @details 
#' 
#' This function computes separate gbm models for each outcome (contained in \code{$models}). From these models, the covariance explained by 
#' each predictor is then computed based on the reduction in covariance between outcomes that results from fitting a single tree to each outcome, one outcome at a time.
#' The reduction in covariance between each pair of outcomes due to splitting on each predictor over all trees is the 'covariance explained' by each predictor, and is recorded in \code{$covex}.
#' 
#' The rows (pairs of outcomes) and the columns (predictors) of \code{$covex} can be clustered so that groups of predictors that explain similar pairs of covariances are closer together (see  \code{mvtb.cluster}). 
#' A simple heatmap of this matrix can be obtained by the function \code{mvtb.heat}. The \code{covex} by each predictor is only unambiguous if the predictors are uncorrelated and interaction.depth = 1. 
#' If predictors are not independent, the decomposition of covariance explained is only approximate (like the decomposition of R^2 by each predictor in a linear model). 
#' If interaction.depth > 1, the following heuristic is used: the covariance explained by the tree is assigned to the predictor with the largest influence in each tree.
#' 
#' (Relative) influences can be retrieved using \code{mvtb.ri}, which are the usual reductions in SSE due to splitting on each predictor. 'Weighted' influences for each predictor are the usual reductions in SSE weighted by
#' the covariance explained in all pairs of outcomes by that predictor. This allows predictor selection to be informed by the covariance explained. Higher weight can be given to explaining variances or covariances by controlling the parameter
#' \code{alpha}:  weight = \code{alpha}(varex) + (1-\code{alpha})(covex). By default, \code{alpha} = .5, equally weighting variance and covariance explained. 
#' Setting \code{alpha} = 0 corresponds to weighting covariance explained only,
#' letting \code{alpha} = 1 corresponds to weighting variance explained only.
#' 
#' The model is tuned jointly by selecting the number of trees that minimize multivariate mean squared error in a test set (by setting \code{trainfrac}) or averaged over k folds in k-fold cross-validation (by setting \code{cv.folds > 1}).
#' The best number of trees is avaiable via \code{$best.trees}. Cross-validation can be parallelized by setting mc.cores > 1.  
#' If both \code{cv.folds} and \code{trainfrac} is specified, cross-validation is carried out within the training set.
#' Cross-validation models are usually discarded but can be saved by setting \code{save.cv = TRUE}. CV models can be accessed from \code{$ocv} of the 
#' output object. 
#' Multivariate mean squared training, test, and cv error are available from \code{$trainerr, $testerr, $cverr} from the output object.
#' Observations can be specifically set for inclusion in the training set by passing a vector of integers indexing the rows to include to \code{s}.
#' If \code{s} is specified, \code{trainfrac} is ignored but cross-validation will be carried out for observations in \code{s}.
#' 
#' Since the output objects can be large, automatic compression is available by setting \code{compress=TRUE}. 
#' All methods that use the \code{mvtb} object automatically uncompress this object if necessary. 
#' The function \code{mvtb.uncomp} is available to manually decompress the object.
#' 
#' Note that trees are grown until a minimum number of observations in each node is reached. 
#' If the number of training samples*bag.fraction is less the minimum number of observations, (which can occur with small data sets), this will throw an error. 
#' Adjust the \code{n.minobsinnode}, \code{tranfrac}, or \code{bag.fraction}.
#' 
#' This is a beta version with details subject to change. Any contributions are welcome.
#' @seealso \code{summary.mvtb}, \code{predict.mvtb}, \code{mvtb.nonlin} to help detect nonlinear effects, \code{plot.mvtb}, \code{mvtb.perspec} for plots, \code{mvtb.cluster}, \code{mvtb.heat} 
#' \code{mvtb.uncomp} to uncompress a compressed output object
#' @references Miller P.J., Lubke G.H, McArtor D.B., Bergeman C.S. (Submitted) Finding structure in data: A data mining alternative to multivariate multiple regression. Psychological Methods.
#' 
#' Ridgeway, G., Southworth, M. H., & RUnit, S. (2013). Package 'gbm'. Viitattu, 10, 2013.
#'
#' Elith, J., Leathwick, J. R., & Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology, 77(4), 802-813.
#'  
#' Friedman, J. H. (2001). Greedy function approximation: a gradient boosting machine. Annals of statistics, 1189-1232.
#' @examples
#' set.seed(123)
#' n <- 1000
#' X <- matrix(rbinom(n*3,size=1,prob=.5),n,3)    # create 3 dichotomous predictors
#' X2 <- cbind(x1x2=X[,1]*X[,2],x2x3=X[,2]*X[,3]) # create 2 interaction terms
#' Xf <- cbind(X,X2)                              # full design matrix, used for data generation
#' E <- matrix(rnorm(n*4),nrow=n,ncol=4)          # independent errors
#' B <- matrix(0,nrow=5,ncol=4)
#' B[4,1] <- 1      # x1x2 interaction has a true effect on outcome 1
#' B[5,3:4] <- 1    # x2x3 interaction has a true effect on outcomes 3 and 4
#' Y <- Xf %*% B + E
#' 
#' 
#' out <- mvtb(
#'  X=X,                   # matrix of predictors
#'  Y=Y,                   # matrix of responses
#'  n.trees=100,           # number of trees
#'  shrinkage=.1,          # shrinkage or learning rate
#'  interaction.depth = 5, # number of splits in each tree
#'  bag.frac = .5          # bagging fraction
#'  )
#' 
#' summary(out)
#' plot(out)
#' mvtb.nonlin(out,X=X,Y=Y)
#' mvtb.perspec(out)
#' mvtb.cluster(out)
#' mvtb.heat(out)

#' @export
#' @importFrom stats cov
mvtb <- function(X,Y,n.trees=100,shrinkage=.01,interaction.depth=1,
                 trainfrac=1,bag.frac=1,cv.folds=1,
                 s=NULL,seednum=NULL,compress=FALSE,save.cv=FALSE,mc.cores=1,samp.iter=FALSE,alpha=.5,cov.discrep=1,weight.type=1,...) {

  Y <- as.matrix(Y)
  params <- c(as.list(environment()))
  ## seeds
  if(!is.null(seednum)){
    #print(c("mvtboost: seednum ",seednum));
    set.seed(seednum)
  }
  stopifnot(nrow(X) == nrow(Y))
  
  ## Data
  n <- nrow(X); k <- ncol(Y); p <- ncol(X);
  if(is.null(colnames(X))) { colnames(X) <- paste("X",1:p,sep="") }
  if(is.null(colnames(Y))) { colnames(Y) <- paste("Y",1:k,sep="") } 
  xnames <- colnames(X)
  ynames <- colnames(Y)
  D <- Rm <- matrix(0,n,k)            
  
  ## sampling
  if(is.null(s)){
    params$s <- s <- sample(1:n,floor(n*trainfrac),replace=F) #force round down if odd
  } 
  ## ss <- matrix(0,nrow=length(s),ncol=n.trees)

  ## 1. bootstrap covex if desired.
  ## for(i in 1:n.trees) {
  ##  if(samp.iter) {
  ##    ss[,i] <-  sample(s,length(s),replace=TRUE) # if replace = FALSE, this just permutes the rows.
  ##  } else {
  ##    ss[,i] <- s
  ##  }
  #}  
  
  ## parameters
  plist <- params
  plist$cov.discrep <- NULL
  plist$alpha       <- NULL
  plist$weight.type <- NULL
  
  ## Checks
  if(any(is.na(Y))){ stop("NAs not allowed in response variables.")}
    
  ## Influence
  wm.raw <- wm.rel <- matrix(0,nrow=n.trees,ncol=k)     #raw, relative
  rel.infl <- w.rel.infl <- array(0,dim=c(p,k,n.trees)) # influences at every iteration
  
  ## Iterations
  trainerr <- testerr <- vector(length=n.trees)
  bestxs   <- matrix(0,nrow=n.trees,ncol=k) 
  final.iter <- FALSE
  
  ## Covex
  Res.cov <- array(0,dim=c(k,k,k,n.trees))
  covex <- array(0,dim=c(p,k*(k+1)/2))
  rownames(covex) <- colnames(X)
  names <- outer(1:k,1:k,function(x,y){paste0(colnames(Y)[x],"-",colnames(Y)[y])})
  colnames(covex) <- names[lower.tri(names,diag=TRUE)]
  
  ## 0. CV?
  if(cv.folds > 1) {
    ocv <- mvtbCV(params=plist)
    best.iters.cv <- ocv$best.iters.cv
    cv.err <- ocv$cv.err
    out.fit <- ocv$models.k[[cv.folds+1]]
   } else {
    plist$mc.cores <- NULL
    out.fit <- do.call("mvtb.fit",args=c(plist))
    best.iters.cv <- NULL
    cv.err <- NULL
    ocv <- NULL
  }
  models <-  out.fit$models
  trainerr <- out.fit$trainerr
  testerr <- out.fit$testerr
  yhat <- out.fit$yhat
  
  init <- unlist(lapply(models,function(m){m$initF}))
  for(m in 1:k) { D[,m] <- Y[,m]-init[m] } # current residuals at iteration 1
  yhat <- array(c(rep(init,each=n),yhat),dim=c(n,k,n.trees+1))
  Dpred <- yhat[,,2:(n.trees+1),drop=F]-yhat[,,1:(n.trees),drop=F] # Dpred is the unique contribution of each tree
  
  # 1. Get trees from each model
  finaltree <- list()
  for(m in 1:k) { 
    finaltree[[m]] <- models[[m]]$trees
  }
  # 2. Compute covariance discrepancy
  for(i in 1:(n.trees)) {        
    
    ## 2.1 From each model get the stuff we need at the current iteration
    for(m in 1:k) {            
      ## 1.2 For each model compute predicted values and influence
      tree.i <- finaltree[[m]][i]
      rel.infl[,m,i] <- ri.one(tree.i,n.trees=1,xnames)
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
    
    ## 7. Compute best iteration criteria: training and test error, sum of residual covariace matrix, weighted sum of residual covariance matrix
    #trainerr[i] <- mean((as.matrix(Rm[s,]))^2,na.rm=TRUE)
    #if(length(unique(s)) == nrow(P)) { # use training error instead of testing
    #  testerr[i] <- trainerr[i]            
    #} else {
    #  testerr[i] <- mean((as.matrix(Rm[-s,]))^2,na.rm=TRUE)
    #}
    
    ## 8. Compute weighted influences
    for(m in 1:k) {
      #rel.infl[,m] <- relative.influence(models[[m]],n.trees=i,scale=FALSE)
      w.rel.infl[,m,i] <- rel.infl[,m,i]*wm.rel[i,m]
    }
    
    if(final.iter) {
      i <- i-1 # don't count the last iteration if all weights are  0
      break;
    }
  }
 
  
  best.trees <- list(best.testerr=which.min(testerr),best.cv=best.iters.cv,last=i)

  covex <- t(covex)
  if(!save.cv){ocv <- NULL}

  fl <- list(models=models, covex=covex,maxiter=i,best.trees=best.trees,
             rel.infl=rel.infl, w.rel.infl=w.rel.infl,params=params,
             trainerr=trainerr[1:i],testerr=testerr[1:i],cv.err=cv.err[1:i],
             bestxs=bestxs,
             resid=D,ocv=ocv,
             wm.raw=matrix(wm.raw[1:i,,drop=FALSE],nrow=i,ncol=k),wm.rel=wm.rel[1:i,,drop=FALSE],
             s=s,n=nrow(X),xnames=colnames(X),ynames=colnames(Y))
  if(compress) {
    # compress each element using bzip2
    fl <- lapply(fl,comp)
  }
  class(fl) <- "mvtb"
  return(fl)
}


mvtb.fit <- function(X,Y,n.trees=100,shrinkage=.01,interaction.depth=1,
                           trainfrac=1,samp.iter=FALSE,bag.frac=1,cv.folds=1,
                           s=NULL,seednum=NULL,compress=FALSE,save.cv=FALSE,...) {
    if(!is.null(seednum)){
      #print(c("mvtboost: seednum ",seednum));
      set.seed(seednum)
    }

    ## 0. Generate random training samples. This is done before the fitting for reproducibility.
    n <- nrow(X) 
    q <- ncol(Y)

    ## 2. Fit the models
    models <- pred <- list()
    for(m in 1:q) {
      models[[m]] <- gbm::gbm.fit(x=as.data.frame(X[s,,drop=FALSE]),y=matrix(Y[s,m,drop=FALSE]),
                              shrinkage=shrinkage,interaction.depth=interaction.depth,
                              n.trees=n.trees,verbose=F,distribution="gaussian",
                              bag.fraction=bag.frac,keep.data=FALSE,...)
    }
    ## Multivariate mse. clumsy, I know.
    yhat <- predict.mvtb(list(models=models),n.trees=1:n.trees,newdata=X,drop=FALSE)
    testerr <- trainerr <- rep(0,length=n.trees)
    for(i in 1:n.trees){
      R <- Y-yhat[,,i]
      trainerr[i] <- mean((as.matrix(R[s,]))^2,na.rm=TRUE)
      testerr[i] <- mean((as.matrix(R[-s,]))^2,na.rm=TRUE)
    }
    
    ## 3. Compute multivariate MSE
    fl <- list(models=models,trainerr=trainerr,testerr=testerr,yhat=yhat,s=s)
    return(fl)
}

ri.one <- function(object,n.trees=1,var.names) {
  get.rel.inf <- function(obj) {
    lapply(split(obj[[6]], obj[[1]]), sum)
  }
  temp <- unlist(lapply(object[1:n.trees], get.rel.inf))
  rel.inf.compact <- unlist(lapply(split(temp, names(temp)), 
                                   sum))
  rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) != 
                                       "-1"]
  rel.inf <- rep(0, length(var.names))
  i <- as.numeric(names(rel.inf.compact)) + 1
  rel.inf[i] <- rel.inf.compact
  return(rel.inf = rel.inf)
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

## generalized least squares. requires precomputed inverse of the sample covariance matrix of the responses.
#gls <- function(Sd,Sinv) { .5 * sum(diag( (Sd %*% Sinv) %*% t(Sd %*% Sinv) ),na.rm=TRUE) }

comp <- function(obj) { memCompress(serialize(obj,NULL),"bzip2") }
uncomp <-function(obj){ unserialize(memDecompress(obj,type="bzip2"))}

## the only tricky thing here is that some models can stop at earlier iterations than others.
## The default error for each is NA, and the average error (rowMeans over k) is computed with na.rm=TRUE.
## Thus the error estimate at later iterations may not be based on all folds.
mvtbCV <- function(params) {
    n <- nrow(params$X)
    cv.folds <- params$cv.folds
    save.cv <- params$save.cv
    
    #if(is.null(params$s)) {
    #    s <- sample(1:n,floor(n*params$trainfrac))
    #}
    sorig <- params$s  # this is the original shuffling/subsetting from mvtb
    cv.groups <- sample(rep(1:cv.folds,length=length(sorig)))

    # construct the new call
    params$trainfrac <- 1
    params$cv.folds <- 1
    params$save.cv <- save.cv
    #params$X <- params$X[sorig,, drop=FALSE] # training fraction
    #params$Y <- params$Y[sorig,, drop=FALSE]

    testerr.k <- matrix(NA,nrow=params$n.trees,ncol=cv.folds)
    out.k <- list()

    runone <- function(k,params,cv.groups,sorig){
      if(any(k %in% cv.groups)) {
        params$s <- sorig[which(cv.groups != k)]
      } else { 
        # since we already subsetted on s, fit to entire training sample
        params$s <- sorig
        params$compress <- FALSE        
      }
      out <- do.call("mvtb.fit",params) 
      return(out)
    }
    # Last fold contains the full sample
    if(params$mc.cores > 1) {
        cores <- params$mc.cores
        #cluster <- parallel::makeCluster(cores)
        params$mc.cores <- NULL
        out.k <- parallel::mclapply(1:(cv.folds+1),runone,params=params,cv.groups=cv.groups,sorig=sorig,mc.cores=cores)
    } else {
      params$mc.cores <- NULL
      out.k <- lapply(1:(cv.folds+1),runone,params=params,cv.groups=cv.groups,sorig=sorig)
    }
        
    for(k in 1:cv.folds) {
      out <- out.k[[k]]
      testerr.k[,k] <- out$testerr
    }
  
    cv.err <- rowMeans(testerr.k,na.rm=TRUE)
    best.iters.cv <- which.min(cv.err)
    
 
    if(params$save.cv) {
        l <- list(models.k=out.k,best.iters.cv=best.iters.cv,cv.err=cv.err,cv.groups=cv.groups)
    } else {
        out.k[1:cv.folds] <- NULL
        models.k <- c(lapply(1:cv.folds,list),out.k)
        l <- list(models.k=models.k,best.iters.cv=best.iters.cv,cv.err=cv.err)
    }
    return(l)
}

### Examples ###

#assign default values of arguments to workspace
## list2env(formals(mvtb)[-c(1:2)],globalenv())
#  params <- formals(mvtb)[-c(length(formals(mvtb)))]
# 
# B <- matrix(c(.5,.5,0,0,.5,.5),ncol=3,byrow=TRUE)
# B
# varx <- 2
# vare <- diag(3)
# #diag(vare) <- 1:3
# vare[lower.tri(vare)] <- vare[upper.tri(vare)] <- 0
# 
# n <- 1000
# X <- matrix(rnorm(n*2,1,2),ncol=2)
# E <- mvrnorm(n,mu=c(0,0,0),Sigma=vare)
# Y <- X %*% B + E
# bhat <- solve(t(X) %*% X) %*% t(X) %*% Y
# cov(Y)
# cov(X %*% B)
# colnames(X) <- paste0("X",1:varx)
# params$X <- X
# params$Y <- Y

#out <- mvtb(X=X,Y=Y,shrinkage=1,n.trees=200,weight.type=2,weighted.predictions=TRUE)
#out$covex
#out$rel.covex

## I think it works!

#' Predicted values
#' @param object mvtb object
#' @param newdata matrix of predictors. 
#' @param n.trees number of trees. If a vector, returns predictions in an array. Defaults to the minimum best tree estimate.
#' @param drop Drop any dimensions of length 1.
#' @param ... unused
#' @return Returns an (array, matrix, vector) of predictions for all outcomes. The third dimension corresponds to the 
#' predictions at a given number of trees, the second dimension corresponds to the number of outcomes.
#' @export
#' 
predict.mvtb <- function(object, n.trees=NULL, newdata, drop=TRUE,...) {
  out <- object
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- mvtb.uncomp(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  K <- length(out$models)
  treedim <- ifelse(length(n.trees) > 1,max(n.trees),1)
  Pred <- array(0,dim=c(nrow(newdata),K,treedim))  
  for(k in 1:K) {                                     
    p <- rep(0,nrow(newdata))        
    p <- gbm::predict.gbm(out$models[[k]],n.trees=n.trees,newdata=newdata)    
    Pred[,k,] <- p
  }
  #if(length(n.trees) == 1) {
  #  Pred <- drop(Pred)
  #}
  if(drop){
    Pred <- drop(Pred)
  }
  return(Pred)
}

#residuals.mvtb <- fucntion(object,n.trees,...){
#  
#}


