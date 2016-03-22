## Author: Patrick Miller
## Purpose: Multiple response (or multivariate) boosting with decision trees (stumps). Also multi-task, or multi-objective. Only for continuous Y.
## RRV: 9/3/2015
  

#' Fitting a Multivariate Tree Boosting Model
#'
#' Builds on \code{gbm} (Ridgeway 2013; Friedman, 2001) to fit a univariate tree model for each outcome, selecting predictors at each iteration that explain (co)variance in the outcomes. The number of trees included in the model can be chosen by minimizing the multivariate mean squared error using cross validation or a test set.
#'
#' @param X vector, matrix, or data.frame of predictors. For best performance, continuous predictors should be scaled to have unit variance. Categorical variables should converted to factors.
#' @param Y vector, matrix, or data.frame for outcome variables with no missing values. To easily compare influences across outcomes and for numerical stability, outcome variables should be scaled to have unit variance.
#' @param n.trees maximum number of trees to be included in the model. Each individual tree is grown until a minimum number observations in each node is reached. 
#' @param shrinkage a constant multiplier for the predictions from each tree to ensure a slow learning rate. Default is .01. Small shrinkage values may require a large number of trees to provide adequate fit.
#' @param interaction.depth fixed depth of trees to be included in the model. A tree depth of 1 corresponds to fitting stumps (main effects only), higher tree depths capture higher order interactions (e.g. 2 implies a model with up to 2-way interactions)
#' @param bag.frac   proportion of the training sample used to fit univariate trees for each response at each iteration. Default: 1
#' @param cv.folds   number of cross validation folds. Default: 1. Runs k + 1 models, where the k models are run in parallel and the final model is run on the entire sample. If larger than 1, the number of trees that minimize the multivariate MSE averaged over k-folds is reported in \code{object$best.trees}
#' @param trainfrac  proportion of the sample used for training the multivariate additive model. If both \code{cv.folds} and \code{trainfrac} are specified, the CV is carried out within the training set.
#' @param s vector of indices denoting observations to be used for the training sample. If \code{s} is given, \code{trainfrac} is ignored.
#' @param seednum integer passed to \code{set.seed}
#' @param compress \code{TRUE/FALSE}. Compress output results list using bzip2 (approx 10\% of original size). Default is \code{FALSE}.
#' @param save.cv  \code{TRUE/FALSE}. Save all k-fold cross-validation models. Default is \code{FALSE}.
#' @param iter.details \code{TRUE/FALSE}. Return training, test, and cross-validation error at each iteration. Default is \code{FALSE}.
#' @param mc.cores Number of cores for cross validation.
#' @param ... additional arguments passed to \code{gbm}. These include \code{distribution}, \code{weights}, \code{var.monotone}, \code{n.minobsinnode}, \code{keep.data}, \code{verbose}, \code{class.stratify.cv}.  Note that other \code{distribution} arguments have not been tested.
#' @return Fitted model. This is a list containing the following elements:
#' 
#' \itemize{
#'   \item \code{models} - list of gbm models for each outcome. Functions from the gbm package (e.g. to compute relative influence, print trees, obtain predictions, etc) can be directly applied to each of these models 
#'   \item \code{best.trees} - A list containing  the number of trees that minimize the multivariate MSE in a test set or by CV, and \code{n.trees}.
#'     Many of the functions in the package default to using the minimum value of the three. 
#'   \item \code{params} - arguments to mvtb
#'   \item \code{trainerr} - multivariate training error at each tree (If \code{iter.details = TRUE})
#'   \item \code{testerr}  - multivariate test error at each tree (if \code{trainfrac < 1} and \code{iter.details = TRUE})
#'   \item \code{cverr}    - multivariate cv error at each tree (if \code{cv.folds > 1} and \code{iter.details = TRUE})
#'   \item \code{ocv} - the CV models if \code{save.cv=TRUE}
#'   \item \code{s} - indices of training sample
#'   \item \code{n} - number of observations
#'   \item \code{xnames}
#'   \item \code{ynames}
#' }
#' 
#' @usage 
#' mvtb(X, Y, 
#'      n.trees = 100,
#'      shrinkage = 0.01, 
#'      interaction.depth = 1,
#'      trainfrac = 1, 
#'      bag.frac = 1, 
#'      cv.folds = 1, 
#'      s = NULL, 
#'      seednum = NULL, 
#'      compress = FALSE, 
#'      save.cv = FALSE,
#'      iter.details = TRUE,
#'      mc.cores = 1, ...)
#'      
#' @details 
#' 
#' This function selects predictors that explain covariance in multivariate outcomes. 
#' This is done efficiently by fitting separate gbm models for each outcome (contained in \code{$models}). 
#' 
#' (Relative) influences can be retrieved using \code{summary} or \code{mvtb.ri}, which are the usual reductions in SSE due to splitting on each predictor.
#' The covariance explained in pairs of outcomes by each predictor can be computed using \code{mvtb.covex}. 
#' Partial dependence plots can be obtained from \code{mvtb.plot}.
#' 
#' The model is tuned jointly by selecting the number of trees that minimize multivariate mean squared error in a test set (by setting \code{trainfrac}) or averaged over k folds in k-fold cross-validation (by setting \code{cv.folds > 1}).
#' The best number of trees is available via \code{$best.trees}.  
#' If both \code{cv.folds} and \code{trainfrac} is specified, cross-validation is carried out within the training set.
#' If \code{s} is specified, \code{trainfrac} is ignored but cross-validation will be carried out for observations in \code{s}.
#' 
#' Cross-validation models are usually discarded but can be saved by setting \code{save.cv = TRUE}. CV models can be accessed from \code{$ocv} of the 
#' output object. Observations can be specifically set for inclusion in the training set by passing a vector of integers indexing the rows to include to \code{s}.
#' Multivariate mean squared training, test, and cv error are available from \code{$trainerr, $testerr, $cverr} from the output object 
#' when \code{iter.details = TRUE}.
#' 
#' Since the output objects can be large, automatic compression is available by setting \code{compress=TRUE}. 
#' All methods that use the \code{mvtb} object automatically uncompress this object if necessary. 
#' The function \code{mvtb.uncomp} is available to manually decompress the object.
#' 
#' Note that trees are grown until a minimum number of observations in each node is reached. 
#' If the number of \code{training samples}*\code{bag.fraction} is less the minimum number of observations, (which can occur with small data sets), this will cause an error. 
#' Adjust the \code{n.minobsinnode}, \code{trainfrac}, or \code{bag.fraction}.
#' 
#' Cross-validation can be parallelized by setting mc.cores > 1. Parallel cross-validation is carried out using \code{parallel::mclapply}, which makes \code{mc.cores} copies of the original environment.
#' For models with many trees (> 100K), memory limits can be reached rapidly. \code{mc.cores} will not work on Windows. 
#' 
#' @seealso \code{summary.mvtb}, \code{predict.mvtb}
#' 
#' \code{mvtb.covex} to estimate the covariance explained in pairs of outcomes by predictors
#' 
#' \code{mvtb.nonlin} to help detect nonlinear effects or interactions
#' 
#' \code{plot.mvtb}, \code{mvtb.perspec} for partial dependence plots
#'
#' \code{mvtb.uncomp} to uncompress a compressed output object
#' @references Miller P.J., Lubke G.H, McArtor D.B., Bergeman C.S. (Accepted) Finding structure in data with multivariate tree boosting.
#' 
#' Ridgeway, G., Southworth, M. H., & RUnit, S. (2013). Package 'gbm'. Viitattu, 10, 2013.
#'
#' Elith, J., Leathwick, J. R., & Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology, 77(4), 802-813.
#'  
#' Friedman, J. H. (2001). Greedy function approximation: a gradient boosting machine. Annals of statistics, 1189-1232.
#' @examples
#' data(wellbeing)
#' Y <- wellbeing[,21:26]
#' X <- wellbeing[,1:20]
#' Ys <- scale(Y)
#' cont.id <- unlist(lapply(X,is.numeric))
#' Xs <- scale(X[,cont.id])
#' 
#' ## Fit the model
#' res <- mvtb(Y=Ys,X=Xs)
#' 
#' ## Interpret the model
#' summary(res)
#' covex <- mvtb.covex(res, Y=Y, X=X)
#' plot(res,predictor.no = 8)
#' predict(res,newdata=Xs)
#' mvtb.cluster(covex)
#' mvtb.heat(t(mvtb.ri(res)),cexRow=.8,cexCol=1,dec=0)
#' @export
mvtb <- function(X,Y,n.trees=100,
                 shrinkage=.01,
                 interaction.depth=1,
                 trainfrac=1,
                 bag.frac=1,
                 cv.folds=1,
                 s=NULL,
                 seednum=NULL,
                 compress=FALSE,
                 save.cv=FALSE,
                 iter.details=TRUE,
                 mc.cores=1,...) {

  if(class(Y) != "matrix") { Y <- as.matrix(Y) }
  if(is.null(ncol(X))){ X <- as.matrix(X)}
  #if(class(X) != "matrix") { X <- data.matrix(X) }
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
  
  ## sampling
  if(is.null(s)){
    params$s <- sample(1:n,floor(n*trainfrac),replace=F) #force round down if odd
  } 

  ## parameters
  plist <- params
  plist$iter.details <- NULL

  
  ## Checks
  if(any(is.na(Y))){ stop("NAs not allowed in outcome variables.")}
  if(shrinkage > 1 | shrinkage <= 0){ stop("shrinkage should be > 0, < 1")}
  if(trainfrac > 1 | trainfrac <= 0){ stop("trainfrac should be > 0, < 1")}
  if(bag.frac > 1 | bag.frac <= 0){ stop("bag.frac should be > 0, < 1")}
  
  ## Iterations
  trainerr <- testerr <- vector(length=n.trees)
  
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

  best.trees <- list(best.testerr=which.min(testerr),best.cv=best.iters.cv,last=n.trees)

  if(!save.cv){ocv <- NULL}
  if(iter.details==T){
    fl <- list(models=models, best.trees=best.trees,params=params,
             trainerr=trainerr,testerr=testerr,cv.err=cv.err,
             ocv=ocv,
             s=params$s,n=nrow(X),xnames=colnames(X),ynames=colnames(Y))
  } else {
    fl <- list(models=models, best.trees=best.trees,params=params,
               s=params$s,ocv=ocv,n=nrow(X),xnames=colnames(X),ynames=colnames(Y))
  }
  if(compress) {
    # compress each element using bzip2
    fl <- lapply(fl,comp)
  }
  class(fl) <- "mvtb"
  return(fl)
}

#' @describeIn mvtb 
#' @usage 
#' mvtb.fit(X,Y,
#'          n.trees=100,
#'          shrinkage=.01,
#'          interaction.depth=1,
#'          trainfrac=1,
#'          samp.iter=FALSE,
#'          bag.frac=1,
#'          s=NULL,
#'          seednum=NULL,
#'          compress=FALSE,...)
#' @export
mvtb.fit <- function(X,Y,
                     n.trees=100,
                     shrinkage=.01,
                     interaction.depth=1,
                     trainfrac=1,
                     samp.iter=FALSE,
                     bag.frac=1,
                     s=NULL,
                     seednum=NULL,
                     compress=FALSE,...) {
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
                              n.trees=n.trees,verbose=F,
                              bag.fraction=bag.frac,keep.data=FALSE,...)
    }
    yhat <- predict.mvtb(list(models=models),n.trees=1:n.trees,newdata=X,drop=FALSE)
    testerr <- trainerr <- rep(0,length=n.trees)
    for(i in 1:n.trees){
      R <- Y-yhat[,,i]
      trainerr[i] <- mean((as.matrix(R[s,]))^2,na.rm=TRUE)
      testerr[i] <- mean((as.matrix(R[-s,]))^2,na.rm=TRUE)
    }
    
    ## 3. Compute multivariate MSE
    fl <- list(models=models,trainerr=trainerr,testerr=testerr,s=s)
    return(fl)
}


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

#' Predicted values
#' @param object \code{mvtb} object
#' @param newdata matrix of predictors. 
#' @param n.trees number of trees. If a vector, returns predictions in an array. Defaults to the minimum number of trees by CV, test, or training error
#' @param drop \code{TRUE/FALSE} Drop any dimensions of length 1
#' @return Returns an (array, matrix, vector) of predictions for all outcomes. The third dimension corresponds to the 
#' predictions at a given number of trees, the second dimension corresponds to the number of outcomes.
#' @export
#' 
predict.mvtb <- function(object, n.trees=NULL, newdata, drop=TRUE) {
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


