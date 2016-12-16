#' Fitting a Multivariate Tree Boosting Model
#'
#' Builds on \code{gbm} (Ridgeway 2013; Friedman, 2001) to fit a univariate tree model for each outcome, selecting predictors at each iteration that explain (co)variance in the outcomes. The number of trees included in the model can be chosen by minimizing the multivariate mean squared error using cross validation or a test set.
#'
#' @param Y vector, matrix, or data.frame for outcome variables with no missing values. To easily compare influences across outcomes and for numerical stability, outcome variables should be scaled to have unit variance.
#' @param X vector, matrix, or data.frame of predictors. For best performance, continuous predictors should be scaled to have unit variance. Categorical variables should converted to factors.
#' @param n.trees maximum number of trees to be included in the model. Each individual tree is grown until a minimum number observations in each node is reached. 
#' @param shrinkage a constant multiplier for the predictions from each tree to ensure a slow learning rate. Default is .01. Small shrinkage values may require a large number of trees to provide adequate fit.
#' @param interaction.depth fixed depth of trees to be included in the model. A tree depth of 1 corresponds to fitting stumps (main effects only), higher tree depths capture higher order interactions (e.g. 2 implies a model with up to 2-way interactions)
#' @param distribution Character vector specifying the distribution of all outcomes. Default is "gaussian" see ?gbm for further details.
#' @param train.fraction  proportion of the sample used for training the multivariate additive model. If both \code{cv.folds} and \code{train.fraction} are specified, the CV is carried out within the training set.
#' @param bag.fraction   proportion of the training sample used to fit univariate trees for each response at each iteration. Default: 1
#' @param cv.folds   number of cross validation folds. Default: 1. Runs k + 1 models, where the k models are run in parallel and the final model is run on the entire sample. If larger than 1, the number of trees that minimize the multivariate MSE averaged over k-folds is reported in \code{object$best.trees}
#' @param keep.data a logical variable indicating whether to keep the data stored with the object.
#' @param s vector of indices denoting observations to be used for the training sample. If \code{s} is given, \code{train.fraction} is ignored.
#' @param seednum integer passed to \code{set.seed}
#' @param compress \code{TRUE/FALSE}. Compress output results list using bzip2 (approx 10\% of original size). Default is \code{FALSE}.
#' @param save.cv  \code{TRUE/FALSE}. Save all k-fold cross-validation models. Default is \code{FALSE}.
#' @param iter.details \code{TRUE/FALSE}. Return training, test, and cross-validation error at each iteration. Default is \code{FALSE}.
#' @param verbose If \code{TRUE}, will print out progress and performance indicators for each model.  Default is \code{FALSE}.
#' @param mc.cores Number of cores for cross validation.
#' @param ... additional arguments passed to \code{gbm}. These include \code{distribution}, \code{weights}, \code{var.monotone}, \code{n.minobsinnode}, \code{keep.data}, \code{verbose}, \code{class.stratify.cv}.  Note that other \code{distribution} arguments have not been tested.
#' @return Fitted model. This is a list containing the following elements:
#' 
#' \itemize{
#'   \item \code{models} - list of gbm models for each outcome. Functions from the gbm package (e.g. to compute relative influence, print trees, obtain predictions, etc) can be directly applied to each of these models 
#'   \item \code{best.trees} - A list containing  the number of trees that minimize the multivariate MSE in a test set or by CV, and \code{n.trees}.
#'     Many of the functions in the package default to using the minimum value of the three. 
#'   \item \code{params} - arguments to mvtb
#'   \item \code{train.err} - multivariate training error at each tree (If \code{iter.details = TRUE})
#'   \item \code{test.err}  - multivariate test error at each tree (if \code{train.fraction < 1} and \code{iter.details = TRUE})
#'   \item \code{cv.err}    - multivariate cv error at each tree (if \code{cv.folds > 1} and \code{iter.details = TRUE})
#'   \item \code{cv.mods} - the CV models if \code{save.cv=TRUE}
#'   \item \code{s} - indices of training sample
#'   \item \code{n} - number of observations
#'   \item \code{xnames}
#'   \item \code{ynames}
#' }
#' 
#' @usage 
#' mvtb(Y, X, 
#'      n.trees = 100,
#'      shrinkage = 0.01, 
#'      interaction.depth = 1,
#'      distribution="gaussian",
#'      train.fraction = 1, 
#'      bag.fraction = 1, 
#'      cv.folds = 1, 
#'      keep.data = FALSE,
#'      s = NULL, 
#'      seednum = NULL, 
#'      compress = FALSE, 
#'      save.cv = FALSE,
#'      iter.details = TRUE,
#'      verbose=FALSE,
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
#' The model is tuned jointly by selecting the number of trees that minimize multivariate mean squared error in a test set (by setting \code{train.fraction}) or averaged over k folds in k-fold cross-validation (by setting \code{cv.folds > 1}).
#' The best number of trees is available via \code{$best.trees}.  
#' If both \code{cv.folds} and \code{train.fraction} is specified, cross-validation is carried out within the training set.
#' If \code{s} is specified, \code{train.fraction} is ignored but cross-validation will be carried out for observations in \code{s}.
#' 
#' Cross-validation models are usually discarded but can be saved by setting \code{save.cv = TRUE}. CV models can be accessed from \code{$ocv} of the 
#' output object. Observations can be specifically set for inclusion in the training set by passing a vector of integers indexing the rows to include to \code{s}.
#' Multivariate mean squared training, test, and cv error are available from \code{$train.error, $test.error, $cverr} from the output object 
#' when \code{iter.details = TRUE}.
#' 
#' Since the output objects can be large, automatic compression is available by setting \code{compress=TRUE}. 
#' All methods that use the \code{mvtb} object automatically uncompress this object if necessary. 
#' The function \code{mvtb.uncomp} is available to manually decompress the object.
#' 
#' Note that trees are grown until a minimum number of observations in each node is reached. 
#' If the number of \code{training samples}*\code{bag.fraction} is less the minimum number of observations, (which can occur with small data sets), this will cause an error. 
#' Adjust the \code{n.minobsinnode}, \code{train.fraction}, or \code{bag.fraction}.
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
#' covex <- mvtb.covex(res, Y=Ys, X=Xs)
#' plot(res,predictor.no = 8)
#' predict(res,newdata=Xs)
#' mvtb.cluster(covex)
#' mvtb.heat(t(mvtb.ri(res)),cexRow=.8,cexCol=1,dec=0)
#' @export
mvtb <- function(Y,X,n.trees=100,
                 shrinkage=.01,
                 interaction.depth=1,
                 distribution="gaussian",
                 train.fraction=1,
                 bag.fraction=1,
                 cv.folds=1,
                 keep.data=FALSE,
                 s=NULL,
                 seednum=NULL,
                 compress=FALSE,
                 save.cv=FALSE,
                 iter.details=TRUE,
                 verbose=FALSE,
                 mc.cores=1,...) {

  n <- nrow(X)
  if(class(Y) != "matrix") { Y <- as.matrix(Y) }
  if(is.null(ncol(X))){ X <- as.matrix(X)}

  params <- c(as.list(environment()),list(...)) # this won't copy y and x
  
  if(!is.null(seednum)) set.seed(seednum)

  stopifnot(nrow(X) == nrow(Y))
  
  n <- nrow(X); k <- ncol(Y); p <- ncol(X);
  if(is.null(colnames(X))) { colnames(X) <- paste("X",1:p,sep="") }
  if(is.null(colnames(Y))) { colnames(Y) <- paste("Y",1:k,sep="") } 
  xnames <- colnames(X)
  ynames <- colnames(Y)
  
  ## sampling
  if(is.null(s)){
    s <- sample(1:n, floor(n*train.fraction), replace=F) #force round down if odd
  } 
  
  ## Checks
  if(shrinkage > 1 | shrinkage <= 0){ stop("shrinkage should be > 0, < 1")}
  if(train.fraction > 1 | train.fraction <= 0){ stop("train.fraction should be > 0, < 1")}
  if(bag.fraction > 1 | bag.fraction <= 0){ stop("bag.fraction should be > 0, < 1")}
  
  ## Iterations
  train.error <- test.error <- vector(length=n.trees)
  
  ## 0. CV?
  if(cv.folds > 1) {
    cv.mods <- mvtbCV(Y=Y, X=X, cv.folds=cv.folds, s=s, save.cv=save.cv, mc.cores=mc.cores,
                  n.trees=n.trees, shrinkage=shrinkage, interaction.depth=interaction.depth, 
                  distribution=distribution, bag.fraction=bag.fraction, verbose=verbose,
                  keep.data=keep.data, seednum=seednum)
    best.iters.cv <- cv.mods$best.iters.cv
    cv.err <- cv.mods$cv.err
    out.fit <- cv.mods$models.k[[cv.folds+1]]
   } else {
    out.fit <- mvtb.fit(Y=Y, X=X,
                        n.trees=n.trees, shrinkage=shrinkage, interaction.depth=interaction.depth,
                        distribution=distribution, bag.fraction=bag.fraction, verbose=verbose,
                        s=s, keep.data=keep.data, seednum=seednum,...)
    best.iters.cv <- NULL
    cv.err <- NULL
    cv.mods <- NULL
  }
  models <-  out.fit$models
  train.err <- out.fit$train.error
  test.err <- out.fit$test.error
  cv <- ifelse(is.null(best.iters.cv), NA, best.iters.cv)
  test <- ifelse(all(is.nan(test.error)), NA, which.min(test.err))

  # can't do oob
  best.trees <- list(train=which.min(train.err), test=test, oob=NA, cv=cv)
  best.trees <- do.call(rbind, lapply(1:k, function(i){data.frame(best.trees)}))
  rownames(best.trees) <- colnames(Y)

  if(!save.cv){cv.mods <- NULL}
  if(!iter.details){train.err <- NULL; test.err <- NULL; cv.err = NULL}
  
  fl <- list(models=models, best.trees=best.trees, params=params,
             train.err=train.err, test.err=test.err, cv.err=cv.err,
             cv.mods=cv.mods,
             s=s,n=nrow(X), xnames=colnames(X), ynames=colnames(Y))
  
  if(compress) {
    # compress each element using bzip2
    fl <- lapply(fl,comp)
  }
  class(fl) <- "mvtb"
  return(fl)
}

#' @describeIn mvtb 
#' @usage 
#' mvtb.fit(Y,X,
#'          n.trees=100,
#'          shrinkage=.01,
#'          interaction.depth=1,
#'          bag.fraction=1,
#'          s=1:nrow(X),
#'          seednum=NULL,...)
#' @export
mvtb.fit <- function(Y,X,
                     n.trees=100,
                     shrinkage=.01,
                     interaction.depth=1,
                     bag.fraction=1,
                     s=1:nrow(X),
                     seednum=NULL,...) {
    if(!is.null(seednum)){
      #print(c("mvtboost: seednum ",seednum));
      set.seed(seednum)
    }
    
    n <- nrow(X) 
    q <- ncol(Y)

    ## 2. Fit the models
    models <- pred <- list()
    for(m in 1:q) {
      models[[m]] <- gbm::gbm.fit(y=matrix(Y[s,m,drop=FALSE]), 
                                  x=as.data.frame(X[s,,drop=FALSE]),
                                  n.trees=n.trees, 
                                  interaction.depth=interaction.depth, 
                                  shrinkage=shrinkage,
                                  bag.fraction=bag.fraction,
                                  ...)
    }
    yhat <- predict.mvtb.array(list(models=models),n.trees=1:n.trees,newdata=X,drop=FALSE)
    test.error <- train.error <- rep(0,length=n.trees)
    for(i in 1:n.trees){
      R <- Y-yhat[,,i]
      train.error[i] <- mean((as.matrix(R[s,]))^2,na.rm=TRUE)
      test.error[i] <- mean((as.matrix(R[-s,]))^2,na.rm=TRUE)
    }
    
    ## 3. Compute multivariate MSE
    fl <- list(models=models, train.error=train.error, test.error=test.error, s=s)
    return(fl)
}


## generalized least squares. requires precomputed inverse of the sample covariance matrix of the responses.
#gls <- function(Sd,Sinv) { .5 * sum(diag( (Sd %*% Sinv) %*% t(Sd %*% Sinv) ),na.rm=TRUE) }

comp <- function(obj) { memCompress(serialize(obj,NULL),"bzip2") }
uncomp <-function(obj){ unserialize(memDecompress(obj,type="bzip2"))}

## the only tricky thing here is that some models can stop at earlier iterations than others.
## The default error for each is NA, and the average error (rowMeans over k) is computed with na.rm=TRUE.
## Thus the error estimate at later iterations may not be based on all folds.
mvtbCV <- function(Y, X, n.trees, cv.folds, save.cv, s, mc.cores, ...) {
    n <- nrow(X)
    
    cv.groups <- sample(rep(1:cv.folds,length=length(s)))
    
    test.error.k <- matrix(NA,nrow=n.trees,ncol=cv.folds)
    out.k <- list()

    runone <- function(k,cv.groups,sorig, x, ...){
      if(any(k %in% cv.groups)) {
        s <- sorig[which(cv.groups != k)]
      } else { 
        # since we already subsetted on s, fit to entire training sample
        s <- sorig
      }
      out <- mvtb.fit(s=s, X=x, ...) # the only thing that changes is s. everything is passed via ..., including Y and X
      return(out)
    }
    
    # Last fold contains the full sample
    ## The 'if' notation is just to make sure it works on windows.
    if(.Platform$OS.type == "unix") { 
      out.k <- parallel::mclapply(1:(cv.folds + 1), FUN=runone, cv.groups=cv.groups, sorig=s, Y=Y, x=X, n.trees=n.trees, ..., mc.cores=mc.cores)
    } else {
      out.k <- parallel::mclapply(1:(cv.folds + 1), FUN=runone, cv.groups=cv.groups, sorig=s, Y=Y, x=X, n.trees=n.trees, ..., mc.cores=1)
    }
        
    for(k in 1:cv.folds) {
      test.error.k[,k] <- out.k[[k]]$test.error
    }
  
    cv.err <- rowMeans(test.error.k,na.rm=TRUE)
    best.iters.cv <- which.min(cv.err)
    
    if(save.cv) {
        l <- list(models.k=out.k,best.iters.cv=best.iters.cv,cv.err=cv.err,cv.groups=cv.groups)
    } else {
        out.k[1:cv.folds] <- NULL
        models.k <- c(lapply(1:cv.folds,list),out.k)
        l <- list(models.k=models.k,best.iters.cv=best.iters.cv,cv.err=cv.err)
    }
    return(l)
}

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

#' Predicted values
#' @param object \code{mvtb} object
#' @param newdata matrix of predictors. 
#' @param n.trees (vector) of the number of trees for each outcome. Defaults to the minimum number of trees by CV, test, or training error for each outcome.
#' @param ... not used
#' @return Returns a matrix or vector of predictions for all outcomes. 
#' @export
predict.mvtb <- function(object, newdata, n.trees=NULL, ...) {
  if(any(unlist(lapply(object,function(li){is.raw(li)})))){
    object <- mvtb.uncomp(object)
  }
  K <- length(object$models)
  if(is.null(n.trees)) { n.trees <- apply(object$best.trees, 1, min, na.rm=T) }
  if(length(n.trees) == 1){ n.trees <- rep(n.trees, K)}
  
  Pred <- matrix(0,nrow(newdata),K)
  for(k in 1:K) {                                     
    p <- rep(0,nrow(newdata))        
    p <- predict(object$models[[k]],n.trees=n.trees[k],newdata=data.frame(newdata))
    Pred[,k] <- p
  }
  return(drop(Pred))
}



