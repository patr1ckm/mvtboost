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
#'      compress = FALSE, 
#'      save.cv = FALSE,
#'      iter.details = TRUE,
#'      verbose=FALSE,
#'      mc.cores = 1, ...)
#'      
#' @details 
#' 
#' This function selects predictors that explain covariance in multivariate outcomes. 
#' This is done efficiently by fitting separate gbm models for each outcome 
#' (contained in \code{$models}). 
#' 
#' (Relative) influences can be retrieved using \code{summary} or \code{mvtb.ri}, which are the usual 
#' reductions in SSE due to splitting on each predictor.
#' The covariance explained in pairs of outcomes by each predictor can be computed using 
#' \code{mvtb.covex}. 
#' Partial dependence plots can be obtained from \code{mvtb.plot}.
#' 
#' The model is tuned by selecting the number of trees that minimize the mean squared error in a test set
#' for each outcome (by setting \code{train.fraction}) or averaged over k folds in k-fold 
#' cross-validation (by setting \code{cv.folds > 1}).
#' The best number of trees is available via \code{$best.trees}.  
#' If both \code{cv.folds} and \code{train.fraction} is specified, cross-validation is carried out 
#' within the training set.
#' If \code{s} is specified, \code{train.fraction} is ignored but cross-validation will be carried out
#'  for observations in \code{s}.
#' 
#' Cross-validation models are usually discarded but can be saved by setting \code{save.cv = TRUE}. 
#' CV models can be accessed from \code{$ocv} of the 
#' output object. Observations can be specifically set for inclusion in the training set by passing
#'  a vector of integers indexing the rows to include to \code{s}.
#' Multivariate mean squared training, test, and cv error are available from \code{$train.error, 
#' $test.error, $cverr} from the output object when \code{iter.details = TRUE}.
#' 
#' Since the output objects can be large, automatic compression is available by setting 
#' \code{compress=TRUE}. 
#' All methods that use the \code{mvtb} object automatically uncompress this object if necessary. 
#' The function \code{mvtb.uncomp} is available to manually decompress the object.
#' 
#' Note that trees are grown until a minimum number of observations in each node is reached. 
#' If the number of \code{training samples}*\code{bag.fraction} is less the minimum number of 
#' observations, (which can occur with small data sets), this will cause an error. 
#' Adjust the \code{n.minobsinnode}, \code{train.fraction}, or \code{bag.fraction}.
#' 
#' Cross-validation can be parallelized by setting mc.cores > 1. Parallel cross-validation is 
#' carried out using \code{parallel::mclapply}, which makes \code{mc.cores} copies of the
#' original environment.
#' For models with many trees (> 100K), memory limits can be reached rapidly.
#' \code{mc.cores} will not work on Windows. 
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
#' @inheritParams  mvtb_joint
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
                 compress=FALSE,
                 save.cv=FALSE,
                 iter.details=TRUE,
                 verbose=FALSE, 
                 mc.cores=1, ...) {
  # Gives Y and X different automatic names. Y and X need to be data frames b/c loops over columns.
  if(!is.data.frame(Y)) Y <- data.frame(Y)  
  if(!is.data.frame(X)) X <- as.data.frame(X)
  n <- nrow(X)
  
  params <- c(as.list(environment()),list(...)) # this won't copy y and x
  ## Checks
  if(shrinkage > 1 | shrinkage <= 0){ stop("shrinkage should be > 0, < 1")}
  if(train.fraction > 1 | train.fraction <= 0){ stop("train.fraction should be > 0, < 1")}
  if(bag.fraction > 1 | bag.fraction <= 0){ stop("bag.fraction should be > 0, < 1")}

  
  if(is.null(s)) { 
    s <- sample(1:n, floor(n*train.fraction), replace=F) #force round down if odd
  }

  do.one <- function(y, x, s, ...){ gbm::gbm.fit(y=y[s], x=as.data.frame(x[s,]), ...) }
  
  get.pred.err <- function(o, y, x, n.trees){
    yhat.iter <- as.matrix(predict(o, newdata=data.frame(x), n.trees=1:n.trees))
    apply(yhat.iter, 2, function(yhat, y){
      mean((y - yhat)^2)
    }, y=y)
  }
  get.test.err <- function(oL, Y, x, s, n.trees){
    test.err <- vector(mode = "list", length=ncol(Y))
    names(test.err) <- colnames(Y)
    for(i in 1:ncol(Y)){
      test.err[[i]] <- get.pred.err(o=oL[[i]], y=Y[-s, i], x=x[-s,], n.trees=n.trees)
    }
    return(test.err)
  }
  
  cv.mods <- cv.mod.err <- vector(mode="list", length(ncol(Y)))
  
  if(cv.folds > 1){
    folds <- sample(1:cv.folds, size=length(s), replace=T)
    for(i in 1:cv.folds){
      train <- s[folds != i]
      test <- s[folds == i]
      cv.mods[[i]] <- parallel::mclapply(Y, FUN=do.one, x=X, n.trees=n.trees, shrinkage=shrinkage, 
                             interaction.depth=interaction.depth, distribution=distribution,
                             bag.fraction=bag.fraction, s=train, 
                             keep.data=keep.data, verbose=verbose, mc.cores=mc.cores,...)
      cv.mod.err[[i]] <- get.test.err(cv.mods[[i]], Y=Y, x=X, n.trees=n.trees, s=test)
    }
  } else {
    folds <- list()
  }
  
  models <- parallel::mclapply(Y, FUN=do.one, x=X, n.trees=n.trees, shrinkage=shrinkage, 
                 interaction.depth=interaction.depth, distribution=distribution,
                 bag.fraction=bag.fraction, s=s, 
                 keep.data=keep.data, verbose=verbose, mc.cores=mc.cores,...)
  
  train.err <- lapply(models, function(o){ o$train.error})
  oob.err <- lapply(models, function(o){ -cumsum(o$oobag.improve)})
  if(cv.folds > 1) {
    cv.err <- lapply(1:ncol(Y), function(i) { rowMeans(sapply(cv.mod.err, "[[", i)) })
  } else {
    cv.err <- lapply(1:ncol(Y), function(i){ rep(NaN, n.trees) })
  }
  names(cv.err) <- colnames(Y)
  if(all(1:n %in% s)){
    test.err <- lapply(models, function(m){ m$valid.error})
  } else {
    test.err <- get.test.err(models, Y=Y, x=X, s=s, n.trees=n.trees)
  }
  
  # which.min can return length 0, return NA instead
  which.min.na <- function(x){
    res <- which.min(x)
    if(length(res) == 0){
      res <- NA
    }
    res
  }
  
  best.trees <- data.frame(train = sapply(train.err, which.min.na), 
                     test = sapply(test.err, which.min.na),
                     oob = sapply(oob.err, which.min.na), 
                     cv = sapply(cv.err, which.min.na))
  rownames(best.trees) <- colnames(Y)
  
  if(!save.cv){cv.mods <- NULL; folds <- NULL}
  if(!iter.details){train.err <- NULL; test.err <- NULL; cv.err = NULL}

  fl <- list(models=models, best.trees=best.trees, params=params,
             train.err=train.err, test.err=test.err, cv.err=cv.err,
             cv.mods=cv.mods, folds=folds,
             s=s,n=nrow(X), xnames=colnames(X), ynames=colnames(Y))

  if(compress) {
    # compress each element using bzip2
    fl <- lapply(fl,comp)
  }

  class(fl) <- "mvtb"
  return(fl)
}

comp <- function(obj) { memCompress(serialize(obj,NULL),"bzip2") }
uncomp <-function(obj){ unserialize(memDecompress(obj,type="bzip2"))}


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

