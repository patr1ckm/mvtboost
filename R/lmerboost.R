
#' Boosted decision trees with random effects
#' 
#' At each iteration, a single decision tree is fit using \code{gbm.fit}, and 
#' the terminal node means are allowed to vary by group using \code{lmer}.
#' 
#' Meta-parameter tuning is handled by passing vectors of possible values for 
#' \code{n.trees}, \code{shrinkage}, \code{indep}, \code{interaction.depth}, 
#' and \code{n.minobsinnode} and setting \code{cv.folds > 1}. Setting 
#' \code{mc.cores > 1} will carry out the tuning in parallel by forking via 
#' \code{mclapply}. Tuning is only done within the training set.
#' 
#' Prediction is most easily carried out by passing the entire \code{X} matrix to
#' \code{lmerboost}, and specifying the training set using \code{subset}. Otherwise, 
#' set \code{save.mods=TRUE} and use \code{predict}.
#' 
#' @param y outcome vector (continuous)
#' @param X matrix or data frame of predictors 
#' @param id name or index of grouping variable
#' @param n.trees the total number of trees to fit (iterations). 
#' @param cv.folds number of cross-validation folds. In addition to the usual fit,
#'  will perform cross-validation over a grid of meta-parameters (see details). 
#' @param interaction.depth The maximum depth of trees. 1 implies a single split
#'  (stump), 2 implies a tree with 2 splits, etc.
#' @param n.minobsinnode minimum number of observations in the terminal nodes
#'  of each tree
#' @param shrinkage a shrinkage parameter applied to each tree. Also known as the 
#' learning rate or step-size reduction.
#' @param bag.fraction the fraction of the training set observations randomly 
#' selected to propose the next tree. This introduces randomnesses into the model 
#' fit. If \code{bag.fraction<1} then running the same model twice will result in
#'  similar but different fits.  Using \code{set.seed} ensures reproducibility.
#' @param train.fraction of sample used for training
#' @param subset index of observations to use for training
#' @param indep whether random effects are independent or allowed to covary
#'  (default is TRUE, for speed)
#' @param save.mods whether the \code{lmer} models fit at each iteration are saved
#'  (required to use \code{predict})
#' @param mc.cores number of parallel cores
#' @param verbose In the final model fit, will print every `10` trees/iterations.
#' @param ... arguments passed to gbm.fit
#' @return An \code{lmerboost} object consisting of the following list elements:
#' \describe{
#'   \item{\code{yhat}}{Vector of predictions at the best iteration (\code{fixed} + \code{ranef})}
#'   \item{\code{ranef}}{Vector of random effects at the best iteration}
#'   \item{\code{fixed}}{Vector of fixed effect predictions at the best iteration}
#'   \item{\code{shrinkange}}{Amount of shrinkage}
#'   \item{\code{subset}}{Vector of observations used for training}
#'   \item{\code{best.trees}}{Best number of trees by training, test, oob, and cv error}
#'   \item{\code{best.params}}{The best set of meta-parameter values given by CV}
#'   \item{\code{params}}{A data frame of all meta-parameter combinations and the corresponding CV error}
#'   \item{\code{sigma}}{The variance due to the grouping variable at each iteration}
#'   \item{\code{xnames}}{Column names of \code{X}}
#'   \item{\code{mods}}{List of \code{lmer} models (if \code{save.mods=TRUE})}
#'   \item{\code{id}}{name or index of the grouping variable}
#'   \item{\code{trees}}{List of trees fit at each iteration}
#'   \item{\code{init}}{initial prediction}
#'   \item{\code{var.type}}{Type of variables (\code{gbm.fit})}
#'   \item{\code{c.split}}{List of categorical splits (\code{gbm.fit})}
#'   \item{\code{train.err}}{Training error at each iteration}
#'   \item{\code{oob.err}}{Out of bag error at each iteration}
#'   \item{\code{test.err}}{Test error at each iteration}
#'   \item{\code{cv.err}}{Cross-validation error at each iteration}
#' }
#' @export
lmerboost <- function(y, X, id, 
                      n.trees=5,
                      interaction.depth=3,
                      n.minobsinnode=20,
                      shrinkage=.01,
                      bag.fraction=.5,
                      train.fraction=NULL,
                      cv.folds=1,
                      subset=NULL,
                      indep=TRUE, 
                      save.mods=FALSE,
                      mc.cores=1, 
                      verbose = TRUE, ...){

  n <- length(y)
  X <- as.data.frame(X)
  params <- c(as.list(environment()),list(...)) # this won't copy y and x
  
  if(is.null(subset)) {
    train <- 1:n
  } else {
    train <- subset
  }
  if(is.logical(subset)){train <- which(subset)}
  
  if(!is.null(train.fraction)){
    train <- sample(train, ceiling(train.fraction*length(train)), replace = F)
  }
  
  if(is.character(id)){
    id <- match(id, colnames(X))
  }
  
  if(cv.folds > 1){
    #cat("cv:", fill = T)
    
    folds <- sample(1:cv.folds, size=length(train), replace = TRUE)
    
    params <- expand.grid(n.trees = n.trees, shrinkage = shrinkage, interaction.depth = interaction.depth, indep = indep)
    conds <- expand.grid(k = 1:(cv.folds), n.trees = n.trees, shrinkage = shrinkage, interaction.depth = interaction.depth, indep = indep)
    conds.ls <- split(conds, 1:nrow(conds))
    conds$id <- rep(1:nrow(params), each = cv.folds)
    
    cv.mods <- parallel::mclapply(conds.ls, function(args, ...){ 
      try(do.call(lmerboost_cv, append(args, list(...))))
    }, y=y, x=X, id=id, train=train, folds = folds, 
      bag.fraction = bag.fraction, stop.threshold = stop.threshold, 
      verbose=FALSE, save.mods=save.mods, mc.cores = mc.cores)

    # average over cv folds for each condition
    if(any(sapply(cv.mods, function(x){methods::is(x, "try-error")}))){
      return(cv.mods)
    }
    fold.err <- lapply(cv.mods, function(o){o$test.err})
    cv.err <- tapply(fold.err, conds$id, function(x){ 
      rowMeans(do.call(cbind, x), na.rm=T)})
    
    # Select the model with the lowest
    cv.err.cond <- lapply(cv.err, min, na.rm = TRUE)
    best.cond <- which.min(cv.err.cond) 
    
    params$err <- cv.err.cond
    best.params <- params[best.cond, ]
    
    cv.err <- cv.err[[best.cond]]
    best_cv_err <- which.min(cv.err)
    
    o <- lmerboost.fit(y = y, X = X, id = id, 
          train.fraction = train.fraction, subset = subset, 
          bag.fraction = bag.fraction, stop.threshold = stop.threshold, 
          verbose = verbose, save.mods=save.mods,
          n.trees = best.params$n.trees, shrinkage = best.params$shrinkage, 
          interaction.depth=best.params$interaction.depth, indep = best.params$indep,
          nt = nt)
    
  } else { # cv.folds = 1
    cv.err <- rep(NA, n.trees)
    best_cv_err <- NA
    params <- best.params <- NULL

    if(length(n.trees) > 1 | length(shrinkage) > 1 | length(interaction.depth) > 1 | length(indep) > 1){ stop("can't specify vector params without cv.folds > 1")}
    
    o <- lmerboost.fit(y = y, X = X, id = id, 
         train.fraction = train.fraction, subset = subset,
         bag.fraction = bag.fraction, stop.threshold = stop.threshold,
         n.trees = n.trees, shrinkage = shrinkage, interaction.depth=interaction.depth, indep = indep, 
         nt = nt, verbose = verbose, save.mods=save.mods)
  }
  if(all(is.na(o$test.err))){ 
    best_test_err <- NA
  } else {
    best_test_err <- which.min(o$test.err)
  }
  if(all(is.na(o$oob.err))){ 
    best_oob_err <- NA
  } else {
    best_oob_err <- which.min(o$oob.err)
  }
  best.trees <- c(train = which.min(o$train.err), test = best_test_err, oob = best_oob_err, cv = best_cv_err)
  bt <-  min(best.trees, na.rm=TRUE)
  
  out <- list(yhat=o$yhat[,bt], ranef=o$ranef[,bt], fixed=o$fixed[,bt], shrinkage=o$shrinkage, subset = subset, 
              best.trees = best.trees, best.params = best.params, params = params,
              sigma=o$sigma, xnames = colnames(X), mods=o$mods, id=id,
              trees = o$trees, init=o$init, var.type=o$var.type, c.split=o$c.split,
              train.err=o$train.err, oob.err=o$oob.err, test.err=o$test.err, cv.err=cv.err)
  class(out) <- "lmerboost"
  return(out)
}

lmerboost_cv <- function(k, folds, y, x, id, train, ...){
  cv_train <- train[folds != k]
  o <- lmerboost.fit(y = y, X = x, id = id, subset = cv_train, ...)
}



#' @describeIn lmerboost Fitting function for \code{lmerboost}
#' @param calc.derivs whether to calculate derivatives at each iteration with \code{lmer} (only assesses convergence)
#' @export
#' @importFrom stats predict
#' 
lmerboost.fit <- function(y, X, id, 
                          n.trees=5,
                          interaction.depth=3,
                          n.minobsinnode=20,
                          shrinkage=.01,
                          bag.fraction=.5,
                          train.fraction=NULL,
                          subset=NULL,
                          indep=TRUE, 
                          save.mods=FALSE,
                          verbose = TRUE, ...){

  init <- mean(y)
  r <- y - init
  n <- length(y)
  
  if(is.null(subset)) {
    train <- 1:n
  } else {
    train <- subset
  }
  if(!is.null(train.fraction)){
    train <- sample(train, ceiling(train.fraction*length(train)), replace = F)
  }
  if(is.logical(subset)){train <- which(subset)}
  
  yhat <- ranef <- fixed <- matrix(0, n, n.trees)
  sigma <- rep(0, n.trees)
  train.err <- oob.err <- test.err <- rep(NA, n.trees)
  trees <- mods <- c.split <- list()
  
  for(i in 1:n.trees){
    # s = training, s.oob = oob, -train = test
    # 2016-10-19: DO NOT STRATIFY SUBSAMPLES BY GROUPS.
    # note that just using sample(x, 1) will fail if x has length 1
    
    s <- sample(train, size = ceiling(length(train)*bag.fraction), replace = FALSE)
    s.oob <- setdiff(train, s)
    
    # fit a tree
    tree <- gbm::gbm.fit(y = r[s], x=X[s, -id, drop=F], interaction.depth=interaction.depth,
                                  shrinkage=1, bag.fraction=1, distribution="gaussian",
                                  verbose=FALSE, n.trees = 1, ...)
    if(i == 1){
      var.type = tree$var.type
    }
    trees[[i]] <- tree$trees[[1]]
    c.split[[i]] <- tree$c.splits
    pt <- gbm::pretty.gbm.tree(tree, 1)
    
    # get gbm predictions for whole sample
    gbm_pred <- predict(tree, newdata = X[,-id,drop=F], n.trees = 1) 
    
    
    # list terminal nodes (-1) first; rownames are are terminal node ids
    # this forces node factor order to have a terminal node as reference, not surrogate
    pt <- pt[order(pt$SplitVar), ]
    
    # prediction determines into which node observations fall
    # factor labels correspond to terminal node id (rows of pt)
    nodes <- droplevels(factor(gbm_pred, 
                    levels=as.character(pt$Prediction+tree$initF), 
                    labels=rownames(pt)))
    
    # design matrix - add intercept via lmer. 
    mm <- stats::model.matrix(~nodes)[,-1, drop=FALSE]
    colnames(mm) <- gsub("nodes", "X", colnames(mm))
    
    # Have to handle missinginess on id as a special case
    # observations that are missing on id will receive gbm predictions
    complete_obs <- intersect(s, which(!is.na(X[,id])))
    
    # check rank. problem is if columns are included for obs not in s via surrogates
    # dropping non-full-rank column assigns these obs to default node (arbitrary)
    # solved by: drop columns myself, replace dropped obs with gbm predictions
    
    keep_cols <- colSums(mm[complete_obs, ,drop=FALSE]) > 0
    dropped_obs  <- rowSums(mm[,!keep_cols, drop=FALSE]) > 0
    
    mm <- mm[,keep_cols, drop = FALSE]
    d <- data.frame(r=r, mm, id=X[,id])
    
    # lmer on training
    addx <- paste0(colnames(mm), collapse = " + ")
    bars <- "||"
    if(!indep) bars <- "|"
    form <- stats::as.formula(paste0("r ~ ", addx, " + (",addx, " ", bars," id)"))
    e <- new.env(parent=globalenv()) 
    e$s <- s
    environment(form) <- e
    
    o <- lme4::lmer(form, data=d, REML=T, subset = s, 
                  control = lme4::lmerControl(calc.derivs = FALSE))
    o@frame <- o@frame[1, ]
     
    if(save.mods) mods[[i]] <- o
    sigma[i] <- sigma_merMod(o) * shrinkage
    
    # 2016-10-19: Timed to show that this was fastest with large n and large ngrps
    
    yhatm <- predict(o, newdata=d, allow.new.levels = TRUE)
    fixedm <- cbind(1, mm) %*% o@beta
      
    zuhat <- yhatm - fixedm
    fixedm[dropped_obs,] <- gbm_pred[dropped_obs]
    yhatm[dropped_obs] <- gbm_pred[dropped_obs]

    # update totals at each iteration
    if(i == 1){
      fixed[,i]  <- fixedm * shrinkage 
      ranef[,i]  <- zuhat  * shrinkage
      yhat[,i]   <- yhatm  * shrinkage 
      
    } else {
      fixed[,i]  <- fixed[,i-1] + fixedm * shrinkage 
      ranef[,i]  <- ranef[,i-1] + zuhat  * shrinkage
      yhat[,i]   <- yhat[,i-1]  + yhatm  * shrinkage
    }
    
    r <- r - yhatm * shrinkage
    
    if(verbose && (i %% 10 == 0)) cat(i, "")
    if(i==1){ 
      train.err[i] <- mean((y[s] - init)^2)
      oob.err[i]   <- mean((y[s.oob] - init)^2)
      test.err[i]  <- mean((y[-train] - init)^2)
    }
    train.err[i] <- mean((yhat[s,i] - (y[s] - init))^2)
    oob.err[i]   <- mean((yhat[s.oob,i] - (y[s.oob] - init))^2)
    test.err[i]  <- mean((yhat[-train,i] - (y[-train] - init))^2)
    
    # 2016-10-19: This was removed because it can stop too early.
    #if((i %% lag == 0) && (abs(test.err[i] - test.err[i - (lag - 1)]) < stop.threshold)){
    #  break;
    #}
  }
  yhat <- yhat + init
  fixed <- fixed + init
  if(!save.mods) mods <- NULL # so you can test for it
  
  out <- list(yhat=yhat, ranef=ranef, fixed=fixed, shrinkage=shrinkage, 
              trees = trees, init=init, var.type=var.type, c.split=c.split,
              mods=mods, sigma=sigma, 
              train.err=train.err, oob.err=oob.err, test.err=test.err)
  class(out) <- "lmerboost"
  return(out)  
}


#' Prediction for lmerboost objects
#' @param object lmerboost object
#' @param newdata data frame of new data
#' @param id column name or index referring to id variable
#' @param n.trees number of trees
#' @param ... unused
#' @export 
#' @importFrom stats model.matrix
predict.lmerboost <- function(object, newdata, id, n.trees=NULL, ...){
  # save trees, lmer objects at each iteration (damn)
  if(is.null(object$mods)) stop("need to save models for predictions in newdata")
  
  if(is.null(n.trees)){ n.trees <- length(object$mods)}
  n <- nrow(newdata)
  yhat <- ranef <- fixed <- matrix(0, n, n.trees)
  if(is.character(id)){
    id <- match(id, colnames(newdata))
  }
  newid <- newdata[,id]
  newdata <- newdata[,-id, drop=F]
  
  shrinkage <- object$shrinkage
  
  for(i in 1:n.trees){
    
    pt <- data.frame(object$trees[[i]])
    names(pt) <- c("SplitVar", "SplitCodePred", "LeftNode", 
                     "RightNode", "MissingNode", "ErrorReduction", "Weight", 
                     "Prediction")
    row.names(pt) <- 0:(nrow(pt) - 1)
    
  
    # coerce lb object to a gbm object
    gbm.obj <- list(initF=object$init, trees=object$trees, 
                c.split = object$c.split, var.type=object$var.type)
    class(gbm.obj) <- "gbm"
    
    # note: when single.tree=TRUE, initF is not included in gbm predicted values
    gbm_pred <- predict(gbm.obj, n.trees=i, single.tree=TRUE, newdata=newdata) 
    
    
    # list terminal nodes (-1) first; rownames are are terminal node ids
    # this forces node factor order to have a terminal node as reference, not surrogate
    pt <- pt[order(pt$SplitVar), ]
    
    # prediction determines into which node observations fall
    # factor labels correspond to terminal node id (rows of pt)
    # note: don't want to droplevels b/c it's valid for obs to fall in just one node
    nodes <- factor(gbm_pred, 
                    levels=as.character(pt$Prediction), 
                    labels=rownames(pt))
    
    # column names of design matrix in new data match tree fit to training data
    mm <- model.matrix(~nodes)[,-1, drop=FALSE]
    colnames(mm) <- gsub("nodes", "X", colnames(mm))
    d <- data.frame(mm, id=newid)
    # no rank check because no subsampling; no dropped obs

    o <- object$mods[[i]]
    
    yhatm <- predict(o, newdata=d, allow.new.levels = TRUE)
    fixedm <- predict(o, newdata=d, re.form=NA)
    zuhat <- yhatm - fixedm

    # update totals at each iteration
    if(i == 1){
      fixed[,i]  <- fixedm * shrinkage 
      ranef[,i]  <- zuhat  * shrinkage
      yhat[,i]   <- yhatm  * shrinkage 
      
    } else {
      fixed[,i]  <- fixed[,i-1] + fixedm * shrinkage 
      ranef[,i]  <- ranef[,i-1] + zuhat  * shrinkage
      yhat[,i]   <- yhat[,i-1]  + yhatm  * shrinkage
    }
  }
  yhat <- yhat + object$init
  fixed <- fixed + object$init
  
  return(list(yhat=yhat[,n.trees], ranef=ranef[,n.trees], fixed=fixed[,n.trees]))
}

sigma_merMod <- function (object, ...) {
  dc <- object@devcomp
  dd <- dc$dims
  if (dd[["useSc"]]) 
    dc$cmp[[if (dd[["REML"]]) 
      "sigmaREML"
      else "sigmaML"]]
  else 1
}