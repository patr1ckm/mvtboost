## RRV: 2016-09-21
## - Added best_iter, a function to compute the best iteration in various ways
## RRV: 2016-09-20
## - Added early stopping (via stop.threshold)
## - Made it 10x faster by computing the random effects Z %*% u fast


# y vector of observations
# X matrix of predictors
# id grouping variable
# indep whether the random effects of nodes are correlated or indep (default FALSE for speed)
# M is number of random effects estimation steps
# lambda is the step size
# cv.folds = number of cv.folds
# nt = number of trees used in gbm/lmer iterations
# depth = depth of trees used in gbm/lmer iteration

# vectors can be passed to M, lambda, indep, and depth for cross validation

#' boosting with random effects
#' 
#' @param y outcome vector (continuous)
#' @param X matrix or data frame of predictors 
#' @param id name or index of grouping variable
#' @param train.fraction of sample used for training
#' @param subset index of observations to use for training
#' @export
#' @importFrom parallel mclapply
lmerboost <- function(y, X, id, 
                      train.fraction=NULL, 
                      subset=NULL, 
                      bag.fraction=.5, 
                      cv.folds=1,
                      indep=TRUE, 
                      M=100, 
                      lambda=.01, 
                      nt=1, 
                      depth=5, 
                      tune=FALSE, 
                      stop.threshold = .001,
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
  if(!is.null(train.fraction)){
    train <- sample(train, ceiling(train.fraction*length(ss)), replace = F)
  }
  if(is.logical(subset)){train <- which(subset)}
  if(is.character(id)){
    id <- match(id, colnames(X))
  }
  
  if(cv.folds > 1){
    #cat("cv:", fill = T)
    
    folds <- sample(1:cv.folds, size=length(train), replace = TRUE)
    
    params <- expand.grid(M = M, lambda = lambda, depth = depth, indep = indep)
    conds <- expand.grid(k = 1:(cv.folds), M = M, lambda = lambda, depth = depth, indep = indep)
    conds.ls <- split(conds, 1:nrow(conds))
    conds$id <- rep(1:nrow(params), each = cv.folds)
    
    cv.mods <- parallel::mclapply(conds.ls, function(args, ...){ 
      do.call(lmerboost_cv, append(args, list(...)))
    }, y=y, x=X, id=id, train=train, folds = folds, 
      bag.fraction = bag.fraction, stop.threshold = stop.threshold,
      verbose = FALSE, mc.cores = mc.cores)

    # average over cv folds for each condition
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
          bag.fraction = bag.fraction, stop.threshold = stop.threshold, verbose = verbose,
          M = best.params$M, lambda = best.params$lambda, depth=best.params$depth, indep = best.params$indep,
          nt = nt,  tune=FALSE)
    
  } else { # cv.folds = 1
    cv.err <- rep(NA, M)
    best_cv_err <- NA
    params <- best.params <- NULL

    if(length(M) > 1 | length(lambda) > 1 | length(depth) > 1 | length(indep) > 1){ stop("can't specify vector params without cv.folds > 1")}
    
    o <- lmerboost.fit(y = y, X = X, id = id, 
         train.fraction = train.fraction, subset = subset,
         bag.fraction = bag.fraction, stop.threshold = stop.threshold,
         M = M, lambda = lambda, depth=depth, indep = indep, 
         nt = nt,  tune=FALSE, verbose = verbose)
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
  
  out <- list(yhat=o$yhat[,bt], ranef=o$ranef[,bt], fixed=o$fixed[,bt], lambda=o$lambda, subset = subset, 
              best.trees = best.trees, best.params = best.params, params = params,
              sigma=o$sigma, xnames = colnames(X), mods=o$mods,
              trees = o$trees, init=o$init, var.type=o$var.type, c.split=o$c.split,
              train.err=o$train.err, oob.err=o$oob.err, test.err=o$test.err, cv.err=cv.err)
  class(out) <- "lmerboost"
  return(out)
}

lmerboost_cv <- function(k, folds, y, x, id, train, ...){
  cv_train <- train[folds != k]
  o <- lmerboost.fit(y = y, X = x, id = id, subset = cv_train, ...)
}


#' @export
#' @importFrom gbm gbm.fit pretty.gbm.tree
#' @importFrom lme4 lmer
lmerboost.fit <- function(y, X, id, train.fraction=NULL, subset=NULL, indep=TRUE, M=100, 
                          lambda=.01, nt=1, depth=5, tune=FALSE, bag.fraction=.5, 
                          calc.derivs=FALSE, stop.threshold = .001, verbose = TRUE, ...){
  
  init <- mean(y)
  r <- y - init
  n <- length(y)
  
  if(is.null(subset)) {
    train <- 1:n
  } else {
    train <- subset
  }
  if(!is.null(train.fraction)){
    train <- sample(train, ceiling(train.fraction*length(ss)), replace = F)
  }
  if(is.logical(subset)){train <- which(subset)}
  
  yhat <- ranef <- fixed <- matrix(0, n, M)
  sigma <- rep(0, M)
  train.err <- oob.err <- test.err <- rep(NA, M)
  trees <- mods <- c.split <- list()
  
  for(i in 1:M){
    # s = training, s.oob = oob, -train = test
    # 2016-10-19: DO NOT STRATIFY SUBSAMPLES BY GROUPS.
    # note that just using sample(x, 1) will fail if x has length 1
    
    s <- sample(train, size = ceiling(length(train)*bag.fraction), replace = FALSE)
    s.oob <- setdiff(train, s)
    
    # fit a tree
    tree <- gbm.fit(y = r[s], x=X[s, -id, drop=F], interaction.depth=depth,
                                  shrinkage=1, bag.fraction=1, distribution="gaussian",
                                  verbose=FALSE, n.trees = 1, ...)
    if(i == 1){
      var.type = tree$var.type
    }
    trees[[i]] <- tree$trees[[1]]
    c.split[[i]] <- tree$c.splits
    pt <- gbm::pretty.gbm.tree(tree, 1)
    
    # get gbm predictions for whole sample
    gbm_pred <- predict(tree, newdata = X, n.trees = 1) 
    
    
    # list terminal nodes (-1) first; rownames are are terminal node ids
    # this forces node factor order to have a terminal node as reference, not surrogate
    pt <- pt[order(pt$SplitVar), ]
    
    # prediction determines into which node observations fall
    # factor labels correspond to terminal node id (rows of pt)
    nodes <- droplevels(factor(gbm_pred, 
                    levels=as.character(pt$Prediction+tree$initF), 
                    labels=rownames(pt)))
    
    # design matrix - add intercept via lmer. 
    mm <- model.matrix(~nodes)[,-1, drop=FALSE]
    colnames(mm) <- gsub("nodes", "X", colnames(mm))
    
    # check rank. problem is if columns are included for obs not in s via surrogates
    # dropping non-full-rank column assigns these obs to default node.
    # solved by: drop columns myself, replace dropped obs with gbm predictions
    keep_cols <- colSums(mm[s, ,drop=FALSE]) > 0
    dropped_obs  <- rowSums(mm[,!keep_cols, drop=FALSE]) > 0
    
    mm <- mm[,keep_cols, drop = FALSE]
    d <- data.frame(r=r, mm, id=X[,id])
    
    # lmer on training
    addx <- paste0(colnames(mm), collapse = " + ")
    bars <- "||"
    if(!indep) bars <- "|"
    form <- as.formula(paste0("r ~ ", addx, " + (",addx, " ", bars," id)"))
    
    mods[[i]] <- o <- lme4::lmer(form, data=d, REML=T, subset = s, 
                    control = lme4::lmerControl(calc.derivs = FALSE))
    sigma[i] <- sigma_merMod(o) * lambda
    
    # 2016-10-19: Timed to show that this was fastest with large n and large ngrps
    
    yhatm <- predict(o, newdata=d, allow.new.levels = TRUE)
    fixedm <- cbind(1, mm) %*% o@beta
      
    zuhat <- yhatm - fixedm
    fixedm[dropped_obs,] <- gbm_pred[dropped_obs]
    yhatm[dropped_obs] <- gbm_pred[dropped_obs]

    # update totals at each iteration
    if(i == 1){
      fixed[,i]  <- fixedm * lambda 
      ranef[,i]  <- zuhat  * lambda
      yhat[,i]   <- yhatm  * lambda 
      
    } else {
      fixed[,i]  <- fixed[,i-1] + fixedm * lambda 
      ranef[,i]  <- ranef[,i-1] + zuhat  * lambda
      yhat[,i]   <- yhat[,i-1]  + yhatm  * lambda
    }
    
    r <- r - yhatm * lambda
    
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
  
  out <- list(yhat=yhat, ranef=ranef, fixed=fixed, lambda=lambda, 
              trees = trees, init=init, var.type=var.type, c.split=c.split,
              mods=mods, sigma=sigma, 
              train.err=train.err, oob.err=oob.err, test.err=test.err)
  class(out) <- "lmerboost"
  return(out)  
}


#' @export
predict.lmerboost <- function(object, newdata, newid, M=NULL){
  # save trees, lmer objects at each iteration (damn)
  
  if(is.null(M)){ M <- length(object$mods)}
  n <- nrow(newdata)
  yhat <- ranef <- fixed <- matrix(0, n, M)
  
  lambda <- object$lambda
  
  for(i in 1:M){
    
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
      fixed[,i]  <- fixedm * lambda 
      ranef[,i]  <- zuhat  * lambda
      yhat[,i]   <- yhatm  * lambda 
      
    } else {
      fixed[,i]  <- fixed[,i-1] + fixedm * lambda 
      ranef[,i]  <- ranef[,i-1] + zuhat  * lambda
      yhat[,i]   <- yhat[,i-1]  + yhatm  * lambda
    }
  }
  yhat <- yhat + object$init
  fixed <- fixed + object$init
  
  return(list(yhat=yhat[,M], ranef=ranef[,M], fixed=fixed[,M]))
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