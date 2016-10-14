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

#' @export
#' @importFrom parallel mclapply
lmerboost <- function(y, X, id, 
                      train.fraction=NULL, subset=NULL, bag.fraction=.5, cv.folds=1,
                      indep=TRUE, M=100, lambda=.01, nt=1, depth=5, 
                      tune=FALSE, stop.threshold = .001,
                      mc.cores=1, verbose = TRUE, ...){

  n <- length(y)
  X <- as.data.frame(X)
  params <- c(as.list(environment()),list(...)) # this won't copy y and x
  
  if(is.null(subset)) {
    ss <- 1:n
  } else {
    ss <- subset
  }
  if(!is.null(train.fraction)){
    ss <- sample(ss, ceiling(train.fraction*length(ss)), replace = F)
  }
  if(is.logical(subset)){ss <- which(subset)}
  
  
  new.levels <- any(!(id %in% id[ss]))
  
  if(cv.folds > 1){
    #cat("cv:", fill = T)
    
    folds <- assign_fold(ss, id = id, cv.folds = cv.folds)
    
    params <- expand.grid(M = M, lambda = lambda, depth = depth, indep = indep)
    conds <- expand.grid(k = 1:(cv.folds), M = M, lambda = lambda, depth = depth, indep = indep)
    conds.ls <- split(conds, 1:nrow(conds))
    conds$id <- rep(1:nrow(params), each = cv.folds)
    
    cv.mods <- parallel::mclapply(conds.ls, function(args, ...){ 
      do.call(lmerboost_cv, append(args, list(...)))
    }, y=y, x=X, id=id, ss=ss, folds = folds, 
      bag.fraction = bag.fraction, stop.threshold = stop.threshold, new.levels = new.levels,
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
          train.fraction = train.fraction, subset = subset, new.levels = new.levels,
          bag.fraction = bag.fraction, stop.threshold = stop.threshold, verbose = verbose,
          M = best.params$M, lambda = best.params$lambda, depth=best.params$depth, indep = best.params$indep,
          nt = nt,  tune=FALSE)
    
  } else { # cv.folds = 1
    cv.err <- rep(NA, M)
    best_cv_err <- NA
    params <- best.params <- NULL

    if(length(M) > 1 | length(lambda) > 1 | length(depth) > 1 | length(indep) > 1){ stop("can't specify vector params without cv.folds > 1")}
    
    o <- lmerboost.fit(y = y, X = X, id = id, 
         train.fraction = train.fraction, subset = subset, new.levels = new.levels,
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
  
  out <- list(yhat = o$yhat, ranef = o$ranef, fixed = o$fixed, lambda = o$lambda, subset = subset, 
              yhatt = o$yhatt, raneft = o$raneft, fixedt = o$fixedt,
              best.trees = best.trees,
              cond.cv.err = params, best.params = best.params, params = params,
              trees = o$trees, xnames = colnames(X),
              train.err=o$train.err, oob.err=o$oob.err, test.err=o$test.err, cv.err=cv.err, 
              s = ss)
  class(out) <- "lmerboost"
  return(out)
}

lmerboost_cv <- function(k, folds, y, x, id, ss, new.levels, ...){
  train <- ss[folds != k]
  o <- lmerboost.fit(y = y, X = x, id = id, subset = train, ...)
}


#' @export
lmerboost.fit <- function(y, X, id, train.fraction=NULL, subset=NULL, indep=TRUE, M=100, 
                          lambda=.01, nt=1, depth=5, tune=FALSE, bag.fraction=.5, 
                          calc.derivs=FALSE, stop.threshold = .001, new.levels = TRUE, verbose = TRUE, ...){
  
  init <- mean(y)
  r <- y - init
  lag <- 5
  
  n <- length(y)
  
  if(is.null(subset)) {
    ss <- 1:n
  } else {
    ss <- subset
  }
  if(!is.null(train.fraction)){
    ss <- sample(ss, ceiling(train.fraction*length(ss)), replace = F)
  }
  if(is.logical(subset)){ss <- which(subset)}
  
  yhat <- ranef <- fixed <- matrix(0, n, M)
  train.err <- oob.err <- test.err <- rep(NA, M)
  trees <- list()
  
  
  for(i in 1:M){
    # s = training, s.oob = oob, -ss = test
    # make sure one observation from each group is present, where the subset has to cover all groups
    #  note that just using sample(x, 1) will fail if x has length 1
    if(bag.fraction < 1){
      s <- get_subsample(ss = ss, id = droplevels(id[ss]), bag.fraction = bag.fraction)
    } else {
      s <- ss
    }
    
    s.oob <- setdiff(ss, s) # the oob observations are the observations in training but not in subset
    datx <- data.frame(r=r[s], X[s, ])
    colnames(datx) <- c("r", colnames(X))
    
    tree <- gbm::gbm(r ~ ., data=datx, n.trees=nt, shrinkage=1, distribution="gaussian", 
                            bag.fraction=1, interaction.depth=depth, ...)
    
    trees[[i]] <- tree$trees[[1]]
    
    if(tune & (nt > 2)){ # errors with nt = 1:2
      tr <- suppressWarnings(gbm.perf(tree, method="OOB", plot.it=F))
    } else {
      tr <- nt
    }
    
    # Get model matrix for train, oob and test
    # In this formulation, surrogate splits are separate nodes.
    yhat_gbm <- predict(o, n.trees=nt)
    node <- factor(yhat_gbm)
    
    # mm <- model.matrix(~node)[,-1,drop=F] default coding adding intercept in lmer
    # cell means coding. works if tree doesn't split and allows easy rank check
    mm <- model.matrix(~node - 1) 
    
    # Quick rank check with missing data:
    #  drop nodes with no observations in the training set. Only occurs with NAs in X,
    #   not gauranteed to have the same patterns of missing values in training and test.
    #  if columns are not dropped, the model matrix is rank deficient 
    #   and lmer drops the column with only a warning, breaking predictions.
    # otherwise, node assignment is assumed to be a full rank operation.
    dropped_cols <- apply(mm[s,], 2, function(col){length(unique(col))}) == 1
    
    # These are the observations in the training + test affected by dropped columns
    # The predictions for these observations will be replaced with predictions from gbm.
    # This side-steps th issue by treating surrogate nodes not in training as fixed.
    dropped_obs <- which(apply(mm[,dropped_cols, drop=F], 1, function(row){ row > 0}))
    mm <- mm[,!dropped_cols, drop = FALSE]
    
    nodes <- ncol(mm)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    dat.mm <- data.frame(r=r, id=id, mm)
    
    # note that in cell means, intercept is dropped for fixed but not random
    if(indep){
      form <- as.formula(paste0("r ~ 0 + ", paste0("X", 1:nodes, collapse = " + "),  " + (",
                                paste0("X", 1:nodes, collapse=" + "), " + 1 || id)"))
    } else {
      form <- as.formula(paste0("r ~ 0 + ", paste0("X", 1:nodes, collapse = " + "),  " + (",
                                paste0("X", 1:nodes, collapse=" + "), " + 1 | id)"))
    }
    
    
    
    # Fit lmer model to training set only, not oob or test
    out.lmer <- lme4::lmer(form, data=dat.mm, REML=T, subset = s, 
                          control = lme4::lmerControl(calc.derivs = calc.derivs))
    dat.mm$r <- NULL
    
    
    ## for the ids in the test but not in training, augment re with 0
    re <- as.matrix(lme4::ranef(out.lmer)[[1]]) #
    if(new.levels){
      # Note that this step might be able to be optimized...? 
      new_re <- as.data.frame(matrix(0, nrow = length(unique(id)), ncol = ncol(re)))
      rownames(new_re) <- unique(as.character(id))
      new_re[rownames(re), ] <- re
      new_re <- as.matrix(new_re)
    } else {
      new_re <- re
    }
    
    # Get random, fixed, and total for train, oob, and test at iteration m
    zuhat <- get_zuhat(new_re, x = mm, id = id)
    #fixedm <- cbind(1, mm) %*% as.matrix(lme4::fixef(out.lmer))
    fixedm <- mm %*% as.matrix(lme4::fixef(out.lmer))
    fixedm[dropped_obs, ] <- yhat_gbm[dropped_obs]
    yhatm <- fixedm + zuhat

    # update totals at each iteration
    if(i == 1){
      fixed[,i]  <- fixedm * lambda 
      ranef[,i]  <- zuhat * lambda
      yhat[,i]   <- yhatm * lambda 
      
    } else {
      fixed[,i]  <- fixed[,i-1] + fixedm * lambda 
      ranef[,i]  <- ranef[,i-1] + zuhat * lambda
      yhat[,i]   <- yhat[,i-1] + yhatm * lambda
    }
    
    r <- r - yhatm * lambda
    
    if(verbose && (i %% 10 == 0)) cat(i, "")
    if(i==1){ 
      train.err[i] <- mean((y[s] - init)^2)
      oob.err[i] <- mean((y[s.oob] - init)^2)
      test.err[i] <- mean((y[-ss] - init)^2)
    }
    train.err[i] <- mean((yhat[s,i] - (y[s] - init))^2)
    oob.err[i] <- mean((yhat[s.oob,i] - (y[s.oob] - init))^2)
    test.err[i] <- mean((yhat[-ss,i] - (y[-ss] - init))^2)
    
    # This was removed because it can stop too early.
    #if((i %% lag == 0) && (abs(test.err[i] - test.err[i - (lag - 1)]) < stop.threshold)){
    #  break;
    #}
  }
  yhat <- yhat + init
  fixed <- fixed + init
  
  out <- list(yhat=yhat[ss, ], ranef=ranef[ss, ], fixed=fixed[ss,], lambda=lambda, 
              yhatt=yhat[-ss, ], raneft=ranef[-ss, ], fixedt=fixed[-ss, ],
              trees = trees,
              train.err=train.err, oob.err=oob.err, test.err=test.err)
  return(out)  
}


# re = ranef(mod)$id
# x = design matrix without intercept
# id = grouping variable factor
get_zuhat <- function(re, x, id){
  Z <- model.matrix(~id + id:x - 1)
  b <- c(re)
  drop(Z %*% b)
}


gbm_mm <- function(o, n.trees=1, ...){
  yhat <- predict(o, n.trees=n.trees, ...)
  node <- factor(yhat)
  model.matrix(~node)[,-1,drop=F]
}

assign_fold <- function(x, id, cv.folds){
  folds <- list()
  ids_by_group <- split(1:length(x), id[x])
  
  assign_fold_1group <- function(x, folds){
    n <- length(x)
    if(n == 1){
      in_fold <- 0 # observation will always be in training set
    } else {
      # If 1 < n_i <= cv.folds, force one observation from group to be in training set
      #  then randomly assign fold ids as usual for the other observations
      in_fold <- c(sample(1:folds, size=min(n, folds), replace = F),
                   sample(1:folds, size=max(0, (n - folds)), replace = T))
    }
    return(in_fold)
  }
  
  folds <- unlist(lapply(ids_by_group, assign_fold_1group, folds = cv.folds), use.names = F)

  return(folds)
}


## This is a hard problem! The current implementation doesn't really admit a separate prediction method.
## without re-fitting the model.
predict.lmerboost <- function(object, newdata, M=NULL){
  # save trees, lmer objects at each iteration (damn)
  
  if(is.null(M)){ M <- length(object$mods)}
  n <-  length(object$yhat)
  yhat <- ranef <- fixed <- matrix(0, n, M)
  
  lambda <- out$lambda
  
  for(m in 1:M){
  
    mm <- gbm_mm(object$trees[[m]], newdata=newdata, n.trees = 1)
    nodes <- ncol(mm)
    out.lmer <- object$mods[[m]] 
    new.mm <- data.frame(id=id, mm)
    
    zucoefs <- as.matrix(ranef(out.lmer)[[1]])
    zuhat <- get_zuhat(zucoefs, mm, id)
    
    yhatm <- predict(out.lmer, newdata=newdata)
    
    ranef[,i] <- ranef[,i-1] + zuhat * lambda
    
    yhat[,i] <- yhat[,i-1] + yhatm * lambda
    
    fixed[,i] <- fixed[,i-1] + (new.mm[,-grep("id", colnames(new.mm))] %*% as.matrix(fixef(out.lmer))) * lambda
  }
  
  return(list(yhat=yhat, ranef=ranef, fixed=fixed))
}


best_iter <- function(x, threshold, lag, smooth = FALSE){
  err <- x$cv.err
  err <- err[!is.na(err)]
  if(smooth) err <- smooth(err)
  
  best.iter <- which(abs(diff(err, lag = lag)) < threshold)
  
  if(length(best.iter) == 0){
    best.iter <- which.min(x$cv.err)
  } else {
    best.iter <- min(best.iter)
  }
}


get_subsample <- function(ss, id, bag.fraction){
  id <- droplevels(id)
  
  # Get one obs from each group
  sid <- tapply(ss, id, function(x){ 
    x[sample(1:length(x), size=1)]
  })
  # compute the number of samples
  size <- ceiling(length(ss) * bag.fraction) - length(sid)
  
  # sample randomly from the rest
  c(sid, sample(ss[!(ss %in% sid)], size = size, replace=F)) 
}

#' @export
plot.lmerboost <- function(x, threshold = .001, lag = 1, ...){
  M <- length(x$train.err)
  ymax <- c(max(x$test.err, x$train.err, x$oob.err, na.rm = T))
  ymin <- c(min(x$test.err, x$train.err, x$oob.err, na.rm = T))
  
  best.iter <- best_iter(x, threshold = threshold, lag = lag)
  
  plot(x = 1:M, y = x$train.err, type = "l", ylim = c(ymin, ymax), ylab = "error")
  lines(x = 1:M, y = x$test.err, col = "red", lty = 2)
  lines(x = 1:M, y = x$oob.err, col = "blue")
  lines(x = 1:M, y = x$cv.err, col = "red")
  
  abline(v = best.iter)
  
  legend("top", legend = c("train", "test", "oob", "cv"), 
         col = c("black", "red", "blue", "red"),
         lty = c(1, 2, 1, 1), bty = "n")
  x$best.params$err <- formatC(signif(x$best.params$err[[1]], digits=3), digits=3,format="fg", flag="#")
  paramstring <- paste0(names(x$best.params), " = ", x$best.params, collapse = ", ")
  title(sub = paramstring)
}

