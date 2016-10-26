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
  
  if(cv.folds > 1){
    #cat("cv:", fill = T)
    
    folds <- sample(1:cv.folds, size=n, replace = TRUE)
    
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
  
  out <- list(yhat = o$yhat, ranef = o$ranef, fixed = o$fixed, lambda = o$lambda, subset = subset, 
              yhatt = o$yhatt, raneft = o$raneft, fixedt = o$fixedt,
              best.trees = best.trees,
              cond.cv.err = params, best.params = best.params, params = params,
              trees = o$trees, sigma=o$sigma, xnames = colnames(X),
              train.err=o$train.err, oob.err=o$oob.err, test.err=o$test.err, cv.err=cv.err, 
              s = train)
  class(out) <- "lmerboost"
  return(out)
}

lmerboost_cv <- function(k, folds, y, x, id, train, ...){
  cv_train <- train[folds != k]
  o <- lmerboost.fit(y = y, X = x, id = id, subset = cv_train, ...)
}


#' @export
#' @importFrom gbm gbm.fit
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
  trees <- list()
  
  for(i in 1:M){
    # s = training, s.oob = oob, -train = test
    # 2016-10-19: DO NOT STRATIFY SUBSAMPLES BY GROUPS.
    # note that just using sample(x, 1) will fail if x has length 1
    
    s <- sample(train, size = ceiling(length(train)*bag.fraction), replace = FALSE)
    s.oob <- setdiff(train, s)
    
    # fit a tree
    tree <- gbm.fit(y = r[s], x=X[s, ,drop=F], interaction.depth=depth,
                                  shrinkage=1, bag.fraction=1, distribution="gaussian",
                                  verbose=FALSE, n.trees = 1, ...)
    trees[[i]] <- tree$trees[[1]]
    # get gbm predictions for whole sample
    gbm_pred <- predict(tree, newdata = X, n.trees = 1) 
    
    # design matrix - add intercept via lmer. 
    # 2016-10-19: BE VERY CAREFUL IF YOU CHANGE THIS
    mm <- model.matrix(~factor(gbm_pred))[,-1]
    nnodes <- ncol(mm)
    colnames(mm) <- paste0("X", 1:nnodes)
    
    # formula
    addx <- paste0(colnames(mm), collapse = " + ")
    bars <- "||"
    if(!indep) bars <- "|"
    
    form <- as.formula(paste0("r ~ ", addx, " + (",addx, " ", bars," id)"))
    
    # check rank. problem is if columns are included for obs not in s via surrogates
    # dropping non-full-rank column assigns these obs to default node.
    # solved by: drop columns myself, replace dropped obs with gbm predictions
    keep_cols <- colSums(mm[s, ]) > 0
    dropped_obs  <- rowSums(mm[,!keep_cols, drop=FALSE]) > 0
    
    mm <- mm[,keep_cols, drop = FALSE]
    d <- data.frame(r=r, mm, id)
    
    # lmer on training
    o <- lme4::lmer(form, data=d, REML=T, subset = s, 
                    control = lme4::lmerControl(calc.derivs = FALSE))
    sigma[i] <- sigma(o) * lambda
    
    # 2016-10-19: Timed to show that this was fastest with large n and large ngrps
    yhatm <- predict(o, newdata=d, allow.new.levels = TRUE)
    fixedm <- cbind(1, mm) %*% o@beta
    zuhat <- yhatm - fixedm
    fixedm[dropped_obs, ] <- gbm_pred[dropped_obs]

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
  
  out <- list(yhat=yhat[train, ], ranef=ranef[train, ], fixed=fixed[train,], lambda=lambda, 
              yhatt=yhat[-train, ], raneft=ranef[-train, ], fixedt=fixed[-train, ],
              trees = trees, sigma=sigma, init=init,
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

