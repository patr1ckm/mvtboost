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
lmerboost <- function(y, X, id, 
                      train.fraction=NULL, subset=NULL, bag.fraction=.5, cv.folds=1,
                      indep=TRUE, M=100, lambda=.01, nt=1, depth=5, 
                      tune=FALSE, stop.threshold = .001,
                      mc.cores=1, ...){

  n <- length(y)
  X <- as.data.frame(X)
  params <- c(as.list(environment()),list(...)) # this won't copy y and x
  if(cv.folds > 1){
    #cat("cv:", fill = T)
    
    if(is.null(subset)) {
      ss <- 1:n
    } else {
      ss <- subset
    }
    if(!is.null(train.fraction)){
      ss <- sample(ss, ceiling(train.fraction*length(ss)), replace = F)
    }
    if(is.logical(subset)){ss <- which(subset)}
    
    folds <- list()
    idg <- split(1:length(ss), id[ss]) #ids by group
    
    folds <- unlist(sapply(idg, assign_fold, folds = cv.folds))
    
    cv.tune <- function(k, folds, y, x, id, ss, ...){
      train <- ss[folds != k]
      test <- ss[folds == k]
      o <- lmerboost.fit(y=y, X=x, id=id, subset = train, ...)
    }
    
    params <- expand.grid(M = M, lambda = lambda, depth = depth, indep = indep)
    conds <- expand.grid(k = 1:(cv.folds), M = M, lambda = lambda, depth = depth, indep = indep)
    conds.ls <- split(conds, 1:nrow(conds))
    conds$id <- rep(1:nrow(params), each = cv.folds)
    
    cv.mods <- mclapply(conds.ls, function(args, ...){ 
      do.call(cv.tune, append(args, list(...)))
    }, y=y, x=X, id=id, ss=ss, folds = folds, stop.threshold = stop.threshold, mc.cores = mc.cores)

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
    
    o <- lmerboost.fit(y = y, X = X, id = id, 
          train.fraction = train.fraction, subset = subset, 
          bag.fraction = bag.fraction, stop.threshold = stop.threshold,
          M = best.params$M, lambda = best.params$lambda, depth=best.params$depth, indep = best.params$indep,
          nt = nt,  tune=FALSE)
    
  } else { # cv.folds = 1
    cv.err <- rep(NA, M)
    params = NULL
    if(length(M) > 1 | length(lambda) > 1 | length(depth) > 1 | length(indep) > 1){ stop("can't specify vector params without cv.folds > 1")}
    
    o <- lmerboost.fit(y = y, X = X, id = id, 
         train.fraction = train.fraction, subset = subset, 
         bag.fraction = bag.fraction, stop.threshold = stop.threshold,
         M = M, lambda = lambda, depth=depth, indep = indep,
         nt = nt,  tune=FALSE)
  }
  
  cat("", fill=T)
  out <- list(yhat=o$yhat, ranef=o$ranef, fixed=o$fixed, lambda=o$lambda, subset=subset, 
              yhatt=o$yhatt, raneft=o$raneft, fixedt=o$fixedt,
              cond.cv.err = params, best.params = best.params, params = params,
              train.err=o$train.err, oob.err=o$oob.err, test.err=o$test.err, cv.err=cv.err, 
              s = ss)
  class(out) <- "lmerboost"
  return(out)
}

#' @export
lmerboost.fit <- function(y, X, id, train.fraction=NULL, subset=NULL, indep=TRUE, M=100, 
                          lambda=.01, nt=1, depth=5, tune=FALSE, bag.fraction=.5, 
                          calc.derivs=FALSE, stop.threshold = .001, ...){
  
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
  
  for(i in 1:M){
    # s = training, s.oob = oob, -ss = test
    # make sure one observation from each group is present, where the subset has to cover all groups
    #  note that just using sample(x, 1) will fail if x has length 1
    if(bag.fraction < 1){
      s <- get_subsample(ss = ss, id = id, bag.fraction = bag.fraction)
    } else {
      s <- ss
    }
    s.oob <- setdiff(ss, s) # the oob observations are the observations in training but not in subset
    datx <- data.frame(r=r[s], X[s, ])
    colnames(datx) <- c("r", colnames(X))
    tree <- gbm::gbm(r ~ ., data=datx, n.trees=nt, shrinkage=1, distribution="gaussian", 
                     bag.fraction=1, interaction.depth=depth, ...)
    if(tune & (nt > 2)){ # errors with nt = 1:2
      tr <- suppressWarnings(gbm.perf(tree, method="OOB", plot.it=F))
    } else {
      tr <- nt
    }
    
    # Get model matrix or train, oob and test
    mm <- gbm.mm(tree, n.trees = tr, newdata = X)
    
    nodes <- ncol(mm)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    dat.mm <- data.frame(r=r, id=id, mm)
  
    if(indep){
      form <- as.formula(paste0("r ~ ", paste0("X", 1:nodes, collapse = " + "),  " + (",
                                paste0("X", 1:nodes, collapse=" + "), " + 1 || id)"))
    } else {
      form <- as.formula(paste0("r ~ ", paste0("X", 1:nodes, collapse = " + "),  " + (",
                                paste0("X", 1:nodes, collapse=" + "), " + 1 | id)"))
    }
    
    # Fit lmer model to training set only, not oob or test
    out.lmer <- lme4::lmer(form, data=dat.mm, REML=T, subset = s, 
                          control = lme4::lmerControl(calc.derivs = calc.derivs))
    dat.mm$r <- NULL
    
    zucoefs <- as.matrix(lme4::ranef(out.lmer)[[1]]) # 
    
    # Get random, fixed, and total for train, oob, and test at iteration m
    zuhat <- get.zuhat(zucoefs, x = mm, id = dat.mm$id)
    fixedm <- cbind(1, mm) %*% as.matrix(lme4::fixef(out.lmer))
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
    
    cat(i, "")
    if(i==2){ 
      # for the first iteration, compute the training, oob, and test error from only
      # the means
      train.err[i-1] <- var(y[s])
      oob.err[i-1] <- var(y[s.oob])
      test.err[i-1] <- var(y[-ss])
    }
    train.err[i] <- var(yhat[s,i] - y[s])
    oob.err[i] <- var(yhat[s.oob,i] - y[s.oob])
    test.err[i] <- var(yhat[-ss,i] - y[-ss])
    
    if(i %% lag == 0 & abs(test.err[i] - test.err[i - (lag - 1)]) < stop.threshold){
      break;
    }
  }
  yhat <- yhat + init
  
  
  out <- list(yhat=yhat[ss, ], ranef=ranef[ss, ], fixed=fixed[ss,], lambda=lambda, 
              yhatt=yhat[-ss, ], raneft=ranef[-ss, ], fixedt=fixed[-ss, ],
              train.err=train.err, oob.err=oob.err, test.err=test.err)
  return(out)  
}


# re = ranef(mod)$id
# x = design matrix without intercept
# id = grouping variable factor
get.zuhat <- function(re, x, id){
  Z <- model.matrix(~id + id:x - 1)
  b <- c(re)
  drop(Z %*% b)
}


gbm.mm <- function(o, n.trees=1, ...){
  yhat <- predict(o, n.trees=n.trees, ...)
  node <- factor(yhat)
  model.matrix(~node)[,-1,drop=F]
}

# x can be an index or value
assign_fold <- function(x, folds){
  n <- length(x)
  in_fold <- 0 # if n_i = 1, that observation will always be in training set
  if(n > 1){ 
    # If 1 < n_i < cv.folds, randomly assign a fold id to each i w/o replacement
    # then randomly assign fold ids to the rest w/ replacement
    in_fold <- c(sample(1:folds, size=min(n, folds), replace = F),
                 sample(1:folds, size=max(0, (n - folds)), replace = T))
  }
  return(in_fold)
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
  
    mm <- gbm.mm(object$trees[[m]], newdata=newdata, n.trees = 1)
    nodes <- ncol(mm)
    out.lmer <- object$mods[[m]] 
    new.mm <- data.frame(id=id, mm)
    
    zucoefs <- as.matrix(ranef(out.lmer)[[1]])
    zuhat <- get.zuhat(zucoefs, mm, id)
    
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
  # Get one obs from each group
  sid <- tapply(ss, droplevels(id[ss]), function(x){ 
    x[sample(1:length(x), size=1)]
  })
  # compute the number of samples
  size <- ceiling(length(ss) * bag.fraction) - length(sid)
  
  # sample randomly from the rest
  c(sid, sample(ss[!(ss %in% sid)], size = size, replace=F)) 
}

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

## Testing code based on evaluation
#put.args(lmerboost)

