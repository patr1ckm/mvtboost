
set.seed(104)
ngroups <- 100
group_size <- 10
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

traiN <- sample(1:800, size = 500, replace = FALSE)

x <- rnorm(n)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
X <- as.data.frame(x)
tol = 5E-7

context("lmerboost.fit")

test_that("lmerboost runs", {
  o <- lmerboost(y = y, X = X, id = id, M = 5, cv.folds = 1, lambda = .1)
  o <- lmerboost(y = y, X = X, id = id, M = 5, cv.folds = 1, lambda = .1, subset = traiN)
  o <- lmerboost(y = y, X = X, id = id, M = 3, cv.folds = 3)
  Xmis <- X
  Xmis[sample(1:n, size = n/2, replace = FALSE),] <- NA
  o <- lmerboost(y = y, X = Xmis, id = id, M = 3, cv.folds = 3)
  expect_is(o, "lmerboost")
})

test_that("lmerboost.fit m = 1, lambda = 1, bag.fraction = 1", {
  set.seed(104)
  o <- lmerboost.fit(y = y, X = X, id = id, 
                 bag.fraction = 1,  indep = TRUE,
                 M = 1, lambda = 1, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
  set.seed(104)
  init <- mean(y)
  r = y - init
  o.gbm <- gbm::gbm.fit(y = r, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
               distribution = "gaussian", interaction.depth = 5, verbose = F)
  mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
  colnames(mm) <- paste0("X", 1:ncol(mm))
  dr <- data.frame(r, mm, id)
  o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                       control = lme4::lmerControl(calc.derivs = FALSE), data = dr)
  
  Zm <- model.matrix(~id + mm:id - 1)
  zuhat <- drop(Zm %*% c(as.matrix(lme4::ranef(o.lmer)[[1]])))
  fixed <- drop(cbind(1, mm) %*% lme4::fixef(o.lmer)) + init
  yhat <- fixed + zuhat
  
  expect_equal(init, o$init)
  expect_equal(o.gbm$trees[[1]], o$trees[[1]])
  expect_equal(unname(zuhat), o$ranef, tolerance = tol)
  expect_equal(unname(fixed), o$fixed, tolerance = tol)
  expect_equal(unname(yhat), o$yhat, tolerance = tol)
  expect_equal(unname(predict(o.lmer)) + init, o$yhat, tolerance = tol)
})

test_that("lmerboost.fit m = 1, lambda = .5, bag.fraction = 1", {
  lambda <- .5
  o <- lmerboost.fit(y = y, X = X, id = id, 
                     bag.fraction = 1,  indep = TRUE,
                     M = 1, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
  set.seed(104)
  init <- mean(y)
  r <- y - init
  o.gbm <- gbm::gbm.fit(y = r, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
                        distribution = "gaussian", interaction.depth = 5, verbose = F)
  mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
  colnames(mm) <- paste0("X", 1:ncol(mm))
  d <- data.frame(r, mm, id)
  o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                       control = lme4::lmerControl(calc.derivs = FALSE), data = d)
  
  Zm <- model.matrix(~id + mm:id - 1)
  zuhat <- drop(Zm %*% c(as.matrix(lme4::ranef(o.lmer)[[1]]))) * lambda
  fixed <- drop(cbind(1, mm) %*% lme4::fixef(o.lmer)) * lambda + init
  yhat <- fixed + zuhat 
  
  expect_equal(unname(zuhat), o$ranef, tolerance = tol)
  expect_equal(unname(fixed), o$fixed, tolerance = tol)
  expect_equal(unname(fixed + zuhat), o$yhat, tolerance = tol)
  expect_equal(unname(predict(o.lmer))*lambda + init, o$yhat, tolerance = tol)
})

test_that("lmerboost.fit m = 10, lambda = .5, bag.fraction = 1", {
  lambda <- .5
  o <- lmerboost.fit(y = y, X = X, id = id, 
                     bag.fraction = 1,  indep = TRUE, verbose = FALSE,
                     M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
  
  M <- 10
  zuhat <- fixed <- yhat <- matrix(0, n, 10)
  init <- mean(y)
  r <- y - mean(y)
  for(i in 1:M){
    set.seed(104)
    o.gbm <- gbm::gbm.fit(y = r, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
                          distribution = "gaussian", interaction.depth = 5, verbose = F)
    mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    d <- data.frame(r, mm, id)
    o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                         control = lme4::lmerControl(calc.derivs = FALSE), data = d)
    
    Zm <- model.matrix(~id + mm:id - 1)
    zuhatm <- drop(Zm %*% c(as.matrix(lme4::ranef(o.lmer)[[1]])))
    fixedm <- drop(cbind(1, mm) %*% lme4::fixef(o.lmer))
    yhatm <- zuhatm + fixedm 
    if(i == 1){
      zuhat[,i] <- zuhatm * lambda
      fixed[,i] <- fixedm * lambda
      yhat[,i] <- yhatm * lambda
    } else {
      zuhat[,i] <- zuhat[,i - 1] + zuhatm * lambda
      fixed[,i] <- fixed[,i - 1] + fixedm * lambda
      yhat[,i] <- yhat[,i-1] + yhatm * lambda
    }
    r <- r - yhatm * lambda
  }
  yhat <- yhat + init
  fixed <- fixed + init
  
  expect_equal(yhat, o$yhat, tolerance = tol)
  expect_equal(zuhat, o$ranef,  tolerance = tol)
  expect_equal(fixed, o$fixed,  tolerance = tol)
})

test_that("lmerboost.fit bag.fraction = .5", {
  bag.fraction = .5
  lambda <- .5
  set.seed(104)
  o <- lmerboost.fit(y = y, X = X, id = id, 
                     bag.fraction = bag.fraction,  indep = TRUE, verbose = FALSE,
                     M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
  set.seed(104)
  M <- 10
  n <- length(y)
  zuhat <- fixed <- yhat <- matrix(0, n, 10)
  init <- mean(y)
  r <- y - mean(y)
  for(i in 1:M){
    
    s <- sample(1:n, size=ceiling(bag.fraction*n), replace=FALSE)
    o.gbm <- gbm::gbm.fit(y = r[s], x = X[s, , drop = F], n.trees = 1, shrinkage = 1, bag.fraction = 1, 
                          distribution = "gaussian", interaction.depth = 5, verbose = F)
    mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    d <- data.frame(r, mm, id)
    o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                         control = lme4::lmerControl(calc.derivs = FALSE), data = d, subset = s)
    
    Zm <- model.matrix(~id + mm:id - 1)
    zuhatm <- drop(Zm %*% c(as.matrix(lme4::ranef(o.lmer)[[1]])))
    fixedm <- drop(cbind(1, mm) %*% lme4::fixef(o.lmer))
    yhatm <- zuhatm + fixedm 
    if(i == 1){
      zuhat[,i] <- zuhatm * lambda
      fixed[,i] <- fixedm * lambda
      yhat[,i] <- yhatm * lambda
    } else {
      zuhat[,i] <- zuhat[,i - 1] + zuhatm * lambda
      fixed[,i] <- fixed[,i - 1] + fixedm * lambda
      yhat[,i] <- yhat[,i-1] + yhatm * lambda
    }
    r <- r - yhatm * lambda
  }
  yhat <- yhat + init
  fixed <- fixed + init
  
  expect_equal(yhat, o$yhat,  tolerance = tol)
  expect_equal(zuhat, o$ranef, tolerance = tol)
  expect_equal(fixed, o$fixed, tolerance = tol)
})

test_that("lmerboost.fit subset, train/oob/test err", {
  bag.fraction = .5
  lambda <- .5
  n <- length(y)

  set.seed(104)
  o <- lmerboost.fit(y = y, X = X, id = id, subset = traiN,
                     bag.fraction = bag.fraction,  indep = TRUE, verbose = FALSE,
                     M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)


  set.seed(104)
  M <- 10
  zuhat <- fixed <- yhat <- matrix(0, n, M)
  train_err <- oob_err <- test_err <- rep(0, M)
  init <- mean(y)
  r <- y - mean(y)
  for(i in 1:M){

    # the only change is to subsample from traiN rather than 1:n, and to subset on id
    # note that s is always an index to observations in the  original data
    s <- sample(traiN, size=ceiling(length(traiN)*bag.fraction), replace=FALSE)
    s.oob <- setdiff(traiN, s)

    o.gbm <- gbm::gbm.fit(y = r[s], x = X[s, , drop = F], n.trees = 1, shrinkage = 1, bag.fraction = 1,
                          distribution = "gaussian", interaction.depth = 5, verbose = F)
    mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    d <- data.frame(r, mm, id)
    o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T,
                         control = lme4::lmerControl(calc.derivs = FALSE), data = d, subset = s)

    Zm <- model.matrix(~id + mm:id - 1)
    
    re <- as.matrix(lme4::ranef(o.lmer)[[1]]) #
    new_re <- as.data.frame(matrix(0, nrow = length(unique(id)), ncol = ncol(re)))
    new_re[rownames(re), ] <- re
    new_re <- as.matrix(new_re)

    zuhatm <- drop(Zm %*% c(new_re))
    fixedm <- drop(cbind(1, mm) %*% lme4::fixef(o.lmer))
    yhatm <- zuhatm + fixedm
    if(i == 1){
      zuhat[,i] <- zuhatm * lambda
      fixed[,i] <- fixedm * lambda
      yhat[,i] <- yhatm * lambda
    } else {
      zuhat[,i] <- zuhat[,i - 1] + zuhatm * lambda
      fixed[,i] <- fixed[,i - 1] + fixedm * lambda
      yhat[,i] <- yhat[,i-1] + yhatm * lambda
    }
    r <- r - yhatm * lambda
    train_err[i] <- mean(((y[s, ] - init) - yhat[s, i])^2)
    oob_err[i]   <- mean(((y[s.oob, ] - init) - yhat[s.oob, i])^2)
    test_err[i]  <- mean(((y[-traiN, ] - init) - yhat[-traiN, i])^2)
    
  }
  yhat <- yhat + init
  fixed <- fixed + init

  expect_equal(yhat[traiN, ], o$yhat)
  expect_equal(zuhat[traiN, ], o$ranef)
  expect_equal(fixed[traiN, ], o$fixed)
  expect_equal(yhat[-traiN, ], o$yhatt)
  expect_equal(zuhat[-traiN, ], o$raneft)
  expect_equal(fixed[-traiN, ], o$fixedt)
  expect_equal(train_err, o$train.err)  
  expect_equal(oob_err, o$oob.err)  
  expect_equal(test_err, o$test.err)  
})

## TODO:

test_that("lmerboost.fit drops rank deficient cols", {
  # Rank deficiency in training can occur with missing data due to surrogate splitting.
  # gbm always defines a surrogate, and in the full data the surrogate
  #  might be used to make a prediction. 
  
  # Since the design matrix is created for unique predictions, the full
  # data (training + test) might have more unique predictions than the training.
  # The training data might have a column where no observations fall (all 0s)
   
  # Since lmer is fit to training data, the design matrix is rank deficient.
  # lmer drops the column with a warning.
  
  # However, this means that the model matrix of the full data and of training
  # don't have the same dimensions; which breaks predictions.
  
  # Test case that demonstrates missing nodes
  set.seed(104)
  x <- rep(0:1, each = 10)
  x2 <- sample(x)
  y <- x + rnorm(20)   # will split on x1, not x2
  x[c(10, 20)] <- NA   # will use x2 for a surrogate for these obs
  train <- c(1:9, 11:19)   # training set does not have missing values
  id <- factor(rep(1:5, each=4))
  lambda = 1
  bag.fraction=1
  i <- 1
  X <- cbind(x, x2)
  s <- sample(train, size = ceiling(length(train)*bag.fraction), replace = FALSE)
  s.oob <- setdiff(train, s)
  
  init <- mean(y)  
  r <- y - init
  d <- data.frame(r, x, x2)
  og <- gbm::gbm(r ~ ., data=d[s,], n.minobsinnode=1, shrinkage=1, n.trees=1, 
           distribution="gaussian", interaction.depth=1, bag.fraction=1)
  gbm_pred <- predict(og, n.trees=1, newdata=d)
  mm <- model.matrix(~factor(gbm_pred))[,-1,drop=F]
  
  keep_cols <- colSums(mm[s,,drop=FALSE ]) > 0
  dropped_obs <- rowSums(mm[,!keep_cols,drop=FALSE]) > 0
  mm <- mm[,keep_cols, drop = FALSE]
  colnames(mm) <- paste0("X", 1:ncol(mm))
  addx <- paste0(colnames(mm), collapse = "+")
  form <- as.formula(paste0("r ~ 1 + ", addx, "+ (1 + ", addx, " | id)"))
  
  # lmer on training
  d <- data.frame(r, mm, id)
  o <- lme4::lmer(form, data=d, REML=T, subset = s, 
                  control = lme4::lmerControl(calc.derivs = FALSE))
  
  # 2016-10-19: Timed to show that this was fastest with large n and large ngrps
  yhatm <- predict(o, newdata=d, allow.new.levels = TRUE)
  fixedm <- cbind(1, mm) %*% o@beta
  zuhat <- yhatm - fixedm
  fixedm[dropped_obs, ] <- gbm_pred[dropped_obs]
  yhatm[dropped_obs] <- gbm_pred[dropped_obs]
  yhatm <- yhatm + init
  fixedm <- fixedm + init
   
  set.seed(104)
  lb <- lmerboost.fit(y=y, X=X, id=id, lambda=1, M=1, depth=1, bag.fraction=1,
                        n.minobsinnode=1, subset=train)
  expect_equal(lb$yhat, unname(yhatm[train]))
  expect_equal(lb$ranef, unname(zuhat[train]))
  expect_equal(lb$fixed, unname(fixedm[train]))
  expect_equal(lb$trees, og$trees)
  expect_equal(lb$yhatt, unname(yhatm[-train]))
  expect_equal(lb$raneft, unname(zuhat[-train]))
  expect_equal(lb$fixedt, unname(fixedm[-train]))
  
})

## TODO: lmerboost.fit logical subset

test_that("lmerboost.fit get_zuhat", {
  o <- lme4::lmer(y ~ x + (1 + x|id))
  re <- as.matrix(lme4::ranef(o)$id)
  zuhat <- mvtboost:::get_zuhat(re = re, x = x, id = id)
  zuhat_lmer <- c(predict(o) - cbind(1, x) %*% lme4::fixef(o))
  expect_equal(unname(zuhat), zuhat_lmer)
})

## TODO: lmerboost.fit train.fraction, stop.threshold, depth, indep

context("lmerboost")

test_that("lmerboost_cv", {
  # now we can use lmerboost.fit

  cv.folds <- 3
  folds <- sample(1:cv.folds, size=n, replace=TRUE)
  
  set.seed(104)
  ocv <- lapply(1:cv.folds, function(k, folds, train){
    ss <- train[folds != k]
    lmerboost.fit(y = y, X = X, id = id, subset = ss, M = 5, verbose = F, lambda = .5)
  }, folds = folds, train=1:n)
  
  set.seed(104)
  o <- lapply(1:cv.folds, mvtboost:::lmerboost_cv, 
           folds = folds, train = 1:n, y = y, x = X, id = id, M = 5, lambda = .5)
  
  expect_equivalent(o, ocv)
})

test_that("lmerboost cv params", {
  
  set.seed(104)
  cv.folds = 3
  folds <- sample(1:cv.folds, size=n, replace=TRUE)
  paramscv <- expand.grid(M = 5, k = 1:cv.folds, lambda = c(.2, .5), depth = c(3, 5), indep = TRUE)
  params <- expand.grid(M = 5, lambda = c(.2, .5), depth = c(3, 5), indep = TRUE)
  paramscv$id <- factor(rep(1:nrow(params), each = cv.folds))
  paramscv.ls <- split(paramscv, 1:nrow(paramscv))
  do_one <- function(args, folds, y, x, id, train, ...){ 
    mvtboost:::lmerboost_cv(k = args$k, depth = args$depth, lambda = args$lambda,
                 folds = folds, y = y, x = x, id = id, train = train , ...)}
  
  ocv <- lapply(X = paramscv.ls, FUN=do_one, folds=folds, train=1:n, y=y, x=X, id=id, 
                M = 5, bag.fraction = 1)

  fold.err <- lapply(ocv, function(o){o$test.err})
  cv.err <- tapply(fold.err, paramscv$id, function(x){ 
    rowMeans(do.call(cbind, x), na.rm=T)})
  
  
  min.cv.err <- lapply(cv.err, min)
  params$err <- min.cv.err
  cv.err.cond <- lapply(cv.err, min, na.rm = TRUE)
  best.cond <- which.min(cv.err.cond)
  bc <- params[best.cond, ]
  
  set.seed(104)
  o <- lmerboost(y = y, X = X, id = id, cv.folds = 3, bag.fraction = 1, subset = 1:n,
                lambda = c(.2, .5), depth = c(3, 5), M = 5, mc.cores = 1)
  
  expect_identical(bc, o$best.params)
  expect_identical(params, o$params)
  
})



## TODO: checks of train.fraction, logical subset, etc

test_that("lmerboost_cv train", {
  # for now, just test that they run
  o <- lmerboost(y=y, X=X, id=id, cv.folds=3, bag.fraction=1, subset=1:(n/2), M=3)
  o <- lmerboost(y=y, X=X, id=id, cv.folds=3, bag.fraction=.5, subset=1:(n/2), M=3)
})

## TODO: combinations of params

test_that("lmerboost influence", {
  X <- data.frame(X1 = x, X2 = rnorm(n))
  ob <- lmerboost(y = y, X = X, id = id, M = 3, cv.folds = 1, lambda = .5)
  inf <- influence(ob)
  expect_gt(inf[1], 0)
})
