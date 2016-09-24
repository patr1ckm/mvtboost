
set.seed(104)
ngroups <- 100
group_size <- 10
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

train <- unlist(tapply(1:n, id, function(x){
  x[sample(1:length(x), size = group_size * .5, replace = F)]}))

x <- rnorm(n)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
X <- as.data.frame(x)

# summary(lme4::lmer(y ~ x + (1 + x|id)))

context("lmerboost.fit")

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
  d <- data.frame(r, mm, id)
  o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                       control = lme4::lmerControl(calc.derivs = FALSE), data = d)
  
  Zm <- model.matrix(~id + mm:id - 1)
  zuhat <- drop(Zm %*% c(as.matrix(ranef(o.lmer)[[1]])))
  fixed <- drop(cbind(1, mm) %*% fixef(o.lmer))
  yhat <- fixed + zuhat + init
  
  expect_equal(unname(zuhat), o$ranef)
  expect_equal(unname(fixed), o$fixed)
  expect_equal(unname(yhat), o$yhat)
  expect_equal(o$yhat, unname(predict(o.lmer)) + init)
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
  zuhat <- drop(Zm %*% c(as.matrix(ranef(o.lmer)[[1]]))) * lambda
  fixed <- drop(cbind(1, mm) %*% fixef(o.lmer)) * lambda
  yhat <- fixed + zuhat + init
  
  expect_equal(unname(zuhat), o$ranef)
  expect_equal(unname(fixed), o$fixed)
  expect_equal(unname(fixed + zuhat + init), o$yhat)
  expect_equal(o$yhat, unname(predict(o.lmer) * lambda + init))
})

test_that("lmerboost.fit m = 10, lambda = .5, bag.fraction = 1", {
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
    zuhatm <- drop(Zm %*% c(as.matrix(ranef(o.lmer)[[1]])))
    fixedm <- drop(cbind(1, mm) %*% fixef(o.lmer))
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
  
  expect_equal(yhat, o$yhat)
  expect_equal(zuhat, o$ranef)
  expect_equal(fixed, o$fixed)
})

test_that("lmerboost.fit get_subsample", {
  bag.fraction <- .5
  s <- mvtboost:::get_subsample(1:n, id, bag.fraction = bag.fraction)
  expect_equal(length(s), ceiling(n * bag.fraction))
  expect_equal(length(unique(s)), length(s))
  expect_true(all(unique(id) %in% unique(id[s]))) # at least one from each group occurs
})

test_that("lmerboost.fit bag.fraction = .5", {
  bag.fraction = .5
  set.seed(104)
  o <- lmerboost.fit(y = y, X = X, id = id, 
                     bag.fraction = bag.fraction,  indep = TRUE, verbose = FALSE,
                     M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
  set.seed(104)
  M <- 10
  zuhat <- fixed <- yhat <- matrix(0, n, 10)
  init <- mean(y)
  r <- y - mean(y)
  for(i in 1:M){
    
    s <- mvtboost:::get_subsample(1:n, id = id, bag.fraction = bag.fraction)
    o.gbm <- gbm::gbm.fit(y = r[s], x = X[s, , drop = F], n.trees = 1, shrinkage = 1, bag.fraction = 1, 
                          distribution = "gaussian", interaction.depth = 5, verbose = F)
    mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    d <- data.frame(r, mm, id)
    o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                         control = lme4::lmerControl(calc.derivs = FALSE), data = d, subset = s)
    
    Zm <- model.matrix(~id + mm:id - 1)
    zuhatm <- drop(Zm %*% c(as.matrix(ranef(o.lmer)[[1]])))
    fixedm <- drop(cbind(1, mm) %*% fixef(o.lmer))
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
  
  expect_equal(yhat, o$yhat)
  expect_equal(zuhat, o$ranef)
  expect_equal(fixed, o$fixed)
})

test_that("lmerboost.fit subset, train/oob/test err", {
  bag.fraction = .5
  # get a training sample from each group of group_size * .5
  
  expect_true(all(unique(id) %in% unique(id[train])))

  set.seed(104)
  o <- lmerboost.fit(y = y, X = X, id = id, subset = train,
                     bag.fraction = bag.fraction,  indep = TRUE, verbose = FALSE,
                     M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)


  set.seed(104)
  M <- 10
  zuhat <- fixed <- yhat <- matrix(0, n, M)
  train_err <- oob_err <- test_err <- rep(0, M)
  init <- mean(y)
  r <- y - mean(y)
  for(i in 1:M){

    # the only change is to subsample from train rather than 1:n, and to subset on id
    # note that s is always an index to observations in the  original data
    s <- mvtboost:::get_subsample(train, id = id[train], bag.fraction = bag.fraction)
    s.oob <- setdiff(train, s)
    expect_true(all(unique(id) %in% unique(id[s])))

    o.gbm <- gbm::gbm.fit(y = r[s], x = X[s, , drop = F], n.trees = 1, shrinkage = 1, bag.fraction = 1,
                          distribution = "gaussian", interaction.depth = 5, verbose = F)
    mm <- mvtboost:::gbm_mm(o.gbm, n.trees = 1, newdata = X)
    colnames(mm) <- paste0("X", 1:ncol(mm))
    d <- data.frame(r, mm, id)
    o.lmer <- lme4::lmer(r ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T,
                         control = lme4::lmerControl(calc.derivs = FALSE), data = d, subset = s)

    Zm <- model.matrix(~id + mm:id - 1)
    zuhatm <- drop(Zm %*% c(as.matrix(ranef(o.lmer)[[1]])))
    fixedm <- drop(cbind(1, mm) %*% fixef(o.lmer))
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
    train_err[i] <- var(y[s, ] - yhat[s, i])
    oob_err[i]   <- var(y[s.oob, ] - yhat[s.oob, i])
    test_err[i]  <- var(y[-train, ] - yhat[-train, i])
    
  }
  yhat <- yhat + init

  expect_equal(yhat[train, ], o$yhat)
  expect_equal(zuhat[train, ], o$ranef)
  expect_equal(fixed[train, ], o$fixed)
  expect_equal(yhat[-train, ], o$yhatt)
  expect_equal(zuhat[-train, ], o$raneft)
  expect_equal(fixed[-train, ], o$fixedt)
  expect_equal(train_err, o$train.err)  
  expect_equal(oob_err, o$oob.err)  
  expect_equal(test_err, o$test.err)  
})

## TODO: lmerboost.fit logical subset

test_that("lmerboost.fit get_zuhat", {
  o <- lme4::lmer(y ~ x + (1 + x|id))
  re <- as.matrix(ranef(o)$id)
  zuhat <- mvtboost:::get_zuhat(re = re, x = x, id = id)
  zuhat_lmer <- c(predict(o) - cbind(1, x) %*% fixef(o))
  expect_equal(unname(zuhat), zuhat_lmer)
})

## TODO: lmerboost.fit train.fraction, stop.threshold, depth, indep

context("lmerboost")

test_that("lmerboost assign_fold", {
  train <- 1:n
  for(cv.folds in c(3, 5, 10)){
    folds <- mvtboost:::assign_fold(train, id = id, cv.folds = cv.folds)
    expect_equal(length(folds), n)
    for(k in 1:cv.folds){
      expect_true(all(unique(id) %in% unique(id[train[folds != k]])))
      expect_true(all(unique(id) %in% unique(id[train[folds == k]]))) 
    }
  }
  
  # Have observations within group = 1, 2, ..., cv.folds, ... 
  id_short <- factor(unlist(lapply(1:10, function(x){rep(x,x)})))
  n_short <- length(id_short)
  train <- seq_along(id_short)
  for(cv.folds in c(3, 5, 10)){
    folds <- mvtboost:::assign_fold(train, id = id_short, cv.folds = cv.folds)
    expect_equal(length(folds), n_short)
    for(k in 1:cv.folds){
      expect_true(all(unique(id_short) %in% unique(id_short[train[folds != k]])))
    }
  }
  
})

test_that("lmerboost_cv", {
  # now we can use lmerboost.fit
  
  cv.folds <- 3
  folds <- mvtboost:::assign_fold(train, id = id, cv.folds = cv.folds)
  
  set.seed(104)
  ocv <- lapply(1:cv.folds, function(k, folds){
    ss <- (1:n)[folds != k]
    lmerboost.fit(y = y, X = X, id = id, subset = ss, M = 5, verbose = F, lambda = .5)
  }, folds = folds)
  
  set.seed(104)
  o <- lapply(1:cv.folds, mvtboost:::lmerboost_cv, 
           folds = folds, ss = 1:n, y = y, x = X, id = id, M = 5, lambda = .5)
  
  expect_equivalent(o, ocv)
})





