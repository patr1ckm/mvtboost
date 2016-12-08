
# These tests are as follows
# 1. Make sure it runs 
# 2. Compare current version to a previous implementation; make sure changes are ok
# 3. Make sure all the buttons work

set.seed(104)
ngroups <- 100
group_size <- 10
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

train <- sample(1:800, size = 500, replace = FALSE)

x <- rnorm(n)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
X <- data.frame(x, id)
tol = 1E-6

# helper function used in older version

context("lmerboost.fit")

test_that("lmerboost runs", {
  o1 <- lmerboost(y = y, X = X, id=2, M = 5, cv.folds = 1, lambda = .1)
  o2 <- lmerboost(y = y, X = X, id="id", M = 5, cv.folds = 1, lambda = .1)
  o <- lmerboost(y = y, X = X, id=2, M = 5, cv.folds = 1, lambda = .1, subset = train)
  o <- lmerboost(y = y, X = X, id=2, M = 3, cv.folds = 3)
  Xmis <- X
  Xmis[sample(1:n, size = n/2, replace = FALSE),] <- NA
  o <- lmerboost(y = y, X = Xmis, id=2, M = 3, cv.folds = 3)
  expect_is(o, "lmerboost")
})


test_that("lmerboost.fit bag.fraction, lambda, subset, train/oob/test err", {
  bag.fraction = .5
  lambda <- .5
  n <- length(y)

  set.seed(104)
  o <- lmerboost.fit(y = y, X = X, id=2, subset = train,
                     bag.fraction = bag.fraction,  indep = TRUE, verbose = FALSE,
                     M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)


  set.seed(104)
  idx <- 2
  id <- X[,idx]
  M <- 10
  zuhat <- fixed <- yhat <- matrix(0, n, M)
  train_err <- oob_err <- test_err <- rep(0, M)
  init <- mean(y)
  r <- y - mean(y)
  for(i in 1:M){

    # the only change is to subsample from train rather than 1:n, and to subset on id
    # note that s is always an index to observations in the  original data
    s <- sample(train, size=ceiling(length(train)*bag.fraction), replace=FALSE)
    s.oob <- setdiff(train, s)

    o.gbm <- gbm::gbm.fit(y = r[s], x = X[s, -idx, drop = F], n.trees = 1, shrinkage = 1, bag.fraction = 1,
                          distribution = "gaussian", interaction.depth = 5, verbose = F)
    pt <- gbm::pretty.gbm.tree(o.gbm, i.tree=1)
    gbm_pred <- predict(o.gbm, newdata=X[,-idx, drop=F], n.trees=1)
    nodes <- droplevels(factor(gbm_pred, levels=as.character(pt$Prediction + o.gbm$initF), 
                               labels=1:nrow(pt)))
    mm <- model.matrix(~nodes)[,-1]
    colnames(mm) <- gsub("nodes", "X", colnames(mm))
    
    keep_cols <- colSums(mm[s, ]) > 0
    dropped_obs <- rowSums(mm[, !keep_cols, drop=F]) > 0
    mm <- mm[,keep_cols]
    
    d <- data.frame(r, mm, id)
    addx <- paste0(colnames(mm), collapse=" + ")
    form <- as.formula(paste0("r ~ ", addx, " + (", addx, "||id)"))
    o.lmer <- lme4::lmer(form, data=d, subset=s, REML = T,
                         control = lme4::lmerControl(calc.derivs = FALSE))
    Zm <- model.matrix(~id + mm:id - 1)
    
    re <- as.matrix(lme4::ranef(o.lmer)[[1]]) #
    new_re <- as.data.frame(matrix(0, nrow = length(unique(id)), ncol = ncol(re)))
    new_re[rownames(re), ] <- re
    new_re <- as.matrix(new_re)

    zuhatm <- drop(Zm %*% c(new_re))
    fixedm <- drop(cbind(1, mm) %*% lme4::fixef(o.lmer))
    fixedm[dropped_obs] <- gbm_pred[dropped_obs]
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
    test_err[i]  <- mean(((y[-train, ] - init) - yhat[-train, i])^2)
    
  }
  yhat <- yhat + init
  fixed <- fixed + init

  expect_equal(yhat[train, ], o$yhat[train, ])
  expect_equal(zuhat[train, ], o$ranef[train, ])
  expect_equal(fixed[train, ], o$fixed[train, ])
  expect_equal(yhat[-train, ], o$yhat[-train, ])
  expect_equal(zuhat[-train, ], o$ranef[-train, ])
  expect_equal(fixed[-train, ], o$fixed[-train, ])
  expect_equal(train_err, o$train.err)  
  expect_equal(oob_err, o$oob.err)  
  expect_equal(test_err, o$test.err)  
})

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
  
  X <- data.frame(x, x2, id)
  lb <- lmerboost.fit(y=y, X=X, id=3, lambda=1, M=1, depth=1, bag.fraction=1,
                        n.minobsinnode=1, subset=train)
  expect_equal(lb$yhat[train], unname(yhatm[train]))
  expect_equal(lb$ranef[train], unname(zuhat[train]))
  expect_equal(lb$fixed[train], unname(fixedm[train]))
  expect_equal(lb$trees, og$trees)
  expect_equal(lb$yhat[-train], unname(yhatm[-train]))
  expect_equal(lb$ranef[-train], unname(zuhat[-train]))
  expect_equal(lb$fixed[-train], unname(fixedm[-train]))
  
})
## TODO: lmerboost.fit logical subset

## TODO: lmerboost.fit train.fraction, stop.threshold, depth, indep

context("lmerboost")

test_that("lmerboost_cv", {
  # now we can use lmerboost.fit
  cv.folds <- 3
  folds <- sample(1:cv.folds, size=n, replace=TRUE)
  
  set.seed(104)
  ocv <- lapply(1:cv.folds, function(k, folds, train){
    ss <- train[folds != k]
    lmerboost.fit(y = y, X = X, id=2, subset = ss, M = 5, verbose = F, lambda = .5)
  }, folds = folds, train=1:n)
  
  set.seed(104)
  o <- lapply(1:cv.folds, mvtboost:::lmerboost_cv, 
           folds = folds, train = 1:n, y = y, x = X, id=2, M = 5, lambda = .5)
  
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
                 folds = folds, y = y, x = x, id=2, train = train , ...)}
  
  ocv <- lapply(X = paramscv.ls, FUN=do_one, folds=folds, train=1:n, y=y, x=X, id=2, 
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
  o <- lmerboost(y = y, X = X, id=2, cv.folds = 3, bag.fraction = 1, subset = 1:n,
                lambda = c(.2, .5), depth = c(3, 5), M = 5, mc.cores = 1)
  
  expect_identical(bc, o$best.params)
  expect_identical(params, o$params)
  
})


test_that("lmerboost err", {
  # error in lmerboost.fit
  y[56] <- NA
  msg <- capture_messages(
    expect_error(o1 <- lmerboost(y = y, X = X, id=2, M = 5, cv.folds = 1, lambda = .1))
  )
  
  msg <- capture_messages(
    o1 <- lmerboost(y = y, X = X, id=2, M = 5, cv.folds = 3, lambda = .1, mc.cores=3)
  )
  expect_true(all(sapply(o1, function(o){is(o, "try-error")})))
  
})

## TODO: checks of train.fraction, logical subset, etc

test_that("lmerboost_cv train", {
  # for now, just test that they run
  o <- lmerboost(y=y, X=X, id=2, cv.folds=3, bag.fraction=1, subset=1:(n/2), M=3)
  o <- lmerboost(y=y, X=X, id=2, cv.folds=3, bag.fraction=.5, subset=1:(n/2), M=3)
})

## TODO: combinations of params

test_that("lmerboost influence", {
  X <- data.frame(X1 = x, X2 = rnorm(n), id=id)
  ob <- lmerboost(y = y, X = X, id=3, M = 3, cv.folds = 1, lambda = .5)
  inf <- influence(ob)
  expect_gt(inf[1], 0)
  expect_equal(length(inf), 2)
})

test_that("lmerboost predict", {
  M = 5
  lb <- lmerboost(y = y, X = X, id=2, M = M, cv.folds = 3, lambda = c(.5, 1),
                 bag.fraction=.5, subset=1:500, save.mods=TRUE)
  yh <- predict(lb, newdata=X, id=2, M=min(lb$best.trees, na.rm=T))
  expect_equal(lb$yhat, yh$yhat)
  expect_equal(lb$fixed, yh$fixed)
  expect_equal(lb$ranef, yh$ranef)
  
  Xnew <- data.frame(x=rnorm(n), id=id)
  Xnewid <- data.frame(x=rnorm(n), id=factor(rep(101, n)))
  
  yh2 <- predict(lb, newdata=Xnew, id=2)
  
  yh3 <- predict(lb, newdata=Xnewid, id=2)
  # a new group has 0 ranef
  expect_equal(yh3$ranef, rep(0, n))
  
})
