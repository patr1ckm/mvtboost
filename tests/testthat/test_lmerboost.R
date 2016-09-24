
set.seed(104)
ngroups <- 100
group_size <- 10
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

x <- rnorm(n)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
X <- as.data.frame(x)

summary(lme4::lmer(y ~ x + (1 + x|id)))

set.seed(104)

context("test lmerboost.fit m = 1, lambda = 1, bag.fraction = 1")

o <- lmerboost.fit(y = y, X = X, id = id, 
               bag.fraction = 1,  indep = TRUE,
               M = 1, lambda = 1, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
set.seed(104)
init <- mean(y)
r = y - init
o.gbm <- gbm::gbm.fit(y = r, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
             distribution = "gaussian", interaction.depth = 5, verbose = F)
mm <- mvtboost:::gbm.mm(o.gbm, n.trees = 1, newdata = X)
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

context("test lmerboost.fit m = 1, lambda = .5, bag.fraction = 1")

lambda <- .5
o <- lmerboost.fit(y = y, X = X, id = id, 
                   bag.fraction = 1,  indep = TRUE,
                   M = 1, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
set.seed(104)
init <- mean(y)
r <- y - init
o.gbm <- gbm::gbm.fit(y = r, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
                      distribution = "gaussian", interaction.depth = 5, verbose = F)
mm <- mvtboost:::gbm.mm(o.gbm, n.trees = 1, newdata = X)
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

context("test lmerboost.fit m = 10, lambda = .5, bag.fraction = 1")

o <- lmerboost.fit(y = y, X = X, id = id, 
                   bag.fraction = 1,  indep = TRUE,
                   M = 10, lambda = lambda, depth = 5, stop.threshold = 0, n.minobsinnode = 10)

M <- 10
zuhat <- fixed <- yhat <- matrix(0, n, 10)
init <- mean(y)
r <- y - mean(y)
for(i in 1:M){
  set.seed(104)
  o.gbm <- gbm::gbm.fit(y = r, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
                        distribution = "gaussian", interaction.depth = 5, verbose = F)
  mm <- mvtboost:::gbm.mm(o.gbm, n.trees = 1, newdata = X)
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
