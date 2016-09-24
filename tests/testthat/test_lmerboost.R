
set.seed(104)
ngroups <- 100
group_size <- 10
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

x <- rnorm(n)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
y <- y - mean(y)
X <- as.data.frame(x)

summary(lme4::lmer(y ~ x + (1 + x|id)))

set.seed(104)
o <- lmerboost.fit(y = y, X = X, id = id, 
               bag.fraction = 1,  indep = TRUE,
               M = 1, lambda = 1, depth = 5, stop.threshold = 0, n.minobsinnode = 10)
set.seed(104)
o.gbm <- gbm::gbm.fit(y = y, x = X, n.trees = 1, shrinkage = 1, bag.fraction = 1, 
             distribution = "gaussian", interaction.depth = 5)
mm <- mvtboost:::gbm.mm(o.gbm, n.trees = 1, newdata = X)
colnames(mm) <- paste0("X", 1:ncol(mm))
d <- data.frame(y, mm, id)
o.lmer <- lme4::lmer(y ~ X1 + X2 + X3 + X4 + X5 + (1 + X1 + X2 + X3 + X4 + X5 || id), REML = T, 
                     control = lme4::lmerControl(calc.derivs = FALSE), data = d)

Zm <- model.matrix(~id + mm:id - 1)
zuhat <- drop(Zm %*% c(as.matrix(ranef(o.lmer)[[1]])))
fixed <- drop(cbind(1, mm) %*% fixef(o.lmer))

expect_equivalent(zuhat, o$ranef)
expect_equivalent(fixed, o$fixed)
expect_equivalent(unname(fixed + zuhat), o$yhat)
expect_equivalent(o$yhat, unname(predict(o.lmer)))
