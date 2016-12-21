context("twostage")
set.seed(104)
ngroups <- 50
group_size <- 3
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

train <- sample(1:n, size=n/2, replace=T)

x <- rnorm(n)
xc <- cut(x, 2)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
X <- data.frame(x, xc, id=id)

test_that("twostage runs and predicts", {
  n.trees <- 3
  o <- twostage(y=y, x=X, id="id",
                n.trees=n.trees, 
                shrinkage=c(.005, 01),
                interaction.depth=1,
                cv.folds=3,
                subset=train,
                distribution="gaussian", verbose=FALSE)
  expect_true(class(o) == "twostage")
  
  d2 <- X[,-match("id", colnames(X)), drop=F]
  
  yhat <- predict(o, newdata=X, id="id", n.trees=n.trees)
  
  ans <- predict(o$o.lmer, newdata=X, allow.new.levels = TRUE) + 
         predict(o$o.gbm, newdata=d2, n.trees=n.trees)
  
  expect_equal(yhat, ans)
})

test_that("twostage influence drops id", {
  ri <- influence(o)  
  expect_true(!all(grepl("id", names(ri))))
})

test_that("missing ids work", {
  X[1:3,"id"] <- NA
  o <- twostage(y=y, x=X, id="id",
                n.trees=2, 
                subset=train,
                distribution="gaussian", verbose=FALSE)
  expect_true(class(o) == "twostage")
  yhat <- predict(o, newdata=X, id="id", n.trees=2)
  expect_equal(length(yhat), length(y))
})



