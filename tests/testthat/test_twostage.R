

set.seed(104)
ngroups <- 50
group_size <- 3
n <- ngroups * group_size
id <- factor(rep(1:ngroups, each = group_size))

train <- sample(1:n, size=n/2, replace=T)

x <- rnorm(n)
Z <- model.matrix(~id + x:id - 1)
u <- rnorm(ncol(Z), 0, 1)
y <- x * .5 + Z %*% u + rnorm(n)
X <- as.data.frame(x)

test_that("twostage runs", {
  o <- twostage(y=y, x=X, id=id,
                n.trees=5, 
                shrinkage=c(.005, 01), 
                cv.folds=3,
                subset=train,
                distribution="gaussian", verbose=FALSE)
  expect_true(class(o) == "twostage")
})

test_that("twostage predict", {
  d <- data.frame(x, id)
  d2 <- data.frame(x)
  yhat <- predict(o, newdata=d, id="id", n.trees=3)
  ans <- predict(o$o.lmer, newdata=d, allow.new.levels = TRUE) + 
    predict(o$o.gbm, newdata=d2, n.trees=3)
  expect_equal(yhat, ans)
})

test_that("twostage influence drops id", {
  ri <- influence(o)  
  expect_true(names(ri) != "id")
})




