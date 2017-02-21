rm(list=ls())
context("test_gbm_cverr")


################################################################################
## Step 1: Generate data
################################################################################

set.seed(123)
n <- 100
p <- 5
r2 <- 0.50

# Three of ten predictors have an effect: X1^2, X2, and X2*X3
x <- sapply(1:p, FUN = function(k){rnorm(n, 0, 1)})
x.efx <- x[,1:3]
x.efx[,1] <- x.efx[,1]^2
x.efx[,3] <- x[,2] * x[,3]
x.efx <- scale(x.efx)
b <- matrix(c(1,1,1))
xb <- scale(x.efx %*% b) * sqrt(r2)
e <- rnorm(n, 0, sqrt(1 - r2))
y <- scale(c(xb + e))

# Three different versions of Y
y.cont <- y
y.cat <- as.numeric(y <= 0)
y.count <- qpois(p = pnorm(y), lambda = 5)

# Make X a data frame
x <- as.data.frame(x)


################################################################################
## Step 2: Test each loss function with a simple case (one set of meta's)
################################################################################
test_that(desc = 'all loss functions yield results', code = {
  
  # Distributions utilizing a continuous outcome
  mm.gaus <- gbm.cverr(x = x, y = y.cont, 
                       cv.folds = 2, 
                       n.cores = 1, 
                       verbose = F, 
                       distribution = 'gaussian', 
                       interaction.depth = 5,
                       shrinkage = 0.02)
  expect_equal(names(mm.gaus), c('cv.err', 'res'))
  expect_equal(length(mm.gaus$res), 8)
  
  mm.lap <- gbm.cverr(x = x, y = y.cont, 
                      cv.folds = 2, 
                      n.cores = 1, 
                      verbose = F, 
                      distribution = 'laplace', 
                      interaction.depth = 5,
                      shrinkage = 0.02)
  expect_equal(names(mm.lap), c('cv.err', 'res'))
  expect_equal(length(mm.lap$res), 8)
  
  # Distributions using a dichotomous outcome
  mm.bern <- gbm.cverr(x = x, y = y.cat, 
                       cv.folds = 2, 
                       n.cores = 1, 
                       verbose = F, 
                       distribution = 'bernoulli', 
                       interaction.depth = 5,
                       shrinkage = 0.02)
  expect_equal(names(mm.bern), c('cv.err', 'res'))
  expect_equal(length(mm.bern$res), 8)
  
  mm.ada <- gbm.cverr(x = x, y = y.cat, 
                      cv.folds = 2, 
                      n.cores = 1, 
                      verbose = F, 
                      distribution = 'adaboost', 
                      interaction.depth = 5,
                      shrinkage = 0.02)
  expect_equal(names(mm.ada), c('cv.err', 'res'))
  expect_equal(length(mm.ada$res), 8)
  
  # Distributions using count outcome
  mm.pois <- gbm.cverr(x = x, y = y.count, 
                       cv.folds = 2, 
                       n.cores = 1, 
                       verbose = F, 
                       distribution = 'gaussian', 
                       interaction.depth = 5,
                       shrinkage = 0.02)
  expect_equal(names(mm.pois), c('cv.err', 'res'))
  expect_equal(length(mm.pois$res), 8)
  
})



################################################################################
## Step 3: Test that grids of metaparameters work
################################################################################
test_that(desc = 'using grids of metaparameters yield results', code = {
  
  mm.grid <- gbm.cverr(x = x, y = y, 
                       distribution = 'gaussian', 
                       n.cores = 1, cv.folds = 2, 
                       verbose = F, 
                       w = list(rep(1, n), runif(n, 0, 1)), 
                       var.monotone = list(rep(0, p), c(0, 1, rep(0, p-2))), 
                       interaction.depth = c(1, 10), 
                       bag.fraction = c(0.25, 0.5), 
                       n.minobsinnode = c(1, 5), 
                       shrinkage = c(0.01, 0.02))
  
  expect_equal(names(mm.grid), c('w', 'var.monotone', 'cv.err', 'res'))
  expect_equal(length(unlist(mm.grid$w)), 2*n)
  expect_equal(length(unlist(mm.grid$var.monotone)), 2*p)
  expect_equal(length(mm.grid$cv.err), 64)
  expect_equal(dim(mm.grid$res), c(64, 10))
  expect_equal(sum(mm.grid$res$best.meta), 1)
  expect_equal(sum(is.na(mm.grid$res)), 0)
  
})


################################################################################
## Step 4: Test that timing cutoffs and different tree selection strategies work
################################################################################

test_that(desc = 'using time limitations works', code = {
  
  mm.time <- gbm.cverr(x = x, y = y, 
                       distribution = 'gaussian', 
                       n.cores = 1, cv.folds = 2, max.mins = .0001, 
                       nt.start = 100, nt.inc = 100,
                       verbose = F, 
                       w = list(rep(1, n), runif(n, 0, 1)), 
                       var.monotone = list(rep(0, p), c(0, 1, rep(0, p-2))), 
                       interaction.depth = c(1, 10), 
                       bag.fraction = c(0.25, 0.5), 
                       n.minobsinnode = c(1, 5), 
                       shrinkage = c(0.01, 0.02))
  
  expect_equal(names(mm.time), c('w', 'var.monotone', 'cv.err', 'res'))
  expect_equal(length(unlist(mm.time$w)), 2*n)
  expect_equal(length(unlist(mm.time$var.monotone)), 2*p)
  expect_equal(length(mm.time$cv.err), 64)
  expect_equal(dim(mm.time$res), c(64, 10))
  expect_equal(sum(mm.time$res$best.meta), 1)
  expect_equal(sum(is.na(mm.time$res)), 0)
})


################################################################################
## Step 5:  Exact MSE results are returned from a very simple test case (i.e.,
##          test seed setting)
################################################################################

test_that(desc = 'same seed yields same results', code = {
  
  # Distributions utilizing a continuous outcome
  mm.seed <- gbm.cverr(x = x, y = y.cont, seed = 111,
                       cv.folds = 2, 
                       n.cores = 1, 
                       verbose = F, 
                       distribution = 'gaussian', 
                       interaction.depth = 5,
                       shrinkage = 0.02)
  expect_equal(round(mm.seed$res$min.cv.error, 5), 0.90163)
})