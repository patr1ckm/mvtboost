context("test_plot")

set.seed(123)
n <- 1000
B <- matrix(0,nrow=5,ncol=4)
B[3,1:2] <- 0
B[2,2:3] <- 0
B[4,1] <- 1
B[5,3:4] <- 1

X <- matrix(rbinom(n*(nrow(B)-2),size=1,prob=.5),n,nrow(B)-2)
X2 <- cbind(x1x2=X[,1]*X[,2],x2x3=X[,2]*X[,3])
Xf <- cbind(X,X2)
E <- matrix(rnorm(n*4),nrow=n,ncol=4)
Y <- Xf %*% B + E

# Right now, just makes sure the plots run with default arguments and no errors.
out <- mvtb(Y=Y,X=X,n.trees=100,shrinkage = .5)
x <- rnorm(1000)
y <- x*5 + rnorm(1000)
o1 <- mvtb(Y=y, X=x)


test_that("mvtb.plot", {
  plot(out)
  plot(o1)
})

test_that("mvtb.perspec", {
  mvtb.perspec(out)
  expect_error(mvtb.perspec(o1,1,1))
})

test_that("mvtb.heat", {
  mvtb.heat(mvtb.covex(out, Y=Y,X=X))
  mvtb.heat(mvtb.covex(out, Y=Y,X=X),clust.method = "complete")
  mvtb.heat(mvtb.covex(o1, Y=y, X=x))
  mvtb.heat(t(mvtb.ri(out)))
  mvtb.heat(t(mvtb.ri(out)),clust.method=NULL)
  col <- colorRampPaletteAlpha(RColorBrewer::brewer.pal(9,"Greys"),100)
  mvtb.heat(mvtb.covex(o1, Y=y, X=x),col=col)
})

test_that("return.grid", {
  expect_equal(dim(plot(out,return.grid=T)),c(100,2))
})
