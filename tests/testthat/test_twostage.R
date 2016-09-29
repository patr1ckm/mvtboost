

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
X <- as.data.frame(x)

o <- twostage(y=y, x=X, id=id, n.trees=1000, shrinkage=.005)
d <- data.frame(y, x, id)
yhat <- predict(o, newdata=d)
mse <- mean((y - yhat)^2)

ri <- influence(o)

o <- twostage(y = y, x = X, id = id, subset = train,  n.trees=1000, shrinkage=.005)
yhat <- predict(o, newdata = d)


