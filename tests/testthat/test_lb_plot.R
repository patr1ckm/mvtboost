
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
tol = 1E-6
Xt <- data.frame(X, X, X)

o <- lmerboost(y = y, X = X, id = id, M = 5, cv.folds = 1, lambda = .1)
plot(o, X=data.frame(X, id=id), i.var=1)
