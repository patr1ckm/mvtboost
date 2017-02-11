
set.seed(104)

ngroups <- 10
group_size <- 50
n <- ngroups * group_size
id <- gl(ngroups, group_size)

train <- sample(n, size = .5*n, replace = FALSE)

x <- rnorm(n)
X <- model.matrix(~x*id)
b <- c(0, 1, 
       rep(0, 9),
       c(1, 1, 1, 1, -2, -2, -2, -2, -2))

ysig <- X %*% b
y <- X %*% b + rnorm(n, 0, .01)
Xm <- matrix(X[,2])
coef(lm(y~Xm*id))[1:2]
tol = 1E-6
Xn <- data.frame(matrix(rnorm(n*2), n, 2))
Xt <- data.frame(Xm, Xn, id=id)
dd <- data.frame(y=y, Xm, id)



