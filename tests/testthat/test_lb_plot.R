
set.seed(104)

ngroups <- 10
group_size <- 50
n <- ngroups * group_size
id <- gl(ngroups, group_size)

train <- sample(n, size = .5*n, replace = FALSE)

x <- rnorm(n)
X <- model.matrix(~x*id)
b <- c(1, .5, rnorm(ncol(X)-2, 0, 1))

ysig <- X %*% b
y <- X %*% b + rnorm(n, 0, .01)
#X <- as.data.frame(x[,2])
Xm <- matrix(X[,2])
coef(lm(y~Xm*id))[1:2]
tol = 1E-6
Xn <- data.frame(matrix(rnorm(n*2), n, 2))
Xt <- data.frame(Xm, Xn, id=id)
dd <- data.frame(y=y, Xm, id)

# marginal plot
par(mfcol=c(1,1), mar=c(4,4,1,3))
plot(y=ysig, x=Xm)
abline(lm(ysig~Xm), col="red")

# plots as a function of id
par(mfcol=c(3, 3), mar=c(3,1,1,1))
for(i in 1:9){plot(y=ysig[id==i], x=Xm[id==i], type="l", main=paste0(i), ylim=c(-2,6))}
par(mfcol=c(1, 1), mar=c(4,4,1,3))
#summary(lme4::lmer(y ~ 1 + Xm + (1+Xm|id)))

pdf("tests/true.pdf")
xyplot(ysig ~ Xm | id, type="l")
dev.off()




# GBM is correct, sort of?
og <- gbm(y ~ ., data = data.frame(y=y, Xt), n.trees=1000, interaction.depth=20,
          distribution="gaussian", cv.folds=3, n.cores=6)

pdf("tests/gbm_test.pdf")
plot(og, i.var=c(1,4), n.trees=gbm.perf(og, plot.it = F))
dev.off()

o <- lmerboost(y = y, X = Xt, id = "id", M = 500, cv.folds = 3, lambda = c(.01, .05, .1), mc.cores=6)
perf.lmerboost(o)

# things we need to test
# for each: including only the variables used for plotting, and all variables
# for each: return.grid=T

# the ids should be given, otherwise the predictions will treat id as unknown (new groups)

# what is the correct answer?

# x1 = factor
# x1 = cont

grid <- plot(o, X=Xt[,c(1,4)], id=2, i.var=c(1), return.grid=T)
plot(o, X=Xt, id=4, i.var=1, return.grid=F)

grid <- plot(o, X=Xt, id=4, i.var=c(1), return.grid=T)
plot(o, X=Xt, id=4, i.var=1, return.grid=F)

# x1 = cont, x2 = cont
# x1 = cont, x2 = factor

Xs <- Xt[id %in% 1:9, ]
Xs$id <- droplevels(Xs$id)

grid <- plot(o, X=Xs, id=4, i.var=c(1, 4), return.grid=T, continuous.resolution = 10)
xyplot(y ~ Var1 | Var2, data=grid, type="l")
plot(o, X=Xs, id=4, i.var=c(1, 4), return.grid=F, continuous.resolution = 10)

grid <- plot(o, X=Xt, id=2, i.var=c(1, 2), return.grid=T)
plot(o, X=Xt, id=2, i.var=c(1, 2), return.grid=F)

# x1 = factor, x2 = cont
# x1 = factor, x2 = factor


