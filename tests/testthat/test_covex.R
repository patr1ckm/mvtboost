k <- 3
#set.seed(100)

# no error, linear
q <- 3
p <- 1
n <- 1000
B <- matrix(0,p,q)
B[1,c(1,3)] <- 1
X <- matrix(rnorm(n),n,1)

X <- ifelse(X < .5,0,1) 
Y <- X %*% B

a <- cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)
for(i in seq_along(shrink)){
  out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=100)
  expect_true(all((out$covex[out$covex > 0] - a) < 1E-12))
}

X <- matrix(rnorm(n),n,1)
Y <- X %*% B

a <- cov(Y)[1,3]
shrink <- c(.2,.5,.9,1)
for(i in seq_along(shrink)){
  out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=100)
  expect_true(all((out$covex[out$covex > 0] - a) < 1E-12))
}
out <- mvtb(Y=Y,X=X,shrinkage=1,n.trees=100)


#colnames(Y) <- paste0("Y",1:k)
#colnames(d$X) <- paste0("X",1:ncol(d$X))

# boostCov <- function(out){
#   cd <- out$covex
#   k  <- length(out$ynames)
#   covmat <- array(0,dim=c(k,k,ncol(cd)))
#   a <- diag(k)
#   for(i in 1:ncol(cd)){
#     a[lower.tri(a,diag=TRUE)] <- cd[,i]
#     a <- a + t(a) - diag(diag(a))
#     covmat[,,i] <- a    
#   }
#   return(covmat)
# }
# out <- mvtb(Y=Y,X=X,shrinkage=1,n.trees=100)
# out <- mvtb(Y=Y,X=X,shrinkage=1,n.trees=1)
# out <- mvtb(Y=Y,X=X,shrinkage=1,n.trees=2)
# out <- mvtb(Y=Y,X=X,shrinkage=1,n.trees=3)
# 
# out <- mvtb(Y=Y,X=X,shrinkage=.5,n.trees=10)
# (tr <- cov(Y)[,1]-out$covex[1:3])
# out <- mvtb(Y=Y,X=X,shrinkage=.5,n.trees=2)
# (tr <- cov(Y)[1,3]-out$covex[1:3])
# out <- mvtb(Y=Y,X=X,shrinkage=.5,n.trees=3)
# 
# bias2 <- matrix(0,100,3)
# for(i in 1:100) {
#   out <- mvtb(Y=Y,X=X,shrinkage=.5,n.trees=i)
#   bias2[i,] <- cov(Y)[,1]-out$covex[1:3]
# }
# 
# bias1 <- matrix(0,100,3)
# for(i in 1:100) {
#   out <- mvtb(Y=Y,X=X,shrinkage=1,n.trees=i)
#   bias1[i,] <- cov(Y)[,1]-out$covex[1:3]
# }
# 
# shrink <- seq(.5,1,by=.01)
# bias <- rep(0,length(shrink))
# for(i in seq_along(shrink)) {
#   out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=1000)
#   bias[i] <- cov(Y)[1,3]-out$covex[3]
# }
# summary(o <- lm(bias~I(shrink^2)))
# var(X)/coef(o)
# 
# q <- 3
# p <- 1
# n <- 1000
# B <- matrix(0,p,q)
# B[1,c(1,3)] <- 1
# 
# 
# simv <- function(v){
#   X <- matrix(rnorm(n,0,v),n,1)
#   Y <- X %*% B
#   shrink <- seq(.5,1,by=.1)
#   bias <- rep(0,length(shrink))
#   for(i in seq_along(shrink)) {
#     out <- mvtb(Y=Y,X=X,shrinkage=shrink[i],n.trees=1000)
#     bias[i] <- cov(Y)[1,3]-out$covex[3]
#   }
#   return(bias)
# }
# simv.res1 <- mclapply(c(1:5),simv,mc.cores=5)
# 
# vars <-rep(1:5,each=6)
# shrink <- rep(seq(.5,1,by=.1),times=5)
# summary(lm(unlist(simv.res1)~I(vars^2)*I(shrink^2)))
# 
# lvec <- function(a,diag=TRUE){a[lower.tri(a,diag=diag)]}
# 
# B <- matrix(0,p,2)
# B[1,] <- 1
# v <- .25
# k <- 2
# 
# set.seed(101)
# X <- matrix(rnorm(n,0,2),n,1)
# Xb <- ifelse(X < .5,0,1)
# shrink <- .25
# Y <- Xb %*% B
# Ys <- scale(Y,center = T,scale = F)
# 
# o1 <- mvtb(Y=Ys,X=Xb,n.trees=1,shrinkage=.25,s=1:1000,interaction.depth = 1)
# o2 <- mvtb(Y=Ys,X=Xb,n.trees=2,shrinkage=.25,s=1:1000,interaction.depth = 1)
# o3 <- mvtb(Y=Ys,X=Xb,n.trees=100,shrinkage=.25,s=1:1000,interaction.depth = 1)
# 
# (covex1 <- o1$covex)
# (covex2 <- o2$covex)
# 
# # iter 1
# y1 <- Ys[,1]
# y2 <- Ys[,2]
# e1 <- y1-predict(lm(y1~Xb))*v
# e2 <- y2-predict(lm(y2~Xb))*v
# sx <- var(Xb)
# var(e1) - sx*(1-v)^2
# cov(e1,y2) - sx*(1-v)
# D1 <- cbind(e1,y2)
# D2 <- cbind(y1,e2)
# cov(Y)-cov(D1)
# cov(Y)-cov(D2)
# lvec((cov(Y)-cov(D1)) + (cov(Y) - cov(D2)))-covex1
# covex <- c(sx*(1-(1-v)^2),sx*v*k,sx*(1-(1-v)^2))
# covex-covex1
# # iter 2
# D <- cbind(e1,e2)
# all(abs(D-o1$resid) < 1E-8)
# 
# e11 <- e1 - predict(lm(e1~Xb))*v
# e22 <- e2 - predict(lm(e2~Xb))*v
# cov(e11,e2)
# D1 <- cbind(e11,e2)
# D2 <- cbind(e1,e22)
# cov(D)-cov(D1)
# cov(D)-cov(D2)
# covex + lvec((cov(D)-cov(D1)) + (cov(D) - cov(D2))) - o2$covex
# covex2 <- covex + lvec((cov(D)-cov(D1)) + (cov(D) - cov(D2)))
# sx*(v-v^2)
# Xm <- cbind(1,Xb)
# H <- Xm %*% solve(t(Xm) %*% Xm) %*% t(Xm)
# t(y1) %*% (I - H*v) %*% y1 / (n-1)
# var(e11)
# 
# 
# library(expm)
# (I - H*v) %^% 2

