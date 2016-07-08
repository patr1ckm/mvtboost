
context("speed")
x <- matrix(rnorm(1000*5), 1000, 5)
y <- x * 5 + rnorm(1000)

o <- replicate(5, system.time(mvtb.sep(Y=y, X=x, n.trees=100))[3])
o2 <- replicate(5, system.time(mvtb(Y=y, X=x, n.trees=100))[3])
expect_lt(mean(o), mean(o2))
cat(paste0("percent improvement: ", round((mean(o2) - mean(o))/mean(o2), 3)*100, "%"), fill=T)
