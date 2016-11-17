

#' @export
plot.lmerboost <- function(x, X, id, i.var=1,
                           n.trees=min(x$best.trees, na.rm=T), 
                           continuous.resolution=100,
                           return.grid=FALSE, ...){
  
  if(all(is.character(i.var))){
    i.var <- match(i.var, x$xnames)
  }
  grid.levels <- vector("list", length(i.var))
  for(i in i.var){
    if(is.numeric(X[,i])){
      grid.levels[[i]] <- seq(min(X[,i]), max(X[,i]), length.out=continuous.resolution)
    } else {
      grid.levels[[i]] <- levels(X[,i])
    }
  }
  newX <- expand.grid(grid.levels)
  
  if(ncol(X) > length(i.var)){
    # average over other predictors
    X <- data.frame(X, X, id)
  
    for(i in 1:nrow(newX)){
      d <- data.frame(newX[rep(i, nrow(X)), ], X[,-i.var])
      yhat[i] <- predict(x, newdata=d, newid=id)
    }
    loc <- rep(1:nrow(newX), each=nrow(X))
    yhat_m <- tapply(yhat, INDEX = loc, FUN = mean)
    grid <- data.frame(newX, value=yhat_m)
  } else {
    # just compute predictions for grid, no averaging
    colnames(newX) <- colnames(X)[i.var]
    yhat <- predict(x, newdata=newX, newid=id, M=n.trees)$yhat
    grid <- data.frame(newX, y=yhat)
  }

  colnames(newdata) <- colnames(X)  
  if(return.grid){
    return(grid)
  }
 
  if(length(i.var) == 1){
    plot(y=grid$value, x=grid$Var1, type="l", ...)
  } else if(length(i.var) == 2){
    d <- grid
    var.names <- colnames(X)[i.var]
    colnames(d) <- c("X1", "X2", "y")
    f.factor <- sapply(d, is.factor(d))
    
    if(!f.factor[1] && !f.factor[2]){
      print(levelplot(y ~ X1 * X2, data = d, xlab = var.names[i.var[1]], 
                      ylab = var.names[i.var[2]], ...))
    }
    if(!f.factor[1] && f.factor[2]){
      print(xyplot(y ~ X1 | X2, data = d, xlab = var.names[i.var[1]], 
                 ylab = paste("f(", var.names[i.var[1]], ",", 
                              var.names[i.var[2]], ")", sep = ""), type = "l", 
                 panel = panel.xyplot, ...))
    }
    if(f.factor[1] && !f.factor[2]){
      print(xyplot(y ~ X2 | X1, data = d, xlab = var.names[i.var[1]], 
                   ylab = paste("f(", var.names[i.var[2]], ",", 
                                var.names[i.var[1]], ")", sep = ""), type = "l", 
                   panel = panel.xyplot, ...))
    }
    if(f.factor[1] && f.factor[2])
    # all factors
    print(stripplot(X1 ~ y | X2, data = d, xlab = var.names[i.var[2]], 
                    ylab = paste("f(", var.names[i.var[1]], ",", 
                                 var.names[i.var[2]], ")", sep = ""), ...))
  }
}

#' @export
perf.lmerboost <- function(x, threshold = 0, lag = 1, ...){
  M <- length(x$train.err)
  ymax <- c(max(x$test.err, x$train.err, x$oob.err, na.rm = T))
  ymin <- c(min(x$test.err, x$train.err, x$oob.err, na.rm = T))
  
  best.iter <- best_iter(x, threshold = threshold, lag = lag)
  
  plot(x = 1:M, y = x$train.err, type = "l", ylim = c(ymin, ymax), ylab = "error")
  lines(x = 1:M, y = x$test.err, col = "red", lty = 2)
  lines(x = 1:M, y = x$oob.err, col = "blue")
  lines(x = 1:M, y = x$cv.err, col = "red")
  
  abline(v = best.iter)
  
  legend("top", legend = c("train", "test", "oob", "cv"), 
         col = c("black", "red", "blue", "red"),
         lty = c(1, 2, 1, 1), bty = "n")
  x$best.params$err <- formatC(signif(x$best.params$err[[1]], digits=3), digits=3,format="fg", flag="#")
  paramstring <- paste0(names(x$best.params), " = ", x$best.params, collapse = ", ")
  title(sub = paramstring)
}

best_iter <- function(x, threshold, lag, smooth = FALSE){
  err <- x$cv.err
  err <- err[!is.na(err)]
  if(smooth) err <- smooth(err)
  
  best.iter <- which(abs(diff(err, lag = lag)) < threshold)
  
  if(length(best.iter) == 0){
    best.iter <- which.min(x$cv.err)
  } else {
    best.iter <- min(best.iter)
  }
}


