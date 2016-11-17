
#' Marginal (partial dependence) plots for lmerboost objects
#' @param x lmerboost object
#' @param X matrix of predictors
#' @param id index or name of id variable
#' @param i.var index or names of variables to plot over (can include id index)
#' @export
plot.lmerboost <- function(x, X, id, i.var=1,
                           n.trees=min(x$best.trees, na.rm=T), 
                           continuous.resolution=100,
                           return.grid=FALSE, ...){
  
  if(all(is.character(i.var))){
    i.var <- match(i.var, colnames(X))
  }
  if(is.character(id)){
    id <- match(id, colnames(X))
  }
  grid.levels <- vector("list", length(i.var))
  for(i in i.var){
    if(is.numeric(X[,i])){
      grid.levels[[i]] <- seq(min(X[,i]), max(X[,i]), length.out=continuous.resolution)
    } else {
      grid.levels[[i]] <- levels(X[,i])
    }
  }
  
  if(ncol(X) > length(i.var)){
    
    # average over other predictors
    
    # This is very slow. The inner loop here is slow (predict.lmerboost) and it
    # gets bad if newX is very large
    #for(i in 1:nrow(newX)){
    #  d <- data.frame(newX[rep(i, nrow(X)), ], X[,-i.var])
    #  yhat[i] <- mean(predict(x, newdata=d, newid=id)$yhat)
    #}
    grid <- expand.grid(grid.levels[1:2])
    grid.levels[[i+1]] <- 1:nrow(X)
    
    newX <- expand.grid(grid.levels)
    nrows <- prod(lengths(grid.levels[1:2]))
    bigX <- dplyr::bind_rows(lapply(1:nrows, 
                   function(i, i.vars){X[,-i.vars,drop=F]}, i.vars=i.vars))
    newid <- rep(id, nrows)    
    newX <- cbind(newX, bigX)    
    colnames(newX) <- c(colnames(X)[i.vars], colnames(X)[-i.vars])
    
    yhat_long <- predict(x, newdata=newX, newid=newid, M=n.trees)$yhat 
    loc <- rep(1:nrow(X), each=nrows)
    yhat <- tapply(yhat_long, INDEX = loc, FUN = mean)
  } else {
    # just compute predictions for grid, no averaging
    grid <- expand.grid(grid.levels)
    colnames(grid) <- colnames(X)[i.var]
    yhat <- predict(x, newdata=grid, newid=id, M=n.trees)$yhat
  }
  grid$y <- yhat

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


