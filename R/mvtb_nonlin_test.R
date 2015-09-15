

## n <- 1000
## set.seed(n)
## x <- rnorm(n)
## x2 <- rnorm(n)
## x3 <- rnorm(n)
## x4 <- rnorm(n)
## e <- rnorm(n)
## y <- x*x2-x3*x4 + e
## y2 <- x3*x2-x4*x+e
## X <- data.frame(x=x,x2=x2,x3=x3,x4=x4)
## Y <- data.frame(y=y,y2=y2)
## d <- data.frame(y=y,y2=y2,x=x,x2=x2,x3=x3,x4=x4)

## out <- mvtb(X=X,Y=Y,niter=1000,shrinkage=.01,interaction.depth=5,bag.frac=1,s=1:nrow(d))

## res <- mvtb.interactions(out=out,X=X,Y=Y,n.trees=out2$gbm.call$best.trees)
## res[[1]]

#' Detect departures from linearity from a multivariate tree boosting model.
#' @param object object of class \code{mvtb}
#' @param X matrix of predictors
#' @param Y matrix of responses
#' @param n.trees number of trees. Defaults to the minimum number of trees given that minimize CV, test, training error.
#' @param detect method for testing possible non-linear effects or interactions. Possible values are "grid", "influence", and "lm". See details.
#' @param scale For method "influence", whether the resulting influences are scaled to sum to 100.
#' @return For each outcome, a list is produced showing the interactions in two forms. The first is $rank.list, which shows the nonlinear effect for each pair of predictors ranked according to the size of the departure from non-linearity. 
#' The second, $interactions, shows the departure from non-linearity for all pairs of predictors.
#' @details This function provides statistics means to detect departures from linearity in the multivariate boosting model for any outcome as a function of pairs of predictors. 
#' These could be interactions between pairs of variables, or more general non-linear effects. Please note that these methods should be interpreted as exploratory only.
#' 
#' Several methods are provided for detecting departures from non-linearity in pairs from pairs of predictors. 
#' The "grid" method computes a grid of the model implied predictions as a function of two predictors, averaging over the others. A linear model predicting the observed outcomes from the predicted values is fit, and the mean squared residuals (times 1000) are reported. Large residuals indicate deviations from linearity.
#' 
#' The "influence" method computes the reductions of SSE attributable to predictors after the first split on the tree. These reductions in sums of squared error (or influnences) indicate to what extent individual predictors capture deviations from linear, main effects.
#' 
#' The "lm" method is the same as the "grid" method, but produces the grid of predicted values by conditioning on the average values of the other predictors rather than averaging over the values of the other predictors (see Elith et al., 2008) . Like the "grid" approach, large residuals from a linear model (times 1000) indiciate departures from linearity. 
#'
#' These methods are not necessarily overlapping, and can produce different results. We suggest using several approaches, and follow up by plotting the model implied effects of the two predictors.
#' The gbm package contains the function interact.gbm to detect interactions. See ?interact.gbm for details of this function, which can be used directly on individual mvtb output models.
#' @seealso \code{interact.gbm}, \code{mvtb.perspec}, \code{plot.gbm}
#' @references 
#' Miller P.J., Lubke G.H, McArtor D.B., Bergeman C.S. (Submitted) Finding structure in data: A data mining alternative to multivariate multiple regression. Psychological Methods.
#' 
#' Elith, J., Leathwick, J. R., & Hastie, T. (2008). A working guide to boosted regression trees. Journal of Animal Ecology, 77(4), 802-813.
#' 
#' Friedman, J. H., & Meulman, J. J. (2003). Multiple additive regression trees with application in epidemiology. Statistics in medicine, 22(9), 1365-1381.
#' @export
mvtb.nonlin <-function(object, X, Y, n.trees=NULL,detect="grid",scale=TRUE) {
  #
  # p. miller, February 2015. Updated for multiple outcome variables
  # j. leathwick, j. elith - May 2007
  out <- object
  if(any(unlist(lapply(out,function(li){is.raw(li)})))){
    out <- mvtb.uncomp(out)
  }
  if(is.null(n.trees)) { n.trees <- min(unlist(out$best.trees)) }
  data <- X
  n.preds <- ncol(data)
  if(!is.null(colnames(data))) { 
    pred.names <- colnames(data)
  } else {
    pred.names <- out$models[[1]]$var.names
  }
  if(!is.null(colnames(Y))){
    col.names <- colnames(Y)
  } else {
    col.names <- out$ynames
  }
  Y <- as.matrix(Y)
  
  if(detect=="grid") {
    detect.function <- 1
  } else if (detect == "influence") {
    detect.function <- 3
  } else {
    detect.function <- 2
  }
  
  doone <- function(which.y,mvtb.out,detect.function=1,data=data,n.preds=n.preds,pred.names=pred.names,n.trees=n.trees,scale=scale) {
    if(detect.function==1) {
      cross.tab <- intx.grid(mvtb.out,num.pred=n.preds,k=which.y, n.trees=n.trees)
    } else if (detect.function==2) {
      cross.tab <- intx.lm(mvtb.out,n.trees=n.trees,which.y=which.y,data=data,n.preds=n.preds,pred.names=pred.names)
    } else {
      cross.tab <- intx.influence(mvtb.out,k=which.y,n.trees=n.trees,scale=scale)
    }
    dimnames(cross.tab) <- list(pred.names,pred.names)
    
    ## create an index of the values in descending order
    
    search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")
    
    n.important <- max(2,round(0.1 * ((n.preds^2)/2),0))
    var1.names <- rep(" ",n.important)
    var1.index <- rep(0,n.important)
    var2.names <- rep(" ",n.important)
    var2.index <- rep(0,n.important)
    int.size <- rep(0,n.important)
    
    for (i in 1:n.important) {
      index.match <- match(i,search.index)
      j <- trunc(index.match/n.preds) + 1
      var1.index[i] <- j
      var1.names[i] <- pred.names[j]
      k <- index.match%%n.preds
      if (k > 0) {   #only do this if k > 0 - otherwise we have all zeros from here on 
        var2.index[i] <- k
        var2.names[i] <- pred.names[k]
        int.size[i] <- cross.tab[k,j]
      }
    }
    
    rank.list <- data.frame(var1.index,var1.names,var2.index,var2.names,nonlin.size=int.size)
    
    return(list(rank.list = rank.list, nonlin.full = cross.tab))
  }
  
  res <- lapply(1:ncol(Y),doone,mvtb.out=out,detect.function=detect.function,n.trees=n.trees,pred.names=pred.names,n.preds=n.preds,data=data,scale=scale)
  names(res) <- colnames(Y)
  return(res)
}

#' @importFrom stats residuals resid lm
intx.grid <- function(mvtb.out,num.pred,k=1,n.trees) {
  #gbm.obj <- convert.mvtb.gbm(r=mvtb.out,k=k)
  gbm.obj <- mvtb.out$models[[k]]
  cross.tab <- matrix(0,num.pred,num.pred)
  #dimnames(cross.tab) <- list(pred.names,pred.names)
  for(i in 1:(num.pred-1)) {
    for(j in (i+1):num.pred) {
      grid <- gbm::plot.gbm(gbm.obj,i.var=c(i,j),n.trees=n.trees,return.grid=TRUE)
      cross.tab[i,j] <- mean(residuals(lm(y~.,data=grid))^2)*1000
      #fi <- rep(tapply(grid$y,list(factor(grid[,1])),mean),times=2)
      #fj <- rep(tapply(grid$y,list(factor(grid[,2])),mean),each=2)
      
      #stat <- sum((Fij - fi - fj + mean(Fij))^2)  
      #pall <- nrow(grid)
      #p1 <- length(unique(grid[,1]))
      #p2 <- length(unique(grid[,2]))
      #Dbar.2D <- matrix(1/pall,pall,pall)
      #Di. <- matrix(1,p1,p1)%x%diag(p2)/p2
      #D.j <- diag(p2)%x%matrix(1,p1,p1)/p1
      #D <- diag(pall) - Di. - D.j + Dbar.2D
      #X <- grid[,1:2]
      
      #cross.tab[i,j] <- t(D %*% Fij) %*% (D %*% Fij)  
    }
  }
  return(cross.tab)
}

#' @importFrom stats residuals resid lm
intx.lm <- function (out,n.trees,which.y,data,n.preds,pred.names) {
  cross.tab <- matrix(0,n.preds,n.preds)
  for (i in 1:(n.preds - 1)) {  # step through the predictor set
    if (is.vector(data[,i])) {  # create a sequence through the range
      x.var <- seq(min(data[,i],na.rm=T),max(data[,i],na.rm=T),length = 20)
    }
    else {                      # otherwise set up simple factor variable
      x.var <- factor(names(table(data[,i])),levels = levels(data[,i]))
    }
    x.length <- length(x.var)
    for (j in (i+1):n.preds) { #create vector or factor data for second variable      
      if (is.vector(data[,j])) {
        y.var <- seq(min(data[,j],na.rm=T),max(data[,j],na.rm=T),length = 20)
      }
      else {
        y.var <- factor(names(table(data[,j])),levels = levels(data[,j]))
      }
      y.length <- length(y.var)
      pred.frame <- expand.grid(list(x.var,y.var))
      names(pred.frame) <- c(pred.names[i],pred.names[j])
      n <- 3 # and add the balance of the variables to it
      for (k in 1:n.preds) {
        if (k != i & k != j) {
          if (is.vector(data[,k])) {  # either with the mean
            pred.frame[,n] <- mean(data[,k],na.rm=T)
          }
          else {   # or the most common factor level
            temp.table <- sort(table(data[,k]),decreasing = TRUE)
            pred.frame[,n] <- rep(names(temp.table)[1],x.length * y.length)
            pred.frame[,n] <- as.factor(pred.frame[,n])
          }
          names(pred.frame)[n] <- pred.names[k]
          n <- n + 1
        }
      }        
      ## form the prediction
      prediction <- predict.mvtb(out,newdata=data.frame(pred.frame),n.trees = n.trees,drop=FALSE)[,which.y,]
      interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[,2]))             
      interaction.flag <- round(mean(resid(interaction.test.model)^2)*1000,2)
      cross.tab[i,j] <- interaction.flag
    }   # end of j loop
  }  # end of i loop
  return(cross.tab)
}


## Purpose: Another way to detect possible interactions. This is done by computing
## the relative influence attributable to splitting on any predictors after the first split.
## The results are printed for each outcome variable in a table where the first split is on the 
## predictor in the column, and the other splits are in the rows.

## Arguments: 
##  mvtb.obj - object from mvtb
##  scale    - influences are scaled to sum to 100

## Value:
##   list of interaction tables for each outcome variable. In each table, we list the 
##   reduction in sums of squared errors attributable to splitting on variables in the row
##   AFTER the first split on the variable in the column. 

intx.influence <- function(object,k=1,n.trees,scale=TRUE) {
  #do.one <- function(k,mvtb.obj,scale) {
  #gbm.object <- convert.mvtb.gbm(mvtb.obj,k)
  mvtb.obj <- object
  gbm.object <- mvtb.obj$models[[k]]
  trees <- gbm.object$trees
  pred.names <- gbm.object$var.names
  get.intx.sse <- function(t) {
    intx_sse <- lapply(split(t[[6]][-1], t[[1]][-1]), sum)
  }
  intx.tab <- matrix(0,nrow=length(pred.names),ncol=length(pred.names))
  rownames(intx.tab) <- colnames(intx.tab) <- pred.names
  for(m in 1:n.trees) {
    var <- trees[[m]][[1]][1] + 1
    intx.sse <- get.intx.sse(trees[[m]])
    intx.sse <- intx.sse[names(intx.sse) != "-1"]
    int.vars <- as.numeric(names(intx.sse))+1
    intx.tab[var,int.vars] <- intx.tab[var,int.vars] + unlist(intx.sse)
  }
  if(scale) {
    intx.tab <- intx.tab/sum(intx.tab)*100
  }
  return(intx.tab)
  #}
  #res <-  lapply(1:length(mvtb.obj$finaltree),do.one,mvtb.obj,scale)
  #names(res) <- colnames(mvtb.obj$ri[[1]])
  #return(res)
}

## GBM Tree structure
## [[1]] - INDEX: vector of indeces of splitting variables. -1 indicates a terminal node. starts from 0.
## [[2]] - prediction: split point, or c.split describing the split
## [[3]] - Node assignments (left). 
## [[4]] - Node assignments (right). 
## [[5]] - Missing node.
## These describe where the where the L and R and M nodes are in the table. For the first round, implementation, this won't be used.
## [[6]] - SSE reduction!
## [[7]] - Total observations in node (if weights = 1)