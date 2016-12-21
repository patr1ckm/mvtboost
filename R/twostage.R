#' Two stage lmer boost
#' @param y outcome
#' @param x predictor
#' @param id grouping variable
#' @param subset subset of observations
#' @param ... arguments passed to gbm
#' @export
twostage <- function(y, x, id, ..., cv.folds=3, mc.cores=1, subset = NULL){

  x <- as.data.frame(x)
  if(is.null(subset)){ 
    s <- 1:nrow(x)
  } else {
    s <- subset
  }
  d <- data.frame(y, x)
  o.grid <- gbm_grid(y=y, x=x, ..., subset=subset, cv.folds = cv.folds,
                     mc.cores=mc.cores)
  o.gbm <- o.grid$gbm
  tr <- gbm.perf(o.gbm, plot.it=FALSE, method="cv")
  res <- y - predict(o.gbm, newdata = x, n.trees = tr)
  df <- data.frame(y=res, x, id=id)
  ol <- lme4::lmer(y ~ 1 + (1 | id), data=df, subset = s)
  res <- list(o.lmer=ol, o.gbm=o.gbm, tr=tr)
  class(res) <-  "twostage"
  return(res)
}


#' Prediction for twostage objects
#' @param object object of class \code{twostage}
#' @param newdata data frame of predictors augmented with grouping variable
#' @param id character or index of grouping variable in \code{newdata}
#' @param n.trees number of trees to use. If NULL, uses the best number of trees in \code{object}
#' @param ... unused
#' @export
predict.twostage <- function(object, newdata, id, n.trees = NULL, ...){
  if(is.null(n.trees)){
    n.trees <- object$tr
  }
  if(is.character(id)){
    id <- match(id, colnames(newdata))
  }
  yhat <- predict(object$o.lmer, newdata=newdata, allow.new.levels = TRUE) + 
    predict(object$o.gbm, newdata=newdata[,-id, drop=F], n.trees=n.trees)
  return(yhat)
}