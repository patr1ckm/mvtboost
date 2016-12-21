#' Two stage lmer boost
#' 
#' Fits a random intercept model with the grouping variable \code{id} to the residuals
#' of \code{gbm} fit to the remaining predictors.
#' 
#' @param y vector of observations of length \code{n}
#' @param X matrix or data frame of predictors
#' @param id name or index of grouping variable in x (can have missing values)
#' @param subset subset of observations
#' @param ... arguments passed to gbm. Arguments can be passed as vectors, and 
#' tuning by cross validation will be carried out over expand.grid(...)
#' @export
#' @seealso \link{gbm_grid} 
#' @importFrom lme4 lmer
#' @details The gbm model is cross-validated over a grid of all meta-parameters, 
#' as in \code{gbm_grid}. The model is fit to all predictors, except the grouping 
#' variable. Subsequently, a simple random intercept model (from \code{lmer}) is
#' fit to the residuals, using \code{id} as a grouping variable, with the formula
#' \code{r ~ 1 + (1|id)}. A \code{predict}, and \code{influence} methods are provided
#' for computing predictions and influence of predictors.
#' 
twostage <- function(y, X, id, ..., cv.folds=3, mc.cores=1, subset = NULL){

  x <- as.data.frame(X)
  if(is.null(subset)){ 
    s <- 1:nrow(x)
  } else {
    s <- subset
  }
  if(is.character(id)){
    id <- match(id, colnames(x))
  }
  idname <- colnames(x)[id]
  
  o.grid <- gbm_grid(y=y, x=x[,-id, drop=FALSE], ..., 
                     subset=subset, 
                     cv.folds = cv.folds,
                     mc.cores=mc.cores)
  o.gbm <- o.grid$gbm
  tr <- gbm.perf(o.gbm, plot.it=FALSE, method="cv")
  
  res <- y - predict(o.gbm, newdata=x[,-id, drop=FALSE], n.trees = tr)
  
  df <- data.frame(y=res, x)
  form <- as.formula(paste0("y ~ 1 + (1|", idname, ")"))
  ol <- lme4::lmer(form, data=df, subset = s)
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