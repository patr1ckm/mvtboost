#' Two stage lmer boost
#' @param y outcome
#' @param x predictor
#' @param id grouping variable
#' @param subset subset of observations
#' @param gbm.first whether gbm is fit first or not
#' @param ... arguments passed to gbm
#' @export
twostage <- function(y, x, id, subset = NULL, ...){

  x <- as.data.frame(x)
  if(is.null(subset)){ 
    s <- 1:nrow(x)
  } else {
    s <- subset
  }
  d <- data.frame(y, x)
  ob <- gbm::gbm(y ~ . , data=d[s, ], distribution="gaussian", ...)
  tr <- suppressWarnings(gbm::gbm.perf(ob, plot.it=F))
  res <- y - predict(ob, newdata = d, n.trees = tr)
  df <- data.frame(y=res, x, id=id)
  form <- stats::as.formula(paste0("y ~ 1 + ", paste0(colnames(x), collapse = " + "),  " + (",
                            paste0(colnames(x), collapse=" + "), " + 1 || id)"))
  ol <- lme4::lmer(form, data=df, subset = s)
  res <- list(o.lmer=ol, o.gbm=ob, tr=tr)
  class(res) <-  "twostage"
  return(res)
}


#' @export
predict.twostage <- function(object, newdata, n.trees = NULL, ...){
  if(is.null(n.trees)){
    n.trees <- object$tr
  }
  yhat <- predict(object$o.lmer, newdata=newdata, allow.new.levels = TRUE) + 
    predict(object$o.gbm, newdata=newdata, n.trees=n.trees)
  return(yhat)
}