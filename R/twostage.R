#' @export
twostage <- function(y, x, id, gbm.first=T, ...){
  x <- as.data.frame(x)
  if(gbm.first){
    d <- data.frame(y, x)
    ob <- gbm::gbm(y ~ . , data=d, distribution="gaussian", ...)
    tr <- suppressWarnings(gbm::gbm.perf(ob, plot.it=F))
    res <- y - predict(ob, newdata = d, n.trees = tr)
    df <- data.frame(y=res, x, id=id)
    form <- as.formula(paste0("y ~ 1 + ", paste0(colnames(x), collapse = " + "),  " + (",
                              paste0(colnames(x), collapse=" + "), " + 1 || id)"))
    ol <- lme4::lmer(form, data=df)
  } else {
    df <- data.frame(y=res, x=x, id=id)
    form <- as.formula(paste0("y ~ 1 + ", paste0(colnames(x), collapse = " + "),  " + (",
                              paste0(colnames(x), collapse=" + "), " + 1 || id)"))
    ol <- lme4::lmer(form, data=df)
    res <- y - fitted(ol)
    df <- data.frame(y=res, x=x, id=id)
    ob <- gbm::gbm(y ~ . -id, data=d, distribution="gaussian", ...)
    tr <- suppressWarnings(gbm::gbm.perf(ob, plot.it=F))
  }
  res <- list(o.lmer=ol, o.gbm=ob, tr=tr)
  class(res) <-  "twostage"
  return(res)
}

#' @export
influence.twostage <- function(object, n.trees = NULL, sort = FALSE, relative = TRUE){
  if(is.null(n.trees)){ n.trees = object$tr}
  inf <- gbm::relative.influence(object$o.gbm, n.trees = n.trees, scale. = FALSE, sort. = sort)
  if(relative){
    inf <- inf / sum(inf) * 100
  }
  return(inf)
}


#' @export
predict.twostage <- function(object, newdata){
  yhat <- predict(object$o.lmer, newdata=newdata, allow.new.levels = TRUE) + 
    predict(object$o.gbm, newdata=newdata, n.trees=object$tr)
  return(yhat)
}