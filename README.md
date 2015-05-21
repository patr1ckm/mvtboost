# mvtboost
Tree boosting for multivariate outcomes in R, built on gbm. Estimates a multivariate additive model of decision trees by iteratively selecting predictors that explain covariance in the outcomes. 

This package can be installed directly from github using the devtools package:

    devtools::install_github("patr1ckm/mvtboost")

Example usage

    ?mvtb
    mvtb(Y=Y,X=X, ... )
