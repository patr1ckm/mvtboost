[![Build Status](https://travis-ci.org/patr1ckm/mvtboost.svg?branch=master)](https://travis-ci.org/patr1ckm/mvtboost)
[![codecov.io](https://codecov.io/github/patr1ckm/mvtboost/coverage.svg?branch=master)](https://codecov.io/github/patr1ckm/mvtboost?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mvtboost)](http://cran.r-project.org/package=mvtboost)

# mvtboost

Extends boosted decision trees to multivariate, longitudinal, and hierarchically 
clustered data. Additionally, functions are provided for easy tuning by cross-validated grid search over `n.trees, shrinkage,interaction.depth`, and `n.minobsinnode`.

The package depends on the most recent version of `gbm`, which includes multi-threaded tree-fitting. It can be installed here (eventually deprecated):

    devtools::install_github("patr1ckm/gbm")
    
The package can be installed as follows:

    devtools::install_github("patr1ckm/mvtboost")

Both packages will eventually be pushed to CRAN. 

# mvtb
Tree boosting for multivariate outcomes. Estimates a multivariate additive model of decision trees by iteratively selecting predictors that explain covariance in the outcomes. 

### Example usage

    library(dplyr)
    data("mpg",package="ggplot2")
    Y <- mpg %>% select(cty, hwy) 
    X <- mpg %>% select(-cty, -hwy) %>% 
           mutate_if(is.character, as.factor)

    out <- mvtb(Y=Y,X=X,           # data
            n.trees=1000,          # number of trees
            shrinkage=.01,         # shrinkage or learning rate
            interaction.depth=3)   # tree or interaction depth
    
    
# metb 

Mixed effects tree boosting, useful for longitudinal or hierarchically clustered data. At
each iteration, the terminal node means of each tree are forced to vary by group and shrunk
proportional to group size using `lme4::lmer`. Tuning is done by passing vectors
of meta-parameters as arguments.

### Example usage

    library(dplyr)
    data("mpg",package="ggplot2")
    y <- mpg$cty
    X <- mpg %>% select(-cty, -hwy) %>% 
           mutate_if(is.character, as.factor)
    
    out <- metb(y=y, X=X, id="manufacturer", 
                     n.trees=100,
                     shrinkage=.01, 
                     interaction.depth=3,
                     num_threads=8)
                     
# Grid Tuning by cross-validation

New functions are provided that allow easy grid tuning by cross validation: `gbm.cverr, mvtb_grid`, and `lmerboost`. The grid is defined as `expand.grid(1:cv.folds, ...)` where `...` contains vectors of 
candidate meta-parameter values passed to `n.trees`, `shrinkage`, `interaction.depth`, and `n.minobsinnode`.

With `gbm.cverr`, tuning the number of trees can be carried out by including trees until 1) the cross validation error is minimized or 2) a maximum amount of computation time is reached. This avoids the problem of not including enough trees, or for including more trees than is necessary. 

### Example usage
    
    out <- gbm.cverr(x = X, y = y, 
               distribution = 'gaussian', 
               cv.folds = 2, 
               
               nt.start = 100, 
               nt.inc = 100, 
               max.time = 1, 
               
               seed = 12345,
               interaction.depth = c(1, 5), 
               shrinkage = 0.01,
               n.minobsinnode = c(5, 50), 
               verbose = TRUE)
               
    mm$gbm.fit
    summary(mm$gbm.fit)
        
    
### Limitations

Currently limited to continuous outcomes (generalized outcomes will be added in the future).

The package is experimental, and the interface is subject to change until version
1.0, but usually maintains the original `gbm` interface (2.1.1 and below). 

                  
### Vignettes

    vignette("mvtboost_wellbeing")
    
    
### References

Miller P.J., Lubke G.H, McArtor D.B., Bergeman C.S. (2015) Finding structure in data: A data mining alternative to multivariate multiple regression. Psychological Methods. [arXiv](https://arxiv.org/abs/1511.02025)

Miller P.J., McArtor D.B., Lubke G.H (2017). Abstract: A Gradient Boosting Machine for Hierarchically 
Clustered Data. Multivariate Behavior Research. [arXiv](https://arxiv.org/abs/1702.03994)
