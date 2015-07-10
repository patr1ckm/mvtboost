# mvtboost
Tree boosting for multivariate outcomes in R, built on gbm. Estimates a multivariate additive model of decision trees by iteratively selecting predictors that explain covariance in the outcomes. 

This package can be installed directly from github using the devtools package:

    devtools::install_github("patr1ckm/mvtboost")

## Example usage

    data("mpg",package="ggplot2")
    Y <- mpg[,c("cty","hwy")]      # use both city and highway mileage as dvs
    X <- mpg[,-c(2,8:9)]           # manufacturer, displacement, year, cylinder, transmission, drive, class

    out <- mvtb(Y=Y,X=X,           # data
            n.trees=1000,          # number of trees
            shrinkage=.01,         # shrinkage or learning rate
            interaction.depth=3)   # tree or interaction depth
    ?mvtb
            
## Interpret the model

    summary(out)                   # best trees, relative influences, and covex
    mvtb.ri(out)                   # relative influences
    cluster.covex(out)             # clustered covariance explained in outcomes by predictors
    
    yhat <- predict(out,newdata=X) # predictions
    
    par(mfcol=c(1,2))              # model implied effects for predictor 2 for cty and hwy
    plot(out,response.no=1,predictor.no=2)
    plot(out,response.no=2,predictor.no=2)
    
    heat.covex(out)                # heat map of the clustered covariance explained matrix
    
    mvtb.nonlin(out,X=X,Y=Y)       # indicators of predictors with nonlinear effects

## Tune the model

    out2 <- mvtb(Y=Y,X=X,
            n.trees=1000, 
            shrinkage=.01,
            interaction.depth=3,
            
            bag.frac=.5,          # fit each tree to a sub sample of this fraction
            trainfrac=.5,         # only fit the model to this fraction of the training set
            cv.folds=3,           # number of cross-validation folds
            mc.cores=3)           # run the cross-validation in parallel (not tested on windows)
