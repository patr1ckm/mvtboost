#' Tune gbm via cross-validation to find the best set of metaparameters
#'
#' This function serves two purposes. First, it computes the cross-validation
#' error across a grid of \code{gbm} metaparameters input by the user, allowing
#' the model to be easily tuned for a given problem. This process can be 
#' executed in parallel on linux-based machines. Second, it alleviates the 
#' burden of selecting a maximum number of trees for a given set of 
#' metaparameters by allowing the algorithm to run until the best
#' number of trees has been selected according to the cross-validation error, in
#' contrast to the standard approach to \code{gbm}, in which a maximum
#' \code{n.trees} must also be tuned. This allows users to avoid two
#' types of problem associated with an inappropriate selection for
#' \code{n.trees}: (1) failing to specify enough trees and therefore using a 
#' sub-optimal model, and (2) specifying far more trees than are necessary, 
#' therefore making the \code{gbm} run for far more time than necessary.
#' 
#' The main output of \code{gbm.cverr} is a data frame with rows corresponding
#' to sets of metaparameters and columns corresponding to (1) the values 
#' defining each set of metaparameters, (2) the minimum cross-validation error
#' corresponding to each row, and (3) the number of trees that yielded the 
#' minimum cross-validation error in each row. These results are intended to 
#' allow users to make informed decisions about the metaparameters passed to 
#' \code{gbm} when fitting the model that will be interpreted and/or used for
#' prediction in the future. 
#' 
#' Note that the metaparamter values passed to \code{w}, \code{var.monotone}, 
#' \code{interaction.depth}, \code{n.minobsinnode}, \code{shrinkage}, and 
#' \code{bag.fraction} will be fully crossed and evaluated.
#'
#' @param x A \eqn{n x p} matrix or data frame of predictors.
#' @param y A \eqn{n x 1} matrix or vector corresponding to the observed 
#' outcome.
#' @param distribution The distribution to use when fitting each \code{gbm}
#' model. For continuous outcomes, the available distributions are "gaussian" 
#' (squared error) and "laplace" (absolute loss). For dichotomous outcomes, the
#' available distributions are  "bernoulli" (logistic regression for 0-1 
#' outcomes) and "adaboost" (the AdaBoost exponential loss for 0-1 outcomes). 
#' Finally, the "poisson" distribution is available for count outcomes.
#' @param cv.folds Number of cross-validation folds to perform.
#' @param fit.best Logical variable indicating whether or not the best set of 
#' metaparameters (estimated according to cross-validation error) will be 
#' utilized to fit and return a \code{\link{gbm.fit}} object to the complete 
#' data.
#' @param nt.start Initial number of trees used to model y.
#' @param nt.inc Number of trees incrementally added until the cross-validation
#' error is minimized or until \code{max.time} is reached (see below).
#' @param verbose If TRUE, then \code{gbm.cverr} will print status information
#' to the console.
#' @param w a vector of weights of the same length as y. NOTE: to
#' evaluate the effect of different weight vectors, a list can be
#' passed to w in which each element follows the structure described
#' above.
#' @param var.monotone an optional vector, the same length as the number of 
#' predictors, indicating which variables have a monotone increasing (+1), 
#' decreasing (-1), or arbitrary (0) relationship with the outcome. NOTE: to
#' evaluate the effect of different monotonicity constraints, a list can be
#' passed to var.monotone in which each element follows the structure described
#' above.
#' @param interaction.depth The maximum depth of variable interactions: 1 builds
#' an additive model, 2 builds a model with up to two-way interactions, etc. 
#' NOTE: Multiple values can be passed in a vector to evaluate the 
#' cross-validation  error using multiple interaction depths.
#' @param n.minobsinnode The minimum number of observations (not total weights) 
#' in the terminal nodes of the trees. NOTE: Multiple values can be passed in a
#' vector to evaluate the cross-validation error using multiple minimum node 
#' sizes.
#' @param shrinkage A shrinkage parameter applied to each tree in the expansion.
#' Also known as the learning rate or step-size reduction. NOTE: Multiple values
#' can be passed in a vector to evaluate the cross-validation error using 
#' multiple shrinkage penalties.
#' @param bag.fraction The fraction of independent training observations (or 
#' patients) randomly selected to propose the next tree in the expansion, 
#' depending on the obs.id vector multiple training data rows may belong to a 
#' single 'patient'. This introduces randomness into the model fit. NOTE: 
#' Multiple values can be passed in a vector to evaluate the cross-validation 
#' error using multiple bag fractions.
#' @param n.cores Number of cores that will be used to estimate cross-validation
#' folds in parallel. Only available on linux-based machines.
#' @param max.time Maximum number of seconds that the model will continue adding
#' trees for a given set of metaparameters. This optional argument allows users 
#' to find the best possible solution in scenarios characterized by limited 
#' computational resources.
#' @param seed Seed that will guarantee \code{gbm.cverr} to produce identical 
#' results across multiple runs. Utilizing \code{set.seed} prior to calling 
#' \code{gbm.cverr} does NOT ensure equal results if \code{bag.fraction < 1}
#'
#' @return An object with 2-5 elements and a summary function. The elements
#' of \code{gbm.cverr.res} are,
#' \item{gbm.fit}{If \code{fit.best} was \code{TRUE}, then this element is the
#' \code{\link{gbm.fit}} object fit to \code{x} and \code{y} using the best set
#' of metaparameters identified by \code{gbm.cverr}.}
#' \item{w}{List of the optional weight vectors provided by the user. Will not
#' be returned if \code{w} was left \code{NULL} when calling \code{gbm.cverr}.}
#' \item{var.montone}{List of the optional monotonoicity parameters proivided
#' by the user. Will not be returned if \code{var.monotone} was left \code{NULL}
#' when calling \code{gbm.cverr}.}
#' \item{cv.err}{A list with length corresponding to the number of metaparameter
#' combinations that were evaluated by \code{gbm.cverr}. Each element is a
#' vector quantifying the cross-validation error across all trees corresponding
#' to the given set of metaparameters.}
#' \item{res}{A data frame with ten columns and as many rows as there were 
#' unique combinations of metaparameters. This data frame is the basis of the
#' summary function for \code{gbm.cverr.res} objects (see below), but it differs
#' from the summary object in two ways: (1) it is not sorted in terms of the
#' minimum cross-validaiton error, but rather according to the order in which
#' the metaparameters were passed to \code{gbm.cverr}, and (2) it contains two
#' additional columns. The column \code{best.meta} is a dummy variable that
#' simply indexes the best set of metaparameters, and \code{timer.end} is a
#' dummy variable indicating whether or not the optimal number of trees was
#' found for a given set of metaparameters (FALSE) or whether the user-specified
#' maximum search time was reached prior to minimizing the cross-validation
#' error (TRUE). If the timer ran out, then the estimated optimal number of 
#' trees is likely underestimated. If this occurred for metaparameter set
#' \code{k}, then in order to evaluate the extent to which the error
#' was still decresaing when the timer ended, we recommend investigating a
#' plot of \code{gbm.cverr.res$cv.err[k]}. The rest of the elements of 
#' \code{res} are discussed below.}
#' 
#' Calling \code{summary(gbm.cverr.res)} produces a data frame with rows 
#' corresponding to sets of metaparameters and columns that denote for each row,
#' \item{min.cv.error}{Minimum cross-validation error resulting from the given
#' set of metaparameters.}
#' \item{w.index}{The index of the (optional) list of weight vectors
#' corresponding to the given set of metaparameters. This will be omitted if
#' a list of weights was not provided to \code{gbm.cverr} through the input
#' parameter \code{w}.}
#' \item{var.monotone.index}{The index of the (optional) list of monotonicity
#' vectors corresponding to the given set of metaparameters. This will be
#' omitted if a list of weights was not provided to \code{gbm.cverr} through
#' the input parameter \code{var.monotone}.}
#' \item{interaction.depth}{The interaction depth corresponding to the
#' given set of metaparameters.}
#' \item{n.minobsinnode}{Minimum number of observations in the terminal
#' nodes of the trees for the given set of metaparameters.}
#' \item{shrinkage}{The shrinkage parameter corresponding to the
#' given set of metaparameters.}
#' \item{bag.fraction}{The fraction of independent training observations 
#' randomly selected to propose the next tree corresponding to
#' the given set of metaparameters.}
#' \item{n.trees}{The optimum number of trees to utilize given the set of
#' metaprameters denoted in the row. Note that entries in this column will be
#' marked with '>=' if the boosting procedure was terminated due to time running
#' out for this set of metaparameters, determined by the user-specified 
#' \code{max.time} passed to \code{gbm.cverr}}
#'
#' In the summary object and output, sets of metaparameters (rows) are ordered 
#' from best (top row) to worst (last row) in terms of the resulting 
#' cross-validation error.
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu)
#'
#' @examples
#'data(wellbeing)
#'y <- wellbeing[,25]
#'x <- wellbeing[,1:20]
#'
#'mm <- gbm.cverr(x = x, y = y, 
#'                distribution = 'gaussian', 
#'                cv.folds = 2, 
#'                
#'                nt.start = 100, 
#'                nt.inc = 100, 
#'                max.time = 1, 
#'                
#'                seed = 12345,
#'                interaction.depth = c(1, 5), 
#'                shrinkage = 0.01,
#'                n.minobsinnode = c(5, 50), 
#'                verbose = TRUE)
#'
#'summary(mm)
#'
#'# Investigate gbm results based on the best set of metaparameters
#'mm$gbm.fit
#'summary(mm$gbm.fit)
#'
#' @export
gbm.cverr <- function(
  # Necessary input
  x, y, distribution = "gaussian", cv.folds = 5, 
  
  # Should the best model be fit to the overall data?
  fit.best = T,
  
  # Find best number of trees
  nt.start = 1000, nt.inc = 1000, verbose = T,
  
  # GBM metaparameters
  w = NULL, # vector or list of vectors
  var.monotone = NULL, # vector or list of vectors
  interaction.depth = 1, # number or vector
  n.minobsinnode = 10, # number or vector
  shrinkage = 0.001, # number or vector
  bag.fraction = 0.5, # number or vector
  
  # Time management
  n.cores = 1, max.time = NULL,
  
  # Reproducible results
  seed = NULL){
  
  ##############################################################################
  ## Set up input and source required packages
  ##############################################################################
  
  # ----------------------------------------------------------------------------
  # Set up grid of metaparameters to evaluate
  # ----------------------------------------------------------------------------
  
  # Weights
  if(!is.list(w)){
    w <- list(w)
  }
  w.inds <- 1:length(w)
  
  # Monotone variance indicators
  return.var <- T
  if(is.null(var.monotone)){return.var <- F}
  if(!is.list(var.monotone)){
    var.monotone <- list(var.monotone)
  }
  var.inds <- 1:length(var.monotone)
  
  # Grid of metaparameters
  meta.grid <- expand.grid(w.index = w.inds,
                           var.monotone.index = var.inds,
                           interaction.depth = interaction.depth,
                           n.minobsinnode = n.minobsinnode,
                           shrinkage = shrinkage,
                           bag.fraction = bag.fraction)
  nmeta <- nrow(meta.grid)
  
  # ----------------------------------------------------------------------------
  # Misc data management
  # ----------------------------------------------------------------------------
  
  x <- as.data.frame(x)
  n <- length(y)
  if(n != nrow(x)){stop('Differing number of observations in x and y')}
  
  # Get seeds for each set of metaparameters
  if(is.null(seed)){seed <- round(stats::runif(1, 0, 1) * .Machine$integer.max)}
  set.seed(seed)
  seeds <- round(stats::runif(nmeta, 0, 1) * .Machine$integer.max)
  
  # Needs weights to evaluate loss function, so use equality if w = NULL
  return.w <- T
  if(is.null(w[[1]])){
    return.w <- F
    w[[1]] <- rep(1, n)
  }
  
  # Check to make sure a valid distribution has been specified
  if(!(distribution %in% 
       c('gaussian', 'adaboost', 'bernoulli', 'laplace', 'poisson'))){
    stop(paste0('Specified distribution unavailable. Please select from:\n',
                'gaussian, adaboost, bernoulli, laplace, poisson'))
  }
  
  if(is.null(max.time)){max.time <- Inf}
  
  ##############################################################################
  ## Program loss functions for each GBM distribution
  ##############################################################################
  loss <- function(fx, yobs, wt, distrb){
    
    err <- NULL
    
    if(distrb == 'gaussian'){
      err <- sum(wt * (yobs - fx)^2) / sum(wt)
    }
    
    if(distrb == 'adaboost'){
      err <- sum(wt * exp(-(2 * yobs - 1) * fx)) / sum(wt)
    }
    
    if(distrb == 'bernoulli'){
      err <- -2 / sum(wt) * sum(wt * (
        yobs * fx - log(1 + exp(fx))
      ))
    }
    
    if(distrb == 'laplace'){
      err <- sum(wt * abs(yobs - fx)) / sum(wt)
    }
    
    if(distrb == 'poisson'){
      err <- -2 * sum(wt * (yobs * fx - exp(fx))) / sum(wt)
    }
    
    return(err)
  }
  
  ##############################################################################
  ## Cross-validate each combination of metaparameters
  ##############################################################################
  res <- lapply(1:nmeta, FUN = function(i){
    if(verbose){
      cat(paste(rep('=', 80), collapse = ''), fill = T)
    }
    # ==========================================================================
    # Step 1: Get current metaparameters
    # ==========================================================================
    w.hold <- w[[meta.grid[i,1]]]
    var.hold <- var.monotone[[meta.grid[i,2]]]
    int.hold <- meta.grid[i,3]
    nmin.hold <- meta.grid[i,4]
    shink.hold <- meta.grid[i,5]
    bag.hold <- meta.grid[i,6]
    
    # ==========================================================================
    # Step 2: Randomize folds
    # ==========================================================================
    set.seed(seeds[i])
    cv.inds <- sample(1:n, replace = F, size = n)
    
    # ==========================================================================
    # Step 3: Fit the model for the current set of metaparameters
    # ==========================================================================
    # Initialize timer
    time.start <- Sys.time()
    
    # --------------------------------------------------------------------------
    # Step 3a: Fit first set of trees
    # --------------------------------------------------------------------------
    # Initialize number of trees
    nt <- nt.start
    if(verbose){
      cat('Fitting trees 1 -', nt, 'using metaparameter set', i,
          'of', nmeta, fill = T)
    }
    
    # Begin cross-validation
    cv.err <- 
      parallel::mclapply(1:cv.folds, mc.cores = n.cores, FUN = function(fld){
        
        # Training and test indices for this fold
        test <- cv.inds[seq(fld, n, by = cv.folds)]
        train <- cv.inds[-seq(fld, n, by = cv.folds)]
        
        # Fit the model to the training data using the initial number of folds
        set.seed(fld)
        mm.cv <- gbm::gbm.fit(x = x[train,], 
                              y = y[train], 
                              distribution = distribution,
                              
                              n.trees = nt,
                              
                              w = w.hold,
                              var.monotone = var.hold,
                              interaction.depth = int.hold,
                              n.minobsinnode = nmin.hold,
                              shrinkage = shink.hold,
                              bag.fraction = bag.hold,
                              
                              nTrain = length(train),
                              keep.data = T,
                              verbose = F)
        
        # Get test error
        ytest.cv <- y[test]
        wtest.cv <- w.hold[test]
        fx.cv <-
          predict(mm.cv, newdata = x[test,], n.trees = 1:nt.start) 
        err <- apply(fx.cv, 2, loss, 
                     yobs = ytest.cv, wt = wtest.cv, distrb = distribution)
        
        return(list(err = err, mm.cv = mm.cv))
        
      })
    
    # Compute CV error across folds
    err <- Reduce('+', lapply(cv.err, FUN = function(err){err[[1]]})) / 
      length(cv.err)
    
    # --------------------------------------------------------------------------
    # Step 3b: Add trees as necessary and as time allows
    # --------------------------------------------------------------------------
    tt <- as.numeric(difftime(Sys.time(), time.start, units = 'secs'))
    keep.going <- tt <= max.time
    while(all((which.min(err) >= (0.90 * length(err))), keep.going)){
      
      # Update status
      if(verbose){
        # cat('Current best tree:', which.min(err), '   ', fill = F)
        cat('Fitting trees', nt + 1, '-', nt + nt.inc, 
            'using metaparameter set', i, 'of', nmeta, fill = T)
      }
      
      # Update number of trees
      nt <- nt + nt.inc
      
      # Update models
      cv.err <- 
        parallel::mclapply(
          1:length(cv.err), mc.cores = n.cores, FUN = function(fld){
            
            test <- cv.inds[seq(fld, n, by = cv.folds)]
            train <- cv.inds[-seq(fld, n, by = cv.folds)]
            
            # Get the current state of the model
            mm.cv <- cv.err[[fld]]$mm.cv
            
            # Add trees
            set.seed(fld + nt)
            mm.cv <- 
              gbm::gbm_more(gbm_fit_obj = mm.cv,
                            num_new_trees = nt.inc, is_verbose = F)
            
            # Get test error for new trees
            ytest.cv <- y[test]
            wtest.cv <- w.hold[test]
            fx.cv <-
              predict(mm.cv, 
                      newdata = x[test,], 
                      n.trees = (length(mm.cv$train.error) - nt.inc + 1):length(
                        mm.cv$train.error), 
                      type = 'link')
            err <- apply(fx.cv, 2, loss, 
                         yobs = ytest.cv, wt = wtest.cv, distrb = distribution)
            
            return(list(err = err, mm.cv = mm.cv))
            
          })
      
      # Update CV error
      err <- c(err,
               Reduce('+', lapply(cv.err, FUN = function(err){err[[1]]})) / 
                 length(cv.err))
      
      # Update timer
      tt <- as.numeric(difftime(Sys.time(), time.start, units = 'secs'))
      keep.going <- tt <= max.time
    }
    
    if(all(verbose, !keep.going)){
      if(verbose){
        cat(paste(rep('-   ', 20), collapse = ''), fill = T)
      }
      cat('Maximum time reached for metaparameter set', i, 'of', 
          nmeta, fill = T)
    }
    
    # --------------------------------------------------------------------------
    # Step 3c: Return CV error for this set of metaparameters
    # --------------------------------------------------------------------------
    return(list(err = err, keep.going = keep.going))
    
  })
  if(verbose){
    cat(paste(rep('=', 80), collapse = ''), fill = T)
  }
  
  # Optimal CV error for each set of metaparameters
  min.cv.error <- unlist(lapply(res, FUN = function(rr){min(rr$err)}))
  min.nt <- unlist(lapply(res, FUN = function(rr){which.min(rr$err)}))
  time.up <- unlist(lapply(res, FUN = function(rr){!rr$keep.going}))
  best.meta <- min.cv.error == min(min.cv.error)
  
  ##############################################################################
  ## Call gbm.fit() on x and y using the best set of metaparameters
  ##############################################################################
  gbm.out <- NULL
  if(fit.best){
    
    w.hold <- w[[meta.grid[best.meta,1]]]
    var.hold <- var.monotone[[meta.grid[best.meta,2]]]
    int.hold <- meta.grid[best.meta,3]
    nmin.hold <- meta.grid[best.meta,4]
    shink.hold <- meta.grid[best.meta,5]
    bag.hold <- meta.grid[best.meta,6]
    
    gbm.out <- gbm::gbm.fit(x = x, 
                            y = y, 
                            distribution = distribution,
                            
                            n.trees = min.nt[best.meta],
                            
                            w = w.hold,
                            var.monotone = var.hold,
                            interaction.depth = int.hold,
                            n.minobsinnode = nmin.hold,
                            shrinkage = shink.hold,
                            bag.fraction = bag.hold,
                            
                            nTrain = n,
                            keep.data = F,
                            verbose = F)
  }
  
  ##############################################################################
  ## Package up results
  ##############################################################################
  
  # Add results to grid of metaparameters
  meta.grid <- cbind(best.meta = best.meta,
                     meta.grid, min.cv.error = min.cv.error, n.trees = min.nt,
                     timer.end = time.up)
  rownames(meta.grid) <- NULL
  
  # Package results
  out <- list(gbm.fit = gbm.out,
              w = w, var.monotone = var.monotone,
              cv.err = lapply(res, FUN = function(rr){rr[[1]]}),
              res = meta.grid)
  
  # Get rid of unnecessary elements (i.e., if people didn't specify weights or 
  # variance monotone parameters)
  if(!fit.best){
    out <- out[-which(names(out) == 'gbm.fit')]
  }
  if(!return.w){
    out <- out[-which(names(out) == 'w')]
    out$res <- out$res[,-which(colnames(out$res) == 'w.index')]
  }
  if(!return.var){
    out <- out[-which(names(out) == 'var.monotone')]
    out$res <- out$res[,-which(colnames(out$res) == 'var.monotone.index')]
  }
  
  # Apply class label
  class(out) <- c('gbm.cverr', class(out))
  
  # Return output
  return(out)
  
}



#' Print gbm.cverr Object
#'
#' \code{print} method for class \code{gbm.cverr}
#'
#' @param x Output from \code{gbm.cverr}
#' @param ... Further arguments passed to or from other methods.
#'
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu)
#'
#'
#' @export
print.gbm.cverr <- function(x, ...){
  out <- x$res
  out <- out[which.min(out$min.cv.error),-c(1, ncol(out))]
  nn.out <- colnames(out)
  
  cat('BEST METAPARAMETERS ', fill = T)
  if('w.index' %in% nn.out){
    cat('   w: input set', out[1, nn.out == 'w.index'], 'of', length(x$w), 
        fill = T)
  }
  if('var.monotone.index' %in% nn.out){
    cat('   var.monotone: input set', out[1, nn.out == 'var.monotone.index'], 
        'of',length(x$var.monotone), fill = T)
  }
  cat('   interaction.depth =', out[1, nn.out == 'interaction.depth'], 
      fill = T)
  cat('   n.minobsinnode =', out[1, nn.out == 'n.minobsinnode'], fill = T)
  cat('   shrinkage =', out[1, nn.out == 'shrinkage'], fill = T)
  cat('   bag.fraction =', out[1, nn.out == 'bag.fraction'], fill = T)
  cat('   n.trees =', out[1, nn.out == 'n.trees'], fill = T)
  cat('RESULTING CROSS-VALIDATION ERROR:', out[1, nn.out == 'min.cv.error'])
  
}


#' Summarizing gbm.cverr Results
#'
#' \code{summary} method for class \code{gbm.cverr}
#'
#' @param object Output from \code{gbm.cverr}
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Calling
#' \code{summary(gbm.cverr)} produces a data frame with rows corresponding to 
#' sets of metaparameters and columns that denote for each row,
#' \item{min.cv.error}{Minimum cross-validation error resulting from the given
#' set of metaparameters.}
#' \item{w.index}{The index of the (optional) list of weight vectors
#' corresponding to the given set of metaparameters. This will be omitted if
#' a list of weights was not provided to \code{gbm.cverr} through the input
#' parameter \code{w}.}
#' \item{var.monotone.index}{The index of the (optional) list of monotonicity
#' vectors corresponding to the given set of metaparameters. This will be
#' omitted if a list of weights was not provided to \code{gbm.cverr} through
#' the input parameter \code{var.monotone}.}
#' \item{interaction.depth}{The interaction depth corresponding to the
#' given set of metaparameters.}
#' \item{n.minobsinnode}{Minimum number of observations in the terminal
#' nodes of the trees for the given set of metaparameters.}
#' \item{shrinkage}{The shrinkage parameter corresponding to the
#' given set of metaparameters.}
#' \item{bag.fraction}{The fraction of independent training observations 
#' randomly selected to propose the next tree corresponding to
#' the given set of metaparameters.}
#' \item{n.trees}{The optimum number of trees to utilize given the set of
#' metaprameters denoted in the row. Note that entries in this column will be
#' marked with '>=' if the boosting procedure was terminated due to time running
#' out for this set of metaparameters, determined by the user-specified 
#' \code{max.time} passed to \code{gbm.cverr}}
#'
#' Sets of metaparameters (rows) are ordered from best (top row) to worst (last 
#' row) in terms of the resulting cross-validation error.
#'
#' @author Daniel B. McArtor (dmcartor@nd.edu)
#'
#' @export
summary.gbm.cverr <- function(object, ...){
  print.res <- object$res[,-1]
  nn <- colnames(print.res)
  print.res <- print.res[,c(which(nn == 'min.cv.error'), 
                            which(nn != 'min.cv.error'))]
  print.res <- out.res <- print.res[order(print.res[,1]),]
  print.res <- print.res[,-ncol(print.res)]
  print.res$n.trees <- paste(print.res$n.trees)

  rownames(print.res) <- rownames(out.res) <- NULL
  # print(print.res, row.names = F)
  # if(any(out.res$timer.end)){
  #   cat('---', fill = T)
  #   cat(
  #     "entries in n.trees starting with '>=' indicate that the maximum allotted",
  #     fill = T)
  #   cat("time expired for that given set of metaparameters, so the true optimum",
  #       fill = T)
  #   cat("number of trees is at least as large as the number provided")
  # }
  class(out.res) <- 'summary.gbmcverr'
  return(out.res)
}

#' @export
print.summary.gbmcverr <- function(x, ...){
  
  nn <- names(x)
  print.res <- data.frame(matrix(unlist(x), ncol = length(x)))
  colnames(print.res) <- nn
  print.res$n.trees <- paste(print.res$n.trees)
  
  for(k in 1:nrow(print.res)){
    if(print.res$timer.end[k] == 1){
      print.res$n.trees[k] <- paste0('>= ', print.res$n.trees[k])
    }
  }
  
  timer <- print.res[,ncol(print.res)]
  print.res <- print.res[,-ncol(print.res)]

  print(print.res, row.names = F)
  if(any(timer == 1)){
    cat('---', fill = T)
    cat(
      "entries in n.trees starting with '>=' indicate that the maximum allotted",
      fill = T)
    cat("time expired for that given set of metaparameters, so the true optimum",
        fill = T)
    cat("number of trees is at least as large as the number provided")
  }
  
}
