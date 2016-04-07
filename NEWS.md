# mvtboost 0.4.4
- `mvtb.covex` now estimates the covariance explained matrix. 
- Arguments in `mvtb` correspond directly to `gbm`
- Better testing
- Revised documentation
- `mvtb.fit` is now exported
- `distribution` can be any supported by `gbm`, but covariance explained can only
be estimated when `distribution=gaussian`

# mvtboost 0.4
- Added `mvtb.covex` that directly returns the covariance explained matrix.
- `mvtb.plot` gains `return.grid` argument to return gbm grid that is used for plotting
- Documentation of functions improved, using more `\code{}` statements, among other things
- Tests for summary functions improved
- Revised package vignettes
- corrected order of variable names for 'age', 'educ' and 'income' in wellbeing data
- `mvtb.heat` gains a `dec` argument to specify the number of digits after the decimal to plot
- example in `?mvtboost` has been updated to produce a cleaner plot

# mvtboost 0.3
- Added wellbeing data and vignette
- Improved documentation for `mvtb`, and other functions
- `mvtb.summary` gains a `covex` argument if covex is not wanted
- `predict` now drops unnecessary dimensions
- `cluster.covex` becomes `mvtb.cluster`, a more generic way to cluster any table. obtains a 'plot' argument, a call to `mvtb.heat`
- `heat.covex` becomes `mvtb.heat`, a generic way to heatmap any table (e.g. covex or relative influence). Has better plotting defaults to allow room for y-axis labels. Gets `cexRow`,`cexCol` arguments, with sensible defaults. Clustering can also be turned off by setting `clust.method=NULL`.
- `mvtb.perspec` plots variable names by default, gets `ylab, xlab, zlab` arguments
- fixed weighting influence bug in `mvtb`

# mvtboost 0.2.2
- Corrected biased estimates for covex. Required a small correction factor: 1+(1-shrinkage)
- Added tests to verify correctness of this estimate

# mvtboost 0.2.0

- Added vignette
- First stable release.