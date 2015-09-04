# mvtboost 0.3
- Improved documentation for \code{mvtb}, and other functions
- \code{predict} now drops unecessary dimensions
- \code{cluster.covex} obtains a 'plot' argument, a call to \code{heat.covex}
- \code{heat.covex} has better plotting defaults to allow room for y-axis labels. Gets a 'mar' argument, with a sensible default
- \code{plot.perspec} plots variable names by default, gets \code{ylab, xlab, zlab} arguments


# mvtboost 0.2.2
- Corrected biased estimates for covex. Required a small correction factor: 1+(1-shrinkage)
- Added tests to verify correctness of this estimate

# mvtboost 0.2.0

- Added vignette
- First stable release.