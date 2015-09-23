# mvtboost 0.3
- Added wellbeing data and vignette
- Improved documentation for \code{mvtb}, and other functions
- \code{mvtb.summary} gains a covex argument if covex is not wanted
- \code{predict} now drops unnecessary dimensions
- \code{cluster.covex} becomes \code{mvtb.cluster}, a more generic way to cluster any table. obtains a 'plot' argument, a call to \code{mvtb.heat}
- \code{heat.covex} becomes \code{mvtb.heat}, a generic way to heatmap any table (e.g. covex or relative influence). Has better plotting defaults to allow room for y-axis labels. Gets 'cexRow','cexCol' arguments, with sensible defaults. Clustering can also be turned off by setting clust.method=NULL.
- \code{mvtb.perspec} plots variable names by default, gets \code{ylab, xlab, zlab} arguments
- fixed weighting influence bug in \code{mvtb}

# mvtboost 0.2.2
- Corrected biased estimates for covex. Required a small correction factor: 1+(1-shrinkage)
- Added tests to verify correctness of this estimate

# mvtboost 0.2.0

- Added vignette
- First stable release.