## Resubmission

This improves documentation and user interface, and updates package to 0.3.

- Added wellbeing data and vignette
- Improved documentation for `mvtb`, and other functions
- `mvtb.summary` gains a `covex` argument if covex is not wanted
- `predict` now drops unnecessary dimensions
- `cluster.covex` becomes `mvtb.cluster`, a more generic way to cluster any table. obtains a 'plot' argument, a call to `mvtb.heat`
- `heat.covex` becomes `mvtb.heat`, a generic way to heatmap any table (e.g. covex or relative influence). Has better plotting defaults to allow room for y-axis labels. Gets `cexRow`,`cexCol` arguments, with sensible defaults. Clustering can also be turned off by setting `clust.method=NULL`.
- `mvtb.perspec` plots variable names by default, gets `ylab, xlab, zlab` arguments
- fixed weighting influence bug in `mvtb`
- Corrected biased estimates for covex. Required a small correction factor: 1+(1-shrinkage)
- Added tests to verify correctness of this estimate

## Test environments

* local OS X 10.10.5 install, R 3.2.2
* Ubuntu 12.04.1 (travis-ci) R 3.2.2
* win-builder (devel and release)

## R CMD check results

Status: OK

R CMD check succeeded

No notes or warnings.