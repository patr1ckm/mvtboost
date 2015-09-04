## Resubmission

This improves documentation and user interface, and updates package to 0.3.

- Improved documentation for \code{mvtb}, and other functions
- \code{predict} now drops unecessary dimensions
- \code{cluster.covex} obtains a 'plot' argument, a call to \code{heat.covex}
- \code{heat.covex} has better plotting defaults to allow room for y-axis labels. Gets a 'mar' argument, with a sensible default
- \code{plot.perspec} plots variable names by default, gets \code{ylab, xlab, zlab} arguments

## Test environments

* local OS X install, R 3.2.2
* Ubuntu 12.04.1 (travis-ci) R 3.2.2
* win-builder (devel and release)

## R CMD check results

Status: OK

R CMD check succeeded

No notes or warnings.