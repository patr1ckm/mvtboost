## Resubmission

This improves documentation and user interface, and updates package to 0.4.

- Added `mvtb.covex` that directly returns the covariance explained matrix.
- `mvtb.plot` gains `return.grid` argument to return gbm grid that is used for plotting
- Documentation of functions improved, using more `\code{}` statements, among other things
- Tests for summary functions improved
- Revised package vignettes
- corrected order of variable names for 'age', 'educ' and 'income' in wellbeing data
- `mvtb.heat` gains a `dec` argument to specify the number of digits after the decimal to plot
- example in `?mvtboost` has been updated to produce a cleaner plot

## Test environments

* local OS X 10.10.5 install, R 3.2.2
* Ubuntu 12.04.1 (travis-ci) R 3.2.2
* win-builder (devel and release)

## R CMD check results

Status: OK

R CMD check succeeded

No notes or warnings.