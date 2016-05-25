## Resubmission

- `mvtb.covex` now estimates the covariance explained matrix, it is not estimated by default 
- Arguments in `mvtb` correspond directly to `gbm`
- Tests are more flexible, maintainable, have improved coverage
- Documentation for several functions has been simplified
- `mvtb.fit` is now exported
- Any `distribution` can be any supported by `mvtb`, but covariance explained can only
be estimated when `distribution=gaussian`

Updates package to 0.5.0


## Test environments

* local OS X 10.11.2 install, R 3.3.0
* Ubuntu 12.04.5 (travis-ci) R 3.3.0
* win-builder (devel and release)

## R CMD check results

Status: OK

R CMD check results
0 errors | 0 warnings | 0 notes

R CMD check succeeded

