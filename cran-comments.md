## Resubmission
This is a resubmission. It addresses errors on CRAN arising from namespace issues, and adds a vignette.

* Revised NAMESPACE to correctly import functions from stats, utils, graphics, and grDevices
* Revised DESCRIPTION so that stats, utils, graphics and grDevices are in 'Imports' rather than 'Depends'
* Revised DESCRIPTION so that knitr, rmarkdown, and ggplot2 can be used in the vignette

## Test environments

* local OS X install, R 3.2.1
* Ubuntu 12.04.1 (travis-ci) R 3.2.1
* win-builder (devel and release)

## R CMD check results

* There was 1 NOTE:

Maintainer: ‘Patrick Miller <patrick.mil10@gmail.com>’
New submission

