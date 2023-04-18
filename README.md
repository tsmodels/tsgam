
[![R-CMD-check](https://github.com/tsmodels/tsgam/workflows/R-CMD-check/badge.svg)](https://github.com/tsmodels/tsgam/actions)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--04--18-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-0.3.0-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsgam)](https://cran.r-project.org/package=tsgam)

# tsgam

Wrapper for the Generalized Additive Model in package mgcv to conform
with tsmodels output and framework. It is a more abstract/general layer
than other models in tsmodels due to the nature of GAMs.

## Installation

The package can be installed from the tsmodels github repo.

``` r
remotes::install_github("tsmodels/tsgam", dependencies = TRUE)
```

A short vignette is available
[here](https://www.nopredict.com/packages/tsgam.html).
