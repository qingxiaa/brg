
<!-- README.md is generated from README.Rmd. Please edit that file -->

# brg

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/qingxiaa/brg.svg?branch=master)](https://travis-ci.com/qingxiaa/brg)
[![R-CMD-check](https://github.com/qingxiaa/brg/workflows/R-CMD-check/badge.svg)](https://github.com/qingxiaa/brg/actions)
<!-- badges: end -->

The goal of brg is to circumvent the limitations of existing
batch-effect correction methods studies involving the longitudinal
collection of high-dimensional ‘omic’ data, we propose Batch effect
coRrectIon of microarray data with Dependent samples usinG an Empirical
Bayes approach (BRIDGE). BRIDGE accounts for within-subject dependency
expected in longitudinal studies with high-dimensional ‘omic’ data and
involves the estimation and correction of additive batch and
multiplicative batch effects when batch effects are confounded with
time. After correcting for batch effects with BRIDGE, adjusted data can
be used for downstream statistical analyses as if all samples are
measured in a single batch. BRIDGE is applicable to different ‘omics’
platforms, such as microarray gene expression, DNA methylation or
neuroimaging studies of neurodevelopment, as long as the transformed
data are high-dimension and approximately normal distribution. \#\#
Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("qingxiaa/brg")
```
