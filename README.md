
<!-- README.md is generated from README.Rmd. Please edit that file -->

# enviromtx

<!-- badges: start -->
<!-- badges: end -->

The goal of `enviromtx` is to answer the question: Does the abundance of
a given species impact another species’ expression of a gene?

## Installation

You can install the development version of enviromtx like so:

``` r
remotes::install_github("statdivlab/enviromtx")
```

## Example

This is a basic example which shows you how to fit the model. More
documentation is coming soon!

``` r
library(enviromtx)
n <- 10
set.seed(3)
xx1 <- rpois(n, lambda=400)
xstar1 <- rpois(n, lambda=400)
beta0 <- 100
beta1 <- 20
yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)

fit_mgx_model(yy = yy1,
              xstar = xstar1,
              xx = xx1,
              replace_zeros=1)
#>             Estimate Non-robust Std Error     Robust Std Error 
#>          19.97137412           0.02456755           0.01916404 
#>    Non-robust Wald p        Robust Wald p       Robust Score p 
#>           0.00000000           0.00000000           0.01413211
```
