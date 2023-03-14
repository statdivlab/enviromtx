#' Fit model for fixed gene k, taxon j, taxon j*
#'
#' We would like to estimate beta1 in the model
#' $Y_{ijk} \sim Poisson(X_{ij} \times c \times (\frac{X_{ij*}}{X_{ij}})^{\beta_1})$
#' where
#' $Y_{ijk}$ refers to data about the abundance of gene k expressed by taxon j in sample i
#' $X_{ij}$ refers to data about the abundance of taxon j in sample i
#' $X_{ij*}$ refers to data about the abundance of taxon j* in sample i
#'
#' We can then interpret, e.g., $1.01^{\beta_1}$ as the multiplicative change in the expression-per-unit-coverage of gene k in species j for 1\% increase in the coverage of species j compared to species j*. With some slightly stronger assumptions about the sampling mechanism, we can also interpret this on the abundance (rather than just coverage) scale
#'
#' @param yy vector of abundances of gene k in taxon j
#' @param xstar vector of abundances of taxon j*
#' @param xx vector of abundances of taxon j
#' @param wts vector of nonnegative weights. Could be sequencing depth to put more emphasis on deeply sequenced samples
#' @param replace_zeros what to do with zeros in the denominator X_{ij}. Options include "minimum" or pseudocount numeric value (eg. 1)
#'
#' @importFrom tibble tibble
#' @importFrom rigr regress
#'
#'
#' @export
fit_mgx_model <- function(yy, xstar, xx, wts = NULL, replace_zeros = "minimum") {

  if (replace_zeros == "minimum") {
    xx[xx == 0] <- min(xx[xx > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    xx[xx == 0] <- replace_zeros # cheap pseudocount - take minimum
  }

  if (any(xx == 0)) stop("There are zeros in xx and you haven't told me what to do with them!")


  # response <- yy / xx
  predictor <- log(xstar / xx)

  # df <- tibble(yy, xstar, xx, response, predictor)
  df <- tibble(yy, xstar, xx, predictor)

  if(any(apply(df, 2, is.infinite))) {
    print(df)
    print(df, n = nrow(df))
    stop("Infinities in df?")
  }

  #
  if (is.null(wts)) {
    wts <- rep(1, nrow(data.frame(df)))
  }

  # Poisson regression
  rigr_out <- suppressWarnings(rigr::regress("rate",
                                             formula = yy ~ predictor,
                                             offset=log(xx),
                                             data=df,
                                             weights=wts,
                                             intercept=TRUE,
                                             exponentiate=FALSE,
                                             robustSE=TRUE))
  ## this gives raw model, untransformed coefficients

  # Grab rob
  simplified_output <- rigr_out$model[ , c("Estimate", "Robust SE", "F stat", "Pr(>F)")]
  colnames(simplified_output) <- c("Estimate", "Robust SE", "Test Statistic", "p-value")

  simplified_output[2, ]

}


