#' Fit model for fixed gene k, taxon j, taxon j*
#'
#' @param yy vector of abundances of gene k in taxon j
#' @param xstar vector of abundances of taxon j*
#' @param xx vector of abundances of taxon j
#' @param wts vector of nonnegative weights. Could be sequencing depth to put more emphasis on deeply sequenced samples
#'
#' @export
fit_mgx_model <- function(yy, xstar, xx, wts = NULL) {

  xx[xx == 0] <- min(xx) # cheap psuedocount - take minimum

  response <- yy / xx
  predictor <- xstar / xx

  df <- tibble(yy, xstar, xx, response, predictor)

  # Poisson regression
  rigr_out <- rigr::regress("rate", response ~ predictor, data=df, weights=wts)

  rigr_out

}
