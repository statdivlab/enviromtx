#' Fit model for fixed gene k, taxon j, taxon j*
#'
#' We would like to estimate beta1 in the model
#' Y_ijk ~ Poisson(X_ij * c * (X_{ij*} / X_{ij})^{beta_1})
#' where
#' $Y_{ijk}$ refers to data about the abundance of gene k expressed by taxon j in sample i
#' $X_{ij}$ refers to data about the abundance of taxon j in sample i
#' $X_{ij*}$ refers to data about the abundance of taxon j* in sample i
#'
#' We can then interpret, e.g., $1.01^{beta_1}$ as the multiplicative change in the expression-per-unit-coverage of gene k in species j for 1\% increase in the coverage of species j compared to species j*. With some slightly stronger assumptions about the sampling mechanism, we can also interpret this on the abundance (rather than just coverage) scale
#'
#' @param yy vector of abundances of gene k in taxon j
#' @param xstar vector of abundances of taxon j*
#' @param xx vector of abundances of taxon j
#' @param formula a formula describing the environmental covariates that should be included in the model. The variable names should be columns in enviro_df
#' @param enviro_df a data frame or tibble with columns containing environmental covariates that should be included in the model,
#' @param wts vector of nonnegative weights. Could be sequencing depth to put more emphasis on deeply sequenced samples
#' @param replace_zeros what to do with zeros in the denominator X_{ij}. Options include "minimum" or pseudocount numeric value (eg. 1)
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @importFrom rigr regress
#' @importFrom lazyeval f_new as_call
#' @importFrom raoBust glm_test
#'
#'
#' @export
fit_mgx_model <- function(yy, xstar, xx,
                          formula = NULL,
                          enviro_df = NULL,
                          wts = NULL, replace_zeros = "minimum") {

  # Setting pseudocounts for xx and xstar
  if (replace_zeros == "minimum") {
    xx[xx == 0] <- min(xx[xx > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    xx[xx == 0] <- replace_zeros # cheap pseudocount - take minimum
  }
  if (any(xx == 0)) stop("There are zeros in xx and you haven't told me what to do with them!")

  if (replace_zeros == "minimum") {
    xstar[xstar == 0] <- min(xstar[xstar > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    xstar[xstar == 0] <- replace_zeros # cheap pseudocount - take minimum
  }

  if (any(xstar == 0)) stop("There are zeros in xstar and you haven't told me what to do with them!")

  ################################################
  ## set up the data for fitting
  ################################################

  predictor <- log(xstar / xx)

  if (!is.null(formula) & !is.null(enviro_df)) {

    ## formula is yy ~ log(x*/x) + salinity + iron
    my_formula <- lazyeval::f_new(lhs=quote(yy),
                                  rhs = lazyeval::as_call(paste("predictor +", as.character(formula)[2])))
    df <- dplyr::bind_cols(tibble::tibble(yy, xstar, xx, predictor), enviro_df)

  } else if (xor(is.null(formula), is.null(enviro_df))) {

    stop("One of `formula` or `enviro_df` is missing. Please correct this and try again.")

  } else if (is.null(formula) & is.null(enviro_df)) {

    df <- tibble::tibble(yy, xstar, xx, predictor)
    my_formula <- yy ~ predictor

  } else {
    stop("Ergh, there is an error and it is Amy's fault.")
  }

  if(any(apply(df, 2, is.infinite))) {
    print(df, n = nrow(df))
    stop("Infinities in df? Amy's fault.")
  }

  if (is.null(wts)) {
    wts <- rep(1, nrow(data.frame(df)))
  }
  df$wts <- wts

  ################################################
  ## fit the model with Poisson regression
  ################################################

  ## model is
  ## log(E[Y]) = log(X) + beta0 + beta1 * log(X^*/X) + beta2 * salinity + beta3 * iron
  ## E[Y]/X = gamma0 * (X^*/X)^beta1 * e^(beta2*salinity + beta3*iron)

  # rigr_out <- suppressWarnings(rigr::regress("rate",
  #                                            formula = my_formula,
  #                                            offset=log(xx),
  #                                            data=df,
  #                                            weights=wts,
  #                                            intercept=TRUE,
  #                                            exponentiate=FALSE,
  #                                            robustSE=TRUE))
  # stop("no")



  # raoBust_out <- coef(summary(glm(formula = my_formula,
  #                      offset=log(xx),
  #                      family=poisson(link="log"),
  #                      data=df,
  #                      weights=wts)))
  raoBust_out <- raoBust::glm_test(formula = my_formula,
                                   offset=log(xx),
                                   family=poisson(link="log"),
                                   data=df,
                                   weights=wts)

  # print(raoBust_out)

  ## this gives raw model, untransformed coefficients

  # Grab rob
  # simplified_output <- rigr_out$model[ , c("Estimate", "Robust SE", "F stat", "Pr(>F)")]
  # colnames(simplified_output) <- c("Estimate", "Robust SE", "Test Statistic", "p-value")



  raoBust_out[-1, ]

}


