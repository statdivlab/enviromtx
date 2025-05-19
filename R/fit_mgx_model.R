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
#' @param replicates vector that contains the same entries for all observations that are technical replicates. Will be coerced to a factor down the line.
#' @param formula a formula describing the environmental covariates that should be included in the model. The variable names should be columns in enviro_df
#' @param enviro_df a data frame or tibble with columns containing environmental covariates that should be included in the model,
#' @param wts vector of nonnegative weights. Could be sequencing depth to put more emphasis on deeply sequenced samples
#' @param replace_zeros what to do with zeros in the denominator X_{ij}. Options include "minimum" or pseudocount numeric value (eg. 1)
#' @param use_jack_se when replicates are given, if TRUE will use jackknife standard errors instead of sandwich standard errors. This is recommended when
#' there is a small number of clusters.
#' @param cluster_corr_coef when replicates are given, estimated value of the within-cluster correlation coefficient. This will only be used when gee estimation in `raoBust::gee_test` fails, and instead
#' estimation is performed with a glm. This is set to NULL by default.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @importFrom rigr regress
#' @importFrom lazyeval f_new as_call
#' @importFrom raoBust glm_test
#' @import stats
#' @import geepack
#' @import geeasy
#'
#'
#' @export
fit_mgx_model <- function(
    yy,
    xstar,
    xx,
    replicates = NULL,
    formula = NULL,
    enviro_df = NULL,
    wts = NULL,
    replace_zeros = "minimum",
    use_jack_se = FALSE,
    cluster_corr_coef = NULL
) {

  # Setting pseudocounts for xx and xstar
  if (replace_zeros == "minimum") {
    xx[xx == 0] <- min(xx[xx > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    xx[xx == 0] <- replace_zeros
  }
  if (any(xx == 0)) stop("There are zeros in xx and you haven't told me what to do with them!")

  if (replace_zeros == "minimum") {
    xstar[xstar == 0] <- min(xstar[xstar > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    xstar[xstar == 0] <- replace_zeros
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
    my_df <- dplyr::bind_cols(tibble::tibble(yy, xstar, xx, predictor), enviro_df)

  } else if (xor(is.null(formula), is.null(enviro_df))) {

    stop("One of `formula` or `enviro_df` is missing. Please correct this and try again.")

  } else if (is.null(formula) & is.null(enviro_df)) {

    my_df <- tibble::tibble(yy, xstar, xx, predictor)
    my_formula <- yy ~ predictor


  } else {
    stop("Ergh, there is an error and it is Amy's fault.")
  }

  if(any(apply(my_df, 2, is.infinite))) {
    print(my_df, n = nrow(my_df))
    stop("Infinities in my_df? Amy's fault.")
  }

  if (is.null(wts)) {
    wts <- rep(1L, nrow(data.frame(my_df)))
  }
  my_df$wts <- wts

  ################################################
  ## fit the model with Poisson regression
  ################################################

  ## model is
  ## log(E[Y]) = log(X) + beta0 + beta1 * log(X^*/X) + beta2 * salinity + beta3 * iron
  ## E[Y]/X = gamma0 * (X^*/X)^beta1 * e^(beta2*salinity + beta3*iron)

  # raoBust_out <- eval(rlang::expr(raoBust::glm_test( TODO )),
  #                     envir = .env)

  if (is.null(replicates)) {
    raoBust_out <- raoBust::glm_test(formula = my_formula,
                                     offset=log(xx),
                                     family=stats::poisson(link="log"),
                                     data=my_df,
                                     weights=wts)$coef_tab
  } else {
    if (!all(wts == 1L)) {
      warning("Run this by Amy; not sure what this is doing off-the-cuff")
    }
    my_df$id <- replicates

    raoBust_out <- raoBust::gee_test(formula = my_formula,
                                     offset=log(xx),
                                     family=stats::poisson(link="log"),
                                     id=id,
                                     # weights=wts,
                                     data=my_df,
                                     use_jack_se = use_jack_se,
                                     cluster_corr_coef = cluster_corr_coef)$coef_tab
  }


  ## this gives raw model, untransformed coefficients

  raoBust_out[-1, ]

}


