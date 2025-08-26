#' Fit model for fixed gene k, taxon j, taxon j*
#'
#' We would like to estimate beta1 in the model
#' \eqn{Y_ijk ~ Poisson(X_ij * c * (X_{ij*} / X_{ij})^{beta_1})}
#' where
#' \eqn{Y_{ijk}} refers to data about the abundance of gene k expressed by taxon j in sample i
#' \eqn{X_{ij}} refers to data about the abundance of taxon j in sample i
#' \eqn{X_{ij*}} refers to data about the abundance of taxon j* in sample i
#'
#' We can then interpret, e.g., \eqn{1.01^{beta_1}} as the multiplicative change in the expression-per-unit-coverage of gene k in species j for 1\% increase in the coverage of species j compared to species j*. With some slightly stronger assumptions about the sampling mechanism, we can also interpret this on the abundance (rather than just coverage) scale
#'
#' @param enviro_df a data frame or tibble with columns containing relevant abundances, environmental covariates that should be included in the model, and optionally replicate information
#' @param yy column of `enviro_df` containing abundances of gene k in taxon j. Default is "yy".
#' @param xstar column of `enviro_df` containing abundances of taxon j*. Default is "xstar".
#' @param xx column of `enviro_df` containing abundances of taxon j. Default is "xx".
#' @param replicates column of `enviro_df` containing the same entries for all observations that are technical replicates. Will be coerced to a factor down the line. Default is `NULL`.
#' @param formula a formula describing the environmental covariates that should be included in the model. The variable names should be columns in `enviro_df`
#' @param wts column of `enviro_df` containing nonnegative weights. Could be sequencing depth to put more emphasis on deeply sequenced samples. Default is `NULL`.
#' @param replace_zeros what to do with zeros in the denominator \eqn{X_{ij}}. Options include "minimum" or pseudocount numeric value (eg. 1)
#' @param use_jack_se when replicates are given, if TRUE will use jackknife standard errors instead of sandwich standard errors. This is recommended when
#' there is a small number of clusters.
#' @param cluster_corr_coef when replicates are given, estimated value of the within-cluster correlation coefficient. This will only be used when gee estimation in `raoBust::gee_test` fails, and instead
#' estimation is performed with a glm. This is set to NULL by default.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @importFrom lazyeval f_new as_call
#' @importFrom raoBust glm_test
#' @import stats
#' @import geepack
#' @import geeasy
#'
#' @examples
#' my_df <- data.frame(xx = rpois(20, lambda = 400),
#'                     xstar = rpois(20, lambda = 400),
#'                     treatment = rep(0:1, each = 10),
#'                     id = rep(1:4, 5))
#' my_df$yy <- rpois(20, my_df$xx * 10 * (my_df$xstar / my_df$xx)^2)
#' fit_mgx_model(enviro_df = my_df, replicates = "id", formula = ~ treatment)
#'
#' @export
fit_mgx_model <- function(
    enviro_df,
    yy = "yy",
    xstar = "xstar",
    xx = "xx",
    replicates = NULL,
    formula = NULL,
    wts = NULL,
    replace_zeros = "minimum",
    use_jack_se = FALSE,
    cluster_corr_coef = NULL
) {

  # first, check that everything is in enviro_df
  if (!(is.character(yy) &
        is.character(xstar) &
        is.character(xx))) {
    stop("`yy`, `xstar`, and `xx` must be character objects giving the name of the column in `enviro_df` containing each relevant vector, not the vectors themselves.")
  }
  missing_args <- which(!(c(yy, xstar, xx) %in% names(enviro_df)))
  if (length(missing_args) > 0) {
    stop(paste0("The column names ", c(yy, xstar, xx)[missing_args], " do not appear in enviro_df. Please include them and then rerun this function. "))
  }
  if (!is.null(replicates)) {
    if (!is.character(replicates)) {
      stop("`replicates` argument must be a column name in `enviro_df` if it is included.")
    } else {
      if (!(replicates %in% names(enviro_df))) {
        stop("`replicates` argument must be a column name in `enviro_df` if it is included.")
      }
    }
  }
  if (!is.null(wts)) {
    if (!is.character(wts)) {
      stop("`wts` argument must be a column name in `enviro_df` if it is included.")
    } else {
      if (!(wts %in% names(enviro_df))) {
        stop("`wts` argument must be a column name in `enviro_df` if it is included.")
      }
    }
  }

  # Setting pseudocounts for xx and xstar
  if (replace_zeros == "minimum") {
    enviro_df[[xx]][enviro_df[[xx]] == 0] <- min(enviro_df[[xx]][enviro_df[[xx]] > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    enviro_df[[xx]][enviro_df[[xx]] == 0] <- replace_zeros
  }
  if (any(enviro_df[[xx]] == 0)) stop("There are zeros in enviro_df$xx and you haven't told me what to do with them!")

  if (replace_zeros == "minimum") {
    enviro_df[[xstar]][enviro_df[[xstar]] == 0] <- min(enviro_df[[xstar]][enviro_df[[xstar]] > 0]) # cheap pseudocount - take minimum
  } else if (is.numeric(replace_zeros)) {
    if (replace_zeros == 0) stop("You've told me to replace zeroes with zero!")
    enviro_df[[xstar]][enviro_df[[xstar]] == 0] <- replace_zeros
  }

  if (any(enviro_df[[xstar]] == 0)) stop("There are zeros in xstar and you haven't told me what to do with them!")

  ################################################
  ## set up the data for fitting
  ################################################

  predictor <- log(enviro_df[[xstar]] / enviro_df[[xx]])

  if (is.null(formula)) {

    ## formula is yy ~ log(x*/x) + salinity + iron
    formula <- lazyeval::f_new(lhs=lazyeval::as_call(yy),
                               rhs = lazyeval::as_call("predictor"))
    my_df <- dplyr::bind_cols(tibble::tibble(predictor), enviro_df)

  } else {

    ## formula is yy ~ log(x*/x) + salinity + iron
    formula <- lazyeval::f_new(lhs=lazyeval::as_call(yy),
                               rhs = lazyeval::as_call(paste("predictor +", as.character(formula)[2])))
    my_df <- dplyr::bind_cols(tibble::tibble(predictor), enviro_df)

  }

  if(any(apply(my_df, 2, is.infinite))) {
    print(my_df, n = nrow(my_df))
    stop("Infinities in my_df? Amy's fault.")
  }

  if (is.null(wts)) {
    my_df$wts <- rep(1L, nrow(my_df))
    wts <- "wts"
  }

  my_df$offset <- log(my_df[[xx]])
  offset <- "offset"

  ################################################
  ## fit the model with Poisson regression
  ################################################

  ## model is
  ## log(E[Y]) = log(X) + beta0 + beta1 * log(X^*/X) + beta2 * salinity + beta3 * iron
  ## E[Y]/X = gamma0 * (X^*/X)^beta1 * e^(beta2*salinity + beta3*iron)

  if (is.null(replicates)) {
    raoBust_out <- raoBust::glm_test(formula = formula,
                                     offset = offset,
                                     family = stats::poisson(link = "log"),
                                     data = my_df,
                                     weights = wts)$coef_tab
  } else {
    if (!all(my_df[[wts]] == 1L)) {
      warning("Run this by Amy; not sure what this is doing off-the-cuff")
    }
    if (replicates != "id") {
      my_df$id <- my_df[[replicates]]
      my_df[[replicates]] <- NULL
      id <- "id"
    } else {
      id <- "id"
    }

    raoBust_out <- raoBust::gee_test(formula = formula,
                                     offset = offset,
                                     family = stats::poisson(link="log"),
                                     id = id,
                                     data = my_df,
                                     use_jack_se = use_jack_se,
                                     cluster_corr_coef = cluster_corr_coef)$coef_tab
  }


  ## this gives raw model, untransformed coefficients

  raoBust_out[-1, ]

}


