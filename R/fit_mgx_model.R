#' Estimate the difference in expression-per-unit-abundance of a gene expressed by a responder taxon as a function of the abundance of a companion taxon
#'
#' We would like to estimate the \eqn{\beta} parameters in the model
#' \eqn{\text{average } Y_{ijk} = X_{ij} \times c \times (X_{ij^*} / X_{ij})^{\beta_1} \times e^{\beta_2 W_{i1} + \beta_3 W_{i2} + \ldots + \beta_p W_{i,p-1}}},
#' where
#' \eqn{Y_{ijk}} refers to the observed abundance of gene \eqn{k} expressed by taxon \eqn{j} in sample \eqn{i},
#' \eqn{X_{ij}} refers to the observed abundance of taxon \eqn{j} in sample \eqn{i},
#' \eqn{X_{ij^*}} refers to the observed abundance of taxon \eqn{j^*} in sample \eqn{i}, and
#' \eqn{W_{im}} is the value of abiotic covariate \eqn{m} observed in sample \eqn{i}.
#' We can interpret \eqn{1.01^{\beta_1}} as the multiplicative change in the expression-per-unit-coverage of gene \eqn{k} in species \eqn{j} for 1\\% increase in the coverage of species \eqn{j} compared to species \eqn{j^*}. With some slightly stronger assumptions about the sampling mechanism, we can also interpret this on the true cell abundance scale (rather than just the coverage scale).
#'
#' @param enviro_df A data frame or tibble with columns containing the relevant variables, and rows denoting the observations. Columns must contain the following data: observed gene expression data, taxon abundance data, environmental covariates, and (optionally) technical replicate information.
#' @param yy The name of the column in `enviro_df` containing the expression data for a gene expressed by the responder taxon. Default is "yy". Expression data could be in the form of coverage, counts, etc. See vignettes for details on acceptable and suggested datatypes and preprocessing ("normalizations"). In the above notation, this corresponds to the expression of gene \eqn{k} expressed by taxon \eqn{j}.
#' @param xx The name of the column in `enviro_df` containing abundances of the responder taxon. Default is "xx". In the above notation, this corresponds to the abundance of taxon \eqn{j}.
#' @param xstar The name of the column in `enviro_df` containing abundances of the companion taxon. Default is "xstar". In the above notation, this corresponds to the abundance of  taxon \eqn{j^*}.
#' @param replicates The name of the column in `enviro_df` containing technical replicate information. Entries in this column should be the same within observations that are technical replicates. Will be coerced to a factor internally. Default is `NULL`.
#' @param formula A formula describing the environmental covariates that should be included in the model. The variable names should be columns in `enviro_df`
#' @param wts The name of the column in `enviro_df` containing nonnegative weights for the observation. Could be sequencing depth to put more emphasis on deeply sequenced samples. Default is `NULL`.
#' @param replace_zeros What to replace zeros with in the denominator \eqn{X_{ij}}. Options include "minimum" (replace with the smallest value) or pseudocount numeric value (eg. 1). Default is "minimum".
#' @param use_jack_se Boolean indicating whether to use jackknife standard errors (TRUE) rather than sandwich standard errors (FALSE). Only relevant if technical replicates are included. Jackknife standard errors are recommended when
#' there are a small number of observations. Defaults to FALSE.
#' @param cluster_corr_coef When technical replicates are given, the estimated value of the within-cluster correlation coefficient. This will only be used when GEE estimation in `raoBust::gee_test` fails, and
#' estimation is performed with a glm. Defaults to NULL by default (no robust score test is returned).
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

  if (any(enviro_df[[xstar]] < 0)) stop("You gave negative abundance data for `xstar`. This doesn't make sense.")
  if (any(enviro_df[[xx]] < 0)) stop("You gave negative abundance data for `xx`. This doesn't make sense.")
  if (any(enviro_df[[yy]] < 0)) stop("You gave negative expression data `yy`. This doesn't make sense.")


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


