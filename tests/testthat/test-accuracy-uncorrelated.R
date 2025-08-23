n <- 20
test_that("reasonably accurate estimates", {

  # accuracy test 1
  set.seed(3)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 100
  beta1 <- 1
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)

  output6 <- fit_mgx_model(yy = yy1,
                           xstar = xstar1,
                           xx = xx1,
                           replace_zeros=1)

  expect_equal(unname(output6["Estimate"]), expected=beta1, tolerance=abs(0.05*beta1))

  # accuracy test 2
  xx1 <- rpois(n, lambda=100)
  xstar1 <- rpois(n, lambda=20)
  beta0 <- 1
  beta1 <- -1
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)
  yy1

  output7 <- fit_mgx_model(yy = yy1,
                           xstar = xstar1,
                           xx = xx1,
                           replace_zeros=1)

  expect_equal(unname(output7["Estimate"]), expected=beta1, tolerance=abs(0.05*beta1))

})


test_that("reasonably accurate estimates with covariates", {

  set.seed(4)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 10
  beta1 <- 2
  beta2 <- 1
  beta3 <- -1
  xx_covariates1 <- rnorm(n)
  xx_covariates2 <- rnorm(n)
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1 * exp(beta2 * xx_covariates1 + beta3 * xx_covariates2))

  yy1
  output8 <- fit_mgx_model(yy = yy1,
                           xstar = xstar1,
                           xx = xx1,
                           formula= ~ xx_covariates1 + xx_covariates2,
                           enviro_df=cbind(xx_covariates1, xx_covariates2),
                           replace_zeros=1)

  expect_equal(output8["predictor", "Estimate"], expected=beta1, tolerance=abs(0.05*beta1))
  expect_equal(output8["xx_covariates1", "Estimate"], expected=beta2, tolerance=abs(0.05*beta2))
  expect_equal(output8["xx_covariates2", "Estimate"], expected=beta3, tolerance=abs(0.05*beta3))

})

test_that("reasonably accurate estimates with covariates and correlation", {

  set.seed(5)
  m <- 3
  n_clust <- 50
  sigma_b <- 0
  n <- n_clust * m
  #### generate observations to be cluster correlated via random effect
  id  <- rep(1:n_clust, each = m)
  b <- rnorm(n_clust, mean = 0, sd = sigma_b)

  set.seed(3)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 1
  beta1 <- 2
  beta2 <- 1
  beta3 <- -1
  xx_covariates1 <- rnorm(n)
  xx_covariates2 <- rnorm(n)
  eta <- xx1 * beta0 * (xstar1/xx1)^beta1 * exp(beta2 * xx_covariates1 + beta3 * xx_covariates2 + b[id])
  # yy1 <- rpois(n, lambda=eta)
  yy1 <- rnbinom(n, mu=eta, size=eta^2)

  my_df <- tibble(yy1, xx1, predictor = log((xstar1 + 1)/(xx1 + 1)), xx_covariates1, xx_covariates2, id)

  output9 <- fit_mgx_model(yy = yy1,
                xstar = xstar1,
                xx = xx1,
                formula= ~ xx_covariates1 + xx_covariates2,
                enviro_df=cbind(xx_covariates1, xx_covariates2),
                replicates=id,
                replace_zeros=1)

  expect_equal(output9["predictor", "Estimate"], expected=beta1, tolerance=abs(0.05*beta1))
  expect_equal(output9["xx_covariates1", "Estimate"], expected=beta2, tolerance=abs(0.05*beta2))
  expect_equal(output9["xx_covariates2", "Estimate"], expected=beta3, tolerance=abs(0.05*beta3))

})


