test_that("reasonably accurate estimates with covariates and correlation", {

  set.seed(7)
  m <- 5
  n_clust <- 200
  sigma_b <- 0.2
  n <- n_clust * m
  #### generate observations to be cluster correlated via random effect
  id  <- rep(1:n_clust, each = m)
  b <- rnorm(n_clust, mean = 0, sd = sigma_b)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 1
  beta1 <- 1
  beta2 <- 0.1
  beta3 <- -0.1
  xx_covariates1 <- rnorm(n)
  xx_covariates2 <- rnorm(n)
  eta <- xx1 * beta0 * (xstar1/xx1)^beta1 * exp(beta2 * xx_covariates1 + beta3 * xx_covariates2 + b[id])
  # yy1 <- rpois(n, lambda=eta)
  yy1 <- rnbinom(n, mu=eta, size=eta^2)

  my_df <- tibble(yy1, xx1, predictor = log((xstar1)/(xx1)), xx_covariates1, xx_covariates2, id)

  output10 <- fit_mgx_model(yy = yy1,
                            xstar = xstar1,
                            xx = xx1,
                            formula= ~ xx_covariates1 + xx_covariates2,
                            enviro_df=cbind(xx_covariates1, xx_covariates2),
                            replicates=id)

  expect_equal(output10["predictor", "Estimate"], expected=beta1, tolerance=0.05)
  expect_equal(output10["xx_covariates1", "Estimate"], expected=beta2, tolerance=0.05)


})


