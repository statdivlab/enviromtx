test_that("with small parameter values, centering does not affect results", {

  n <- 20

  # no replicates
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
  output1 <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1,
                                      xx_covariates1,
                                      xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1)
  output2 <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1,
                                      xx_covariates1,
                                      xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1,
                           control = list(center = FALSE))

  expect_true(all.equal(output1, output2))

  # with replicates
  set.seed(5)
  m <- 3
  n_clust <- 50
  sigma_b <- 0
  n <- n_clust * m
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
  yy1 <- rnbinom(n, mu=eta, size=eta^2)

  my_df <- tibble(yy1, xx1, predictor = log((xstar1 + 1)/(xx1 + 1)), xx_covariates1, xx_covariates2, id)

  output3 <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1, xx_covariates1, xx_covariates2, id),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replicates="id",
                           replace_zeros=1,
                           control = list(center = TRUE))
  output4 <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1, xx_covariates1, xx_covariates2, id),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replicates="id",
                           replace_zeros=1)
  expect_true(all.equal(output3, output4, tolerance = 1e-4))

})

test_that("with large parameter values, centering improves estimation", {

  # turns out this does not improve estimation, estimation is the same but less stable in some of these
  # cases after centering

  skip("Skip longer sim study in automatic testing")

  n <- 20

  # large B values
  set.seed(4)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 10
  beta1 <- 20
  beta2 <- 30
  beta3 <- -1
  xx_covariates1 <- rnorm(n)
  xx_covariates2 <- rnorm(n)
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1 * exp(beta2 * xx_covariates1 + beta3 * xx_covariates2))
  output1 <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1,
                                      xx_covariates1,
                                      xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1,
                           control = list(center = TRUE))
  output2 <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1,
                                      xx_covariates1,
                                      xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1)

  expect_true(all.equal(output1, output2))


  # large covariate values
  n <- 20
  set.seed(3)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 1.5
  beta1 <- -1
  beta2 <- 1
  beta3 <- 0
  xx_covariates1 <- rnorm(n) + 20
  xx_covariates2 <- rnorm(n)
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1 * exp(beta2 * xx_covariates1 + beta3 * xx_covariates2))
  output_center <- fit_mgx_model(data.frame(yy = yy1,
                                            xstar = xstar1,
                                            xx = xx1,
                                            xx_covariates1,
                                            xx_covariates2),
                                 formula= ~ xx_covariates1 + xx_covariates2,
                                 replace_zeros=1,
                                 control = list(center = TRUE))
  output_nocenter <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1,
                                      xx_covariates1,
                                      xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1)

  expect_true(all.equal(output1, output2))

})

n <- 20
test_that("reasonably accurate estimates", {

  # accuracy test 1
  set.seed(3)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 100
  beta1 <- 20
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)

  output6 <- fit_mgx_model(data.frame(yy = yy1,
                           xstar = xstar1,
                           xx = xx1),
                           replace_zeros=1,
                           control = list(center = TRUE))
  output6_noncenter <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1),
                           replace_zeros=1,
                           control = list(center = FALSE))

  expect_equal(output6, output6_noncenter)

  # accuracy test 2
  xx1 <- rpois(n, lambda=100)
  xstar1 <- rpois(n, lambda=20)
  beta0 <- 1
  beta1 <- -4
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)
  yy1

  output7 <- fit_mgx_model(data.frame(yy = yy1,
                           xstar = xstar1,
                           xx = xx1),
                           replace_zeros=1,
                           control = list(center = TRUE))
  output7_nocenter <- fit_mgx_model(data.frame(yy = yy1,
                                      xstar = xstar1,
                                      xx = xx1),
                           replace_zeros=1,
                           control = list(center = FALSE))

  expect_true(all.equal(output7, output7_nocenter))

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
  output8 <- fit_mgx_model(cbind(data.frame(yy = yy1,
                           xstar = xstar1,
                           xx = xx1), xx_covariates1, xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1,
                           control = list(center = TRUE))
  output8_nocenter <- fit_mgx_model(cbind(data.frame(yy = yy1,
                                            xstar = xstar1,
                                            xx = xx1), xx_covariates1, xx_covariates2),
                           formula= ~ xx_covariates1 + xx_covariates2,
                           replace_zeros=1,
                           control = list(center = FALSE))

  expect_equal(output8, output8_nocenter)

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
  beta2 <- 5
  beta3 <- -5
  xx_covariates1 <- rnorm(n)
  xx_covariates2 <- rnorm(n)
  eta <- xx1 * beta0 * (xstar1/xx1)^beta1 * exp(beta2 * xx_covariates1 + beta3 * xx_covariates2 + b[id])
  # yy1 <- rpois(n, lambda=eta)
  yy1 <- rnbinom(n, mu=eta, size=eta^2)

  my_df <- tibble(yy1, xx1, xstar1, xx_covariates1, xx_covariates2, id)

  output9 <- fit_mgx_model(yy = "yy1",
                           xstar = "xstar1",
                           xx = "xx1",
                           formula= ~ xx_covariates1 + xx_covariates2,
                           enviro_df=my_df,
                           replicates="id",
                           replace_zeros=1,
                           control = list(center = TRUE))
  output9_nocenter <- fit_mgx_model(yy = "yy1",
                           xstar = "xstar1",
                           xx = "xx1",
                           formula= ~ xx_covariates1 + xx_covariates2,
                           enviro_df=my_df,
                           replicates="id",
                           replace_zeros=1,
                           control = list(center = FALSE))

  expect_equal(output9, output9_nocenter, tolerance = 1e-4)

})
