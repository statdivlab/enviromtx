test_that("runs with replicates", {

  set.seed(1)
  n <- 20
  xx1 <- rpois(n, lambda=100)
  xstar1 <- rpois(n, lambda=20)
  beta0 <- 1
  beta1 <- -4
  reps <- rep(1:4, each = 5)

  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)
  df <- data.frame(yy = yy1, xstar = xstar1, xx = xx1, reps = reps)

  ## works fine
  expect_type(fit_mgx_model(df,
                            replace_zeros=1),
              "double")

  expect_type(fit_mgx_model(df,
                            replicates = "reps",
                            replace_zeros=1),
              "list")

  # test using jackknife standard errors instead of sandwich standard errors
  expect_type(fit_mgx_model(df,
                            replicates = "reps",
                            replace_zeros=1,
                            use_jack_se = TRUE),
              "list")

})

