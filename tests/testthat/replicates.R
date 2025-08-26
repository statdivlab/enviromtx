# install.packages("/Users/adwillis/software/raoBust_0.0.2.1.tar.gz", repos = NULL, type = "source")

test_that("runs with replicates", {

  set.seed(1)
  n <- 20
  xx1 <- rpois(n, lambda=100)
  xstar1 <- rpois(n, lambda=20)
  beta0 <- 1
  beta1 <- -4
  reps <- rep(1:4, each = 5)

  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)

  ## works fine
  expect_type(fit_mgx_model(data.frame(yy = yy1,
                xstar = xstar1,
                xx = xx1),
                replace_zeros=1),
              "double")

  ### does not work
  # Error in eval(cl$data) : object 'my_df' not found
  expect_type(fit_mgx_model(data.frame(yy = yy1,
                xstar = xstar1,
                xx = xx1, reps),
                replicates = "reps",
                replace_zeros=1),
              "list")

})

