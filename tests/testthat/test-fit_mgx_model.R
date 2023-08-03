n <- 10
test_that("default model runs", {

  set.seed(1)
  output1 <- fit_mgx_model(yy = rpois(n, lambda=100),
                           xstar = rpois(n, lambda=100),
                           xx = rpois(n, lambda=100))
  output1
  expect_type(output1, "double")

})

test_that("weights work", {

  set.seed(2)
  output2 <- fit_mgx_model(yy = rpois(n, lambda=100),
                           xstar = rpois(n, lambda=100),
                           xx = rpois(n, lambda=100),
                           wts = rpois(n, lambda=1000))
  output2
  expect_type(output2, "double")

})



test_that("psuedocounts work", {

  set.seed(2)
  xx1 <- rpois(n, lambda=10)
  xx1[2] <- 0

  output2 <- fit_mgx_model(yy = rpois(n, lambda=100),
                           xstar = rpois(n, lambda=100),
                           xx = xx1,
                           replace_zeros=1,
                           wts = rpois(n, lambda=100))

  expect_type(output2, "double")

  output3 <- fit_mgx_model(yy = rpois(n, lambda=100),
                           xstar = rpois(n, lambda=100),
                           xx = xx1,
                           wts = rpois(n, lambda=100))

  expect_type(output2, "double")

  yy1 <- rpois(n, lambda=100)
  xstar1 <- rpois(n, lambda=100)
  output4 <- fit_mgx_model(yy1,
                           xstar1,
                           xx = xx1,
                           replace_zeros="minimum")

  output5 <- fit_mgx_model(yy = yy1,
                           xstar = xstar1,
                           xx = xx1,
                           replace_zeros=min(xx1[xx1>0]))

  expect_equal(output4, output5)

})


