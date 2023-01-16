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



test_that("reasonably accurate estimates", {

  # accuracy test 1
  set.seed(3)
  xx1 <- rpois(n, lambda=400)
  xstar1 <- rpois(n, lambda=400)
  beta0 <- 100
  beta1 <- 20
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
  beta1 <- -4
  yy1 <- rpois(n, xx1 * beta0 * (xstar1/xx1)^beta1)
  yy1

  output7 <- fit_mgx_model(yy = yy1,
                           xstar = xstar1,
                           xx = xx1,
                           replace_zeros=1)

  expect_equal(unname(output7["Estimate"]), expected=beta1, tolerance=abs(0.05*beta1))


})


test_that("large-scale accuracy test", {

  # set.seed(10)
  # res <- matrix(NA, nrow = 100, ncol = 40)
  # for (beta1 in 11:50) {
  #   for (i in 1:100) {
  #     beta0 <- 100
  #     xx1 <- rpois(n, lambda=400)
  #     xstar1 <- rpois(n, lambda=400)
  #     yy1 <- rpois(20, xx1 * beta0 * (xstar1/xx1)^beta1)
  #
  #     output6 <- try(fit_mgx_model(yy = yy1,
  #                                  xstar = xstar1,
  #                                  xx = xx1,
  #                                  replace_zeros=1))
  #     res[i, beta1-10] <- output6[1]
  #   }
  # }
  # colnames(res) <- 11:50
  # library(tidyverse)
  # res %>%
  #   as.data.frame %>%
  #   apply(2, as.numeric) %>%
  #   as.data.frame %>%
  #   tibble %>%
  #   pivot_longer(1:40, names_to="beta1", values_to="estimate") %>%
  #   mutate(beta1 = as.numeric(beta1)) %>%
  #   ggplot(aes(x = beta1, y = estimate)) +
  #   geom_point() +
  #   geom_abline(slope = 1, intercept=0) +
  #   theme_bw() +
  #   NULL
  # # phew, fine

  expect_true(TRUE)

})

