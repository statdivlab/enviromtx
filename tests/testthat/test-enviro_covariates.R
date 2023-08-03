test_that("environmental covariates work", {
  yy_i_k <- c(6, 6, 0, 15.41786, 9.36569, 9.65386, 1, 26.68171, 8.46619,
              0, 3.26028, 4.27604, 21.7861, 18.81282, 33.24853886644, 0, 0,
              0, 0, 0, 2.6214, 1, 2, 2, 0, 0, 0, 1, 6, 4)
  nn <- length(yy_i_k)

  xstar_i <- c(201.68944, 127.3064219, 92.152581, 237.559910010157, 158.1550348169, 312.489719,
               880.352877, 899.822609, 1593.640589, 1503.400459, 1583.012547, 2105.019572,
               600.5897, 750.13492, 674.74847586644, 234.585574,
               312.5290412931, 226.34653, 195.221137077546, 278.31226,
               233.5987703, 70.57603, 111.66915903, 123.100246,
               205.701248, 243.516225, 217.50233, 206.28941, 320.42961, 307.036039)

  xx_i <- xstar_i + rexp(nn, 1)

  temp <- rnorm(nn, mean=23, sd = 2)
  salinity <- rnorm(nn, mean=35, sd = 5)

  out_i_k <- fit_mgx_model(yy_i_k, xstar_i, xx_i,
                           formula = ~ temp + salinity,
                           enviro_df = cbind(temp, salinity))

  expect_type(out_i_k, "double")

  predictor <- log(xstar_i / xx_i)
  glm_coefs <- suppressWarnings(coef(glm(yy_i_k ~ predictor + temp + salinity,
                                         family=poisson(link="log"), offset=log(xx_i)))[-1])
  glm_coefs
  expect_true(all(out_i_k[, 1] == glm_coefs))

})
