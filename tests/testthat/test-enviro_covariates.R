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

  df <- data.frame(temp, salinity,
              yy = yy_i_k,
              xstar = xstar_i,
              xx = xx_i)

  expect_silent(out_i_k <- fit_mgx_model(df,
                                         formula = ~ temp + salinity))

  expect_type(out_i_k[1,1], "double")

  predictor <- log(xstar_i / xx_i)
  glm_coefs <- suppressWarnings(coef(glm(yy_i_k ~ predictor + temp + salinity,
                                         family=poisson(link="log"), offset=log(xx_i)))[-1])
  glm_coefs
  expect_true(all(out_i_k[1:3, 1] == glm_coefs))


  expect_warning(fit_mgx_model(df))


  ### test with replicates
  df$replicates <- rep(LETTERS[1:10], each = 3)
  expect_silent(out_i_k_rep <- fit_mgx_model(df,
                                             formula = ~ temp + salinity,
                                             replicates="replicates"))
  expect_type(out_i_k_rep[1,1], "double")


  # test with infinite covariate value
  new_df <- df
  new_df$temp[3] <- Inf
  expect_error(fit_mgx_model(new_df, formula = ~ temp + salinity),
               "At least one value of covariate temp is infinite. Please fix this and rerun.")

  # test with infinite replicate value
  new_df$id <- rep(1:5, 6)
  new_df$id[5] <- -Inf
  expect_error(fit_mgx_model(new_df, formula = ~ temp + salinity, replicates = "id"),
               "At least one replicate provided is infinite. Please fix this and then rerun.")

  new_df <- df
  new_df$id <- rnorm(30)
  expect_error(fit_mgx_model(new_df, formula = ~ temp + salinity + id),
               "You have a covariate called 'id' in your formula. This is a protected term in `fit_mgx_model`, please use a different name for this covariate")

})
