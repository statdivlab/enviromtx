my_df <- data.frame(xx = rpois(20, lambda = 400),
                    xstar = rpois(20, lambda = 400),
                    treatment = rep(0:1, each = 10),
                    id = rep(1:4, 5),
                    wt = 1:20)
my_df$yy <- rpois(20, my_df$xx * 10 * (my_df$xstar / my_df$xx)^2)

test_that("when a value of wts or replicates are missing, error occurs", {

  my_df$wt[1:2] <- NA
  expect_error(fit_mgx_model(enviro_df = my_df, replicates = "id",
                             formula = ~ treatment, wts = "wt"),
               "At least one of your weights provided is NA. Please provide a value for this weight, do not use weights, or remove observations with missing weights.")

  my_df$wt[1:2] <- 1:2
  my_df$id[5] <- NA
  expect_error(fit_mgx_model(enviro_df = my_df, replicates = "id",
                             formula = ~ treatment),
               "At least one of your replicates provided is NA. Please provide a value for this replicate, do not use replicates, or remove observations with missing replicates.")

})

test_that("when a value of yy, xx, xstar, or a relevant covariate is missing, that row is dropped", {

  # missing yy value
  my_df_yy <- my_df
  my_df_yy$yy[2] <- NA

  glm_res <- fit_mgx_model(my_df_yy)
  glm_res_no2 <- fit_mgx_model(my_df_yy[-2, ])
  expect_true(all.equal(glm_res, glm_res_no2))

  gee_res <- fit_mgx_model(my_df_yy, replicates = "id")
  gee_res_no2 <- fit_mgx_model(my_df_yy[-2, ], replicates = "id")
  expect_true(all.equal(gee_res, gee_res_no2))

  # missing xx value
  my_df_xx <- my_df
  my_df_xx$xx[2] <- NA

  glm_res <- fit_mgx_model(my_df_xx)
  glm_res_no2 <- fit_mgx_model(my_df_xx[-2, ])
  expect_true(all.equal(glm_res, glm_res_no2))

  gee_res <- fit_mgx_model(my_df_xx, replicates = "id")
  gee_res_no2 <- fit_mgx_model(my_df_xx[-2, ], replicates = "id")
  expect_true(all.equal(gee_res, gee_res_no2))

  # missing xstar value
  my_df_xstar <- my_df
  my_df_xstar$xstar[2] <- NA

  glm_res <- fit_mgx_model(my_df_xstar)
  glm_res_no2 <- fit_mgx_model(my_df_xstar[-2, ])
  expect_true(all.equal(glm_res, glm_res_no2))

  gee_res <- fit_mgx_model(my_df_xstar, replicates = "id")
  gee_res_no2 <- fit_mgx_model(my_df_xstar[-2, ], replicates = "id")
  expect_true(all.equal(gee_res, gee_res_no2))

  # missing covariate value
  my_df_cov <- my_df
  my_df_cov$cov <- rnorm(nrow(my_df))
  my_df_cov$cov[2] <- NA

  glm_res <- fit_mgx_model(my_df_cov, formula = ~ cov)
  glm_res_no2 <- fit_mgx_model(my_df_cov[-2, ], formula = ~ cov)
  expect_true(all.equal(glm_res, glm_res_no2))

  gee_res <- fit_mgx_model(my_df_cov, replicates = "id", formula = ~ cov)
  gee_res_no2 <- fit_mgx_model(my_df_cov[-2, ], replicates = "id", formula = ~ cov)
  expect_true(all.equal(gee_res, gee_res_no2))

  # missing irrelevant covariate value
  my_df_irr_cov <- my_df
  my_df_irr_cov$cov <- rnorm(nrow(my_df))
  my_df_irr_cov$irr_cov <- rnorm(nrow(my_df))
  my_df_irr_cov$irr_cov[2] <- NA

  glm_res <- fit_mgx_model(my_df_irr_cov, formula = ~ cov)
  glm_res_oth <- fit_mgx_model(subset(my_df_irr_cov, select = -irr_cov),
                               formula = ~ cov)
  expect_true(all.equal(glm_res, glm_res_oth))

  gee_res <- fit_mgx_model(my_df_irr_cov, formula = ~ cov,
                           replicates = "id")
  gee_res_oth <- fit_mgx_model(subset(my_df_irr_cov, select = -irr_cov),
                               formula = ~ cov, replicates = "id")
  expect_true(all.equal(gee_res, gee_res_oth))
})

