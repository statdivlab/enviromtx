test_that("debug example", {
  yy_i_k <- c(6, 6, 0, 15.41786, 9.36569, 9.65386, 1, 26.68171, 8.46619,
    0, 3.26028, 4.27604, 21.7861, 18.81282, 33.24853886644, 0, 0,
    0, 0, 0, 2.6214, 1, 2, 2, 0, 0, 0, 1, 6, 4)

  xstar_i <- c(raw_counts_sum.G2ALL.S02C1.15m_A = 201.68944, raw_counts_sum.G2ALL.S02C1.15m_B = 127.3064219,
               raw_counts_sum.G2ALL.S02C1.15m_C = 92.152581, raw_counts_sum.G2ALL.S05C1.15m_A = 237.559910010157,
               raw_counts_sum.G2ALL.S05C1.15m_B = 158.1550348169, raw_counts_sum.G2ALL.S05C1.15m_C = 312.489719,
               raw_counts_sum.G2ALL.S06C1.15m_A = 880.352877, raw_counts_sum.G2ALL.S06C1.15m_B = 899.822609,
               raw_counts_sum.G2ALL.S06C1.15m_C = 1593.640589, raw_counts_sum.G2ALL.S07C1.15m_A = 1503.400459,
               raw_counts_sum.G2ALL.S07C1.15m_B = 1583.012547, raw_counts_sum.G2ALL.S07C1.15m_C = 2105.019572,
               raw_counts_sum.G2ALL.S09C1.15m_A = 600.5897, raw_counts_sum.G2ALL.S09C1.15m_B = 750.13492,
               raw_counts_sum.G2ALL.S09C1.15m_C = 674.74847586644, raw_counts_sum.G2ALL.S11C1.15m_A = 234.585574,
               raw_counts_sum.G2ALL.S11C1.15m_B = 312.5290412931, raw_counts_sum.G2ALL.S11C1.15m_C = 226.34653,
               raw_counts_sum.G2ALL.S15C1.15m_A = 195.221137077546, raw_counts_sum.G2ALL.S15C1.15m_B = 278.31226,
               raw_counts_sum.G2ALL.S15C1.15m_C = 233.5987703, raw_counts_sum.G2ALL.S16C1.15m_A = 70.57603,
               raw_counts_sum.G2ALL.S16C1.15m_B = 111.66915903, raw_counts_sum.G2ALL.S16C1.15m_C = 123.100246,
               raw_counts_sum.G2ALL.S17C1.15m_A = 205.701248, raw_counts_sum.G2ALL.S17C1.15m_B = 243.516225,
               raw_counts_sum.G2ALL.S17C1.15m_C = 217.50233, raw_counts_sum.G2ALL.S18C1.15m_A = 206.28941,
               raw_counts_sum.G2ALL.S18C1.15m_B = 320.42961, raw_counts_sum.G2ALL.S18C1.15m_C = 307.036039
  )

  xx_i <- xstar_i

  out_i_k <- fit_mgx_model(
    yy = yy_i_k,
    xstar = xstar_i,
    xx = xx_i,
    wts = NULL
  )
  out_i_k

  expect_type(out_i_k, "double")

})
