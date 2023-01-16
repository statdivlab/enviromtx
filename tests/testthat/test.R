
n <- 20
set.seed(1)
fit_mgx_model(yy = rpois(n, lambda=1),
              xstar = rpois(n, lambda=1),
              xx = rpois(n, lambda=1),
              wts = rpois(n, lambda=1000))
