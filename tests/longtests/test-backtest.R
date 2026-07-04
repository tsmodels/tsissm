test_that("backtest", {
    data("spy", package = "tsissm")
    y <- as.xts(spy)
    spec_dynamic <- issm_modelspec(y["2014/"], slope = TRUE, seasonal = FALSE, ar = 1, ma = 0, xreg = NULL,
                                   lambda = 0, variance = "dynamic", distribution = "jsu")
    b <- tsbacktest(spec_dynamic, start = 2000, end = 2327, rolling = TRUE, h = 1, estimate_every = 30, trace = FALSE, seed = 100)
    expect_s3_class(b, "tsissm.backtest")
})

test_that("profile", {
    data("spy", package = "tsissm")
    y <- as.xts(spy)
    spec_dynamic <- issm_modelspec(y["2014/"], slope = TRUE, seasonal = FALSE, ar = 1, ma = 0,
                                   lambda = 0, variance = "dynamic", distribution = "jsu")
    mod <- estimate(spec_dynamic)
    p <- tsprofile(mod, h = 10, nsim = 100, seed = 100)
    expect_s3_class(p, "tsissm.profile")
})
    