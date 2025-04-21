test_that("backtest", {
    data("spy", package = "tsissm")
    y <- as.xts(spy)
    xreg <- auto_regressors(y["2014/"], frequency = 1, lambda = 0, sampling = "days", method = "full",
                            check.rank = TRUE, discard.cval = 3.5, maxit.iloop = 10, maxit.oloop = 10, types = "LS")
    exc <- which(xreg$xreg["2020-02-03"] == 1)
    xreg$xreg <- xreg$xreg[,-exc]
    xreg$init <- xreg$init[-exc]
    spec_dynamic <- issm_modelspec(y["2014/"], slope = TRUE, seasonal = FALSE, ar = 1, ma = 0, xreg = xreg$xreg,
                                   lambda = 0, variance = "dynamic", distribution = "jsu")
    b <- tsbacktest(spec_dynamic, start = 2000, end = 2327, rolling = TRUE, h = 1, estimate_every = 30, trace = FALSE, seed = 100)
    expect_s3_class(b, "tsissm.backtest")
})

test_that("profile", {
    data("spy", package = "tsissm")
    y <- as.xts(spy)
    xreg <- auto_regressors(y["2014/"], frequency = 1, lambda = 0, sampling = "days", method = "full",
                            check.rank = TRUE, discard.cval = 3.5, maxit.iloop = 10, maxit.oloop = 10, types = "LS")
    exc <- which(xreg$xreg["2020-02-03"] == 1)
    xreg$xreg <- xreg$xreg[,-exc]
    xreg$init <- xreg$init[-exc]
    spec_dynamic <- issm_modelspec(y["2014/"], slope = TRUE, seasonal = FALSE, ar = 1, ma = 0, xreg = xreg$xreg,
                                   lambda = 0, variance = "dynamic", distribution = "jsu")
    mod <- estimate(spec_dynamic)
    p <- tsprofile(mod, h = 10, nsim = 100, seed = 100)
    expect_s3_class(p, "tsissm.profile")
})
    