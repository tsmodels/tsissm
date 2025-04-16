test_that("ensemble", {
    p <- predict(mod_auto, h = 20, nsim = 100)
    p_ens <- tsensemble(p, weights = c(0.5, 0.5), start = 0)
    tsd <- p_ens$decomposition$Level$distribution + p_ens$decomposition$Slope$distribution + 
        p_ens$decomposition$Seasonal12$distribution + p_ens$decomposition$Irregular$distribution
    tsd <- mod_auto$models[[1]]$spec$transform$inverse(tsd, lambda = 1)
    tsd_mean <- as.numeric(apply(tsd, 2, mean))
    p_mean <- as.numeric(apply(p_ens$distribution, 2, mean))
    expect_equal(tsd_mean, p_mean, tolerance = 1e-9)
})
