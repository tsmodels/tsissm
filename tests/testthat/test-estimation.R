test_that("constant variance benchmark", {
    cf_b <- c(coef(mod_constant_benchmark), "sigma" = sigma(mod_constant_benchmark), "LogLik" = as.numeric(logLik(mod_constant_benchmark)))
    expect_equal(constant_benchmark, cf_b, tolerance = 1e-5)
})

test_that("dynamic variance benchmark", {
    cf_b <- c(coef(mod_dynamic_benchmark), "LogLik" = as.numeric(logLik(mod_dynamic_benchmark)))
    expect_equal(dynamic_benchmark, cf_b, tolerance = 1e-5)
})

test_that("box-cox benchmark", {
    cf_b <- c(coef(mod_bc_benchmark), "sigma" = sigma(mod_bc_benchmark), "LogLik" = as.numeric(logLik(mod_bc_benchmark)))
    expect_equal(bc_benchmark, cf_b, tolerance = 1e-5)
})

test_that("vcov test (NW)", {
    v_nw <- vcov(mod_dynamic_benchmark, type = "NW", kernel = "Quadratic Spectral")
    expect_true(is.matrix(v_nw))
    expect_true(isSymmetric(v_nw))
})

test_that("dynamic variance benchmark (tmb vs eigen filter)", {
    spec <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, variance = "dynamic", seasonal_frequency = 12, seasonal_harmonics = 5)
    mod <- estimate(spec, debug_mode = TRUE)
    sig_tmb <- mod$tmb$report(coef(mod))$sigma
    sig_filt <- mod$model$sigma
    expect_equal(sig_tmb, sig_filt, tolerance = 1e-8)
})


test_that("constant variance benchmark (tmb vs eigen filter)", {
    spec <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, variance = "constant", seasonal_frequency = 12, seasonal_harmonics = 5)
    mod <- estimate(spec, debug_mode = TRUE)
    error_tmb <- mod$tmb$report(coef(mod))$error[-1]
    error_filt <- mod$model$error
    expect_equal(error_tmb, error_filt, tolerance = 1e-8)
})

    