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
