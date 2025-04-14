test_that("constant variance prediction (h=1)", {
    p <- predict(mod_constant_benchmark, h = 1, nsim = 100, seed = 100)
    expect_equal(NROW(p$distribution), 100L)
    expect_equal(NCOL(p$distribution), 1L)
    expect_s3_class(p$analytic_mean, "zoo")
})


test_that("constant variance prediction (h>1)", {
    p <- predict(mod_constant_benchmark, h = 10, nsim = 100, seed = 100)
    expect_equal(NROW(p$distribution), 100L)
    expect_equal(NCOL(p$distribution), 10L)
    expect_s3_class(p$analytic_mean, "zoo")
    
    tsx <- tsdecompose(p, start = 0)
    aggregated_components <- tsx$Level$distribution + tsx$Slope$distribution + tsx$Seasonal12$distribution + tsx$Irregular$distribution
    p_unt <- matrix(mod_constant_benchmark$spec$transform$transform(as.numeric(p$distribution), lambda = mod_constant_benchmark$parmatrix[parameters == "lambda"]$value), ncol = 10, byrow = F)
    b1 <- as.numeric(apply(aggregated_components, 2, mean))
    b2 <- as.numeric(apply(p_unt, 2, mean))
    expect_equal(b1, b2, tolerance = 1e-9)
})


test_that("dynamic variance prediction (h=1)", {
    p <- predict(mod_dynamic_benchmark, h = 1, nsim = 100, seed = 100)
    expect_equal(NROW(p$distribution), 100L)
    expect_equal(NCOL(p$distribution), 1L)
    expect_s3_class(p$analytic_mean, "zoo")
})


test_that("dynamic variance prediction (h>1)", {
    p <- predict(mod_dynamic_benchmark, h = 10, nsim = 100, seed = 100)
    expect_equal(NROW(p$distribution), 100L)
    expect_equal(NCOL(p$distribution), 10L)
    expect_s3_class(p$analytic_mean, "zoo")
    
    tsx <- tsdecompose(p, start = 0)
    aggregated_components <- tsx$Level$distribution + tsx$Slope$distribution + tsx$Seasonal12$distribution + tsx$Irregular$distribution
    p_unt <- matrix(mod_dynamic_benchmark$spec$transform$transform(as.numeric(p$distribution), lambda = mod_dynamic_benchmark$parmatrix[parameters == "lambda"]$value), ncol = 10, byrow = F)
    b1 <- as.numeric(apply(aggregated_components, 2, mean))
    b2 <- as.numeric(apply(p_unt, 2, mean))
    expect_equal(b1, b2, tolerance = 1e-9)
    
})




test_that("bc  prediction (h=1)", {
    p <- predict(mod_bc_benchmark, h = 1, nsim = 100, seed = 100)
    expect_equal(NROW(p$distribution), 100L)
    expect_equal(NCOL(p$distribution), 1L)
    expect_s3_class(p$analytic_mean, "zoo")
})


test_that("bc  prediction (h>1)", {
    p <- predict(mod_bc_benchmark, h = 10, nsim = 100, seed = 100)
    expect_equal(NROW(p$distribution), 100L)
    expect_equal(NCOL(p$distribution), 10L)
    expect_s3_class(p$analytic_mean, "zoo")
})


test_that("bc  prediction taylor approximation", {
    p <- predict(mod_bc_benchmark, h = 10, nsim = 100, seed = 100)
    tsx <- tsdecompose(p, start = 0)
    aggregated_components <- tsx$Level$distribution + tsx$Slope$distribution + tsx$Seasonal12$distribution + tsx$Irregular$distribution
    p_unt <- matrix(mod_bc_benchmark$spec$transform$transform(as.numeric(p$distribution), lambda = mod_bc_benchmark$parmatrix[parameters == "lambda"]$value), ncol = 10, byrow = F)
    b1 <- as.numeric(apply(aggregated_components, 2, mean))
    b2 <- as.numeric(apply(p_unt, 2, mean))
    expect_equal(b1, b2, tolerance = 1e-9)
})

