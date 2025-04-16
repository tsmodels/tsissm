test_that("constant variance simulation (validation)", {
    z <- as.numeric(residuals(mod_constant_benchmark, standardize = TRUE, transformed = TRUE))
    sim <- simulate(mod_constant_benchmark, innov = matrix(z, nrow = 1), innov_type = "z", h = nrow(y), nsim = 1, sim_dates = index(y), init_states = mod_constant_benchmark$model$xseed, seed = 100)
    out <- cbind(as.numeric(y), as.numeric(sim$distribution))
    expect_equal(out[,1], out[,2], tolerance = 1e-8)
})

test_that("dynamic variance simulation (validation)", {
    z <- as.numeric(residuals(mod_dynamic_benchmark, standardize = TRUE, transformed = TRUE))
    sim <- simulate(mod_dynamic_benchmark, innov = matrix(z, nrow = 1), innov_type = "z", h = nrow(y), nsim = 1, sim_dates = index(y), 
                    init_states = mod_dynamic_benchmark$model$xseed, init_sigma = mod_dynamic_benchmark$model$initial_variance^0.5, 
                    init_res = mod_dynamic_benchmark$model$initial_variance^0.5,
                    seed = 100)
    out <- cbind(as.numeric(y), as.numeric(sim$distribution))
    expect_equal(out[,1], out[,2], tolerance = 1e-8)
})



test_that("constant variance simulation", {
    sim <- simulate(mod_constant_benchmark, h = 100, nsim = 2, seed = 100)
    expect_equal(NCOL(sim$distribution), 100L)
    expect_equal(NROW(sim$distribution), 2L)
    
})

test_that("dynamic variance simulation", {
    sim <- simulate(mod_dynamic_benchmark, h = 100, nsim = 2, seed = 100)
    expect_equal(NCOL(sim$distribution), 100L)
    expect_equal(NROW(sim$distribution), 2L)
})

