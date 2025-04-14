test_that("single specification",{
    spec <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, seasonal_frequency = 12,
                           seasonal_harmonics = 5, ar = 1, ma = 1, xreg = NULL, 
                           variance = "constant", lambda = NA, distribution = "std")
    # Level(1), Slope(1), Seasonal (2), AR(1), MA(1), lambda(1) shape(1)
    esimation_parameters <- 1 + 1 + 2 + 1 + 1 + 1 + 1
    expect_equal(esimation_parameters, sum(spec$parmatrix$estimate))
})

test_that("Returns tsissm.spec object when auto is FALSE", {
    spec <- issm_modelspec(y, auto = FALSE)
    expect_true(inherits(spec, "tsissm.spec"))
})

test_that("Returns tsissm.autospec object when auto is TRUE", {
    spec <- issm_modelspec(y, auto = TRUE)
    expect_true(inherits(spec, "tsissm.autospec"))
})

test_that("Errors when y is not an xts object", {
    # Passing a numeric vector should trigger an error
    expect_error(issm_modelspec(as.numeric(1:10)))
})

test_that("Sets seasonal_frequency parameter correctly", {
    spec <- issm_modelspec(y, seasonal = TRUE, seasonal_frequency = c(7, 15), seasonal_harmonics = c(2,4))
    # Assuming the returned object stores seasonal_frequency as provided
    expect_equal(spec$seasonal$seasonal_frequency, c(7, 15))
    expect_equal(spec$seasonal$seasonal_harmonics, c(2, 4))
    
})

