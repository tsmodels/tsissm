test_that("constant variance benchmark (filter)", {
    newy <- xts(551844 + 50, as.Date("2025-03-31"))
    f <- tsfilter(mod_constant_benchmark, y = newy)
    expect_equal(mod_constant_benchmark$model$fitted, f$model$fitted[-length(f$model$fitted)], tolerance = 1e-8)
})


test_that("dynamic variance benchmark (filter)", {
    newy <- xts(551844 + 50, as.Date("2025-03-31"))
    f <- tsfilter(mod_dynamic_benchmark, y = newy)
    expect_equal(mod_dynamic_benchmark$model$fitted, f$model$fitted[-length(f$model$fitted)], tolerance = 1e-8)
})


