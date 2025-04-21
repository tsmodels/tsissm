test_that("hresiduals", {
    hr <- hresiduals(mod_dynamic_benchmark, h = 3, transformed = TRUE, index_start = 0)
    rr <- residuals(mod_dynamic_benchmark, transformed = TRUE)
    z1 <- as.numeric(hr$`1`)
    z2 <- as.numeric(rr)
    expect_equal(z1, z2, tolerance = 1e-8)
})


