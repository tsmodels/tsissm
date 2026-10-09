test_that("tsbacktest: simulation draws are seeded and reproducible", {
    spec <- issm_modelspec(y, slope = TRUE, seasonal = TRUE, seasonal_frequency = 12, seasonal_harmonics = 5)
    oplan <- future::plan("sequential")
    on.exit(future::plan(oplan), add = TRUE)
    warns <- character(0)
    set.seed(42)
    withCallingHandlers(
        b1 <- tsbacktest(spec, start = NROW(y) - 3, end = NROW(y), h = 1, rolling = FALSE, trace = FALSE),
        warning = function(w) {
            warns <<- c(warns, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    set.seed(42)
    b2 <- tsbacktest(spec, start = NROW(y) - 3, end = NROW(y), h = 1, rolling = FALSE, trace = FALSE)
    expect_false(any(grepl("UNRELIABLE VALUE", warns)))
    expect_identical(b1$table, b2$table)
})
