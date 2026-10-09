test_that("solver_status maps native nloptr codes to the unified status", {
    # nloptr: 1-4 success -> 0, 5-6 limit -> 1, <0 failure -> -1
    for (s in 1:4) {
        expect_identical(tsissm:::solver_status(list(status = s), "nloptr"), 0L)
    }
    for (s in 5:6) {
        expect_identical(tsissm:::solver_status(list(status = s), "nloptr"), 1L)
    }
    for (s in -(1:5)) {
        expect_identical(tsissm:::solver_status(list(status = s), "nloptr"), -1L)
    }
})
