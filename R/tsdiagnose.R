tsdiagnose.tsissm.estimate <- function(object, plot = FALSE, ...)
{
    if (sum(object$spec$arma$order) > 0) {
        cat("\nARMA roots")
        cat("\n------------------------------------------\n")
        armav <- coef(object)
        armav <- armav[grepl("theta|psi", names(armav))]
        rt <- .armaroots(armav)
        if (object$spec$arma$order[1] > 0) {
            cat("Inverse AR roots:", 1/abs(rt$root$ar))
            cat("\n")
        }
        if (object$spec$arma$order[2] > 0) {
            cat("Inverse MA roots:", 1/abs(rt$root$ma))
            cat("\n")
        }
    } else {
        rt <- NULL
    }
    cat("\nForecastability (D roots)")
    cat("\n------------------------------------------\n")
    e <- abs(Re(eigen(object$model$D)$values))
    cat("Real Eigenvalues (D):", round(e,3))
    cat("\n")
    cat("\nWeighted Ljung-Box Test [scaled residuals]")
    cat("\n------------------------------------------\n")
    df <- sum(object$spec$arma$order)
    r <- residuals(object, scaled = TRUE)
    if (sum(object$spec$arma$order) > 0 ) j = 0 else j = 1
    b1 <- weighted_box_test(r, lag = 1, type = "Ljung-Box", fitdf = 0)
    b2j <- pmax(2 * df + df - 1, 1 + df + j)
    b2 <- weighted_box_test(r, lag = b2j, type = "Ljung-Box", fitdf = df)
    b3j <- pmax(2 * df + df - 1, 2 + df + j)
    b3 <- weighted_box_test(r, lag = b3j, type = "Ljung-Box", fitdf = df)
    b4j <- pmax(2 * df + df - 1, 3 + df + j)
    b4 <- weighted_box_test(r, lag = b4j, type = "Ljung-Box", fitdf = df)
    lbsr <- data.table(Lag =  c("Lag[1]", paste0("Lag[",b2j,"]"), paste0("Lag[",b3j,"]"), paste0("Lag[",b4j,"]")),
                       statistic = c(b1$statistic[[1]], b2$statistic[[1]], b3$statistic[[1]],b4$statistic[[1]]),
                       pvalue = c(b1$p.value[[1]], b2$p.value[[1]],b3$p.value[[1]], b4$p.value[[1]]))
    print(lbsr, row.names = FALSE, digits = 3)
    rtest <- rosnerTest(as.numeric(residuals(object)), k = 10)
    if (any(rtest$all.stats$Outlier)) {
        out.index <- object$spec$target$index[rtest$all.stats$Obs.Num[rtest$all.stats$Outlier]]
        cat("\nOutlier Diagnostics (based on Rosner Test)")
        cat("\n------------------------------------------")
        cat("\nOutliers:", as.character(out.index))
    } else {
        out.index <- NULL
    }
    if (plot) {
        par(mfrow = c(2,2), mar = c(3,3,3,3))
        if (df > 0) .plotarmaroots(.armaroots(armav))
        acf(as.numeric(r), type = "correlation", main = "Scaled Residuals Autocorrelation")
        hist(r, breaks = "fd", main = "Scaled Residuals Histogram", probability = T)
        box()
        qqnorm(r)
        qqline(r, col = 2)
    }
    L <- list(armaroots = rt, D.eigenvalues = e, lb_test = lbsr, outliers = rtest$all.stats)
    return(invisible(L))
}
