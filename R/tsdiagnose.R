#' Model Diagnostics
#'
#' @description Creates a short summary of model based diagnostics.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param plot whether to generate diagnostic plots to accompany summary.
#' @param ... not currently used.
#' @returns A list of tests including the weighted Ljung-Box test for residual 
#' autocorrelation, system forecastability test, and outlier dates using the 
#' Rosner test (\code{\link[EnvStats]{rosnerTest}}). Optionally generates a plot
#' of relevant diagnostics.
#' @aliases tsdiagnose
#' @method tsdiagnose tsissm.estimate
#' @rdname tsdiagnose
#' @export
#'
#'
tsdiagnose.tsissm.estimate <- function(object, plot = FALSE, ...)
{
    if (sum(object$spec$arma$order) > 0) {
        armav <- coef(object)
        armav <- armav[grepl("theta|psi", names(armav))]
        rt <- .armaroots(armav)
    } else {
        rt <- NULL
    }
    e <- abs(eigen(object$model$D, symmetric = FALSE)$values)
    df <- sum(object$spec$arma$order)
    sigma <- sigma(object)
    r <- as.numeric(na.omit(residuals(object, transformed = TRUE)/sigma))
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
    rtest <- .rosner_test(as.numeric(na.omit(residuals(object, transformed = TRUE))), k = 10)
    if (any(rtest$Outlier)) {
        outliers_index <- object$spec$target$index[which(object$spec$good == 1)][rtest$Obs.Num[rtest$Outlier]]
    } else {
        outliers_index <- NULL
    }
    if (plot) {
        opar <- par(no.readonly = TRUE)
        on.exit(par(opar))
        par(mfrow = c(2,2), mar = c(3,3,3,3))
        if (df > 0) .plotarmaroots(.armaroots(armav))
        acf(as.numeric(r), type = "correlation", main = "Scaled Residuals Autocorrelation")
        hist(r, breaks = "fd", main = "Scaled Residuals Histogram", probability = TRUE)
        box()
        qqnorm(r)
        qqline(r, col = 2)
    }
    L <- list(arma_test = rt, 
              stability_test = e, 
              weighted_box_test = lbsr, 
              rosner_test = rtest, 
              outliers_index = outliers_index)
    class(L) <- "tsissm.diagnose"
    return(L)
}


#' Model Diagnostics Print method
#'
#' @description Print method for class \dQuote{tsissm.diagnose}
#' @param x an object of class \dQuote{tsissm.duagnose} generated from
#' calling \code{\link[tsissm]{tsdiagnose}}.
#' @param ... not currently used.
#' @returns Invisibly returns the original object and prints the output to console.
#' @aliases print.tsissm.diagnose
#' @method print tsissm.diagnose
#' @rdname print.tsissm.diagnose
#' @export
#'
#'
print.tsissm.diagnose <- function(x, ...)
{
    arma_test <- x$arma_test
    if (!is.null(arma_test)) {
        cat("\nARMA roots (<1)")
        cat("\n------------------------------------------\n")
        if (!is.null(arma_test$root$ar)) {
            cat("Inverse AR roots:", 1/abs(arma_test$root$ar))
            cat("\n")
        }
        if (!is.null(arma_test$root$ma)) {
            cat("Inverse MA roots:", 1/abs(arma_test$root$ma))
            cat("\n")
        }
    }
    cat("\nForecastability")
    cat("\n------------------------------------------\n")
    cat("Real Eigenvalues (D):", round(x$stability_test,3))
    cat("\n")
    cat("\nWeighted Ljung-Box Test [scaled residuals]")
    cat("\n------------------------------------------\n")
    print(x$weighted_box_test, row.names = FALSE, digits = 3)
    if (any(x$rosner_test$Outlier)) {
        cat("\nOutlier Diagnostics (based on Rosner Test)")
        cat("\n------------------------------------------")
        cat("\nOutliers:", as.character(x$outliers_index))
    }
    return(invisible(x))
}
