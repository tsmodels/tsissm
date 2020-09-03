#################################################################
# from rugarch
.armaroots = function(coefs)
{
    Names <- names(coefs)
    ar = ma = NULL
    if (any(substr(Names, 1, 5) == "theta")) {
        ar <- which(substr(Names, 1, 5) == "theta")
        armap <- length(ar)
        arcoef <- coefs[ar]
        zar <- polyroot(c(1,-arcoef))
        rezar <- Re(zar)
        imzar <- Im(zar)
        nmr <- paste("ar", 1:armap, sep = "")
    } else {
        zar <- NULL
        rezar <- NULL
        imzar <- NULL
        nmr <- NULL
    }
    if (any(substr(Names, 1, 3) == "psi")) {
        ma <- which(substr(Names, 1, 3) == "psi")
        armaq <- length(ma)
        macoef <- coefs[ma]
        zma <- polyroot(c(1,macoef))
        rezma <- Re(zma)
        imzma <- Im(zma)
        nma <- paste("ma", 1:armaq, sep = "")
    } else {
        zma <- NULL
        rezma <- NULL
        imzma <- NULL
        nma <- NULL
    }
    root <- list()
    root$ar <- zar
    root$ma <- zma
    root$realar <- rezar
    root$imagar <- imzar
    root$realma <- rezma
    root$imagma <- imzma
    amp <- list()
    atan <- list()
    degree <- list()
    if (!is.null(zar)) {
        # Modulus
        amp$ar <- apply(cbind(root$realar, root$imagar), 1, FUN = function(x) sqrt(x[1]^2 + x[2]^2))
        atan$ar <- apply(cbind(amp$ar, root$realar),1 , FUN = function(x) atan2(x[1], x[2]))
        degree$ar <- atan$ar * 57.29577951
    }
    if (!is.null(zma)) {
        # Modulus
        amp$ma <- apply(cbind(root$realma, root$imagma), 1, FUN = function(x) sqrt(x[1]^2 + x[2]^2))
        atan$ma <- apply(cbind(amp$ma, root$realma),1 , FUN = function(x) atan2(x[1], x[2]) )
        degree$ma <- atan$ma * 57.29577951
    }
    res <- list(root = root, amp = amp, atan = atan, deg = degree)
    return(res)
}

.plotarmaroots = function(x, ...)
{
    arroot <- x$root$ar
    maroot <- x$root$ma
    if (!is.null(arroot)) {
        plot(1/arroot, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "", ylab = "", pch = 23, col = 2, ...)
        if (!is.null(maroot)) points(1/maroot, pch = 21, col = 3, ...)
        xy <- (2*pi/360)*(0:360)
        lines(sin(xy), cos(xy), col = "darkgrey")
        abline(h = 0, col = "darkgrey")
        abline(v = 0, col = "darkgrey")
        title("Inverse Roots and Unit Circle\n", xlab = "Real Part", ylab = "Imaginary Part")
        mtext(c("AR = Red | MA = Green"), side = 3, outer = F, cex = 0.6)
    } else {
        if (!is.null(maroot)) {
            plot(1/maroot, xlim = c(-1, 1), ylim = c(-1, 1), xlab = "", ylab = "", pch = 23, ...)
            xy <- (2*pi/360)*(0:360)
            lines(sin(xy), cos(xy), col = "darkgrey")
            abline(h = 0, col = "darkgrey")
            abline(v = 0, col = "darkgrey")
            title("Inverse Roots and Unit Circle\n", xlab = "Real Part", ylab = "Imaginary Part")
        }
    }
    invisible(x)
}

.pointsarmaroots = function(x, ...)
{
    arroot <- x$root$ar
    maroot <- x$root$ma
    if (!is.null(arroot)) points(1/arroot, pch = 23,  ...)
    if (!is.null(maroot)) points(1/maroot, pch = 21, ...)
    invisible(x)
}

##########################################################################################
# Direct Import of Weighted Tests of FISHER and GALLAGHER (WeightedPortTest package)
weighted_box_test = function(x, lag = 1, type = c("Box-Pierce", "Ljung-Box", "Monti"),
                              fitdf = 0, sqrd.res = FALSE, log.sqrd.res = FALSE, abs.res = FALSE,
                              weighted = TRUE)
{
    if (lag < (2 * fitdf + fitdf - 1)) stop("\nLag must be equal to a minimum of 2*fitdf+fitdf-1")
    if (NCOL(x) > 1) stop("\nx is not a vector or univariate time series")
    if (lag < 1) stop("\nLag must be positive")
    if (fitdf < 0) stop("\nFitdf cannot be negative")
    if ((sqrd.res && log.sqrd.res) || (sqrd.res && abs.res) || (log.sqrd.res && abs.res)) stop("Only one option of: sqrd.res, log.sqrd.res or abs.res can be selected")
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)
    if (abs.res) {
        x <- abs(x)
    }
    if (sqrd.res || log.sqrd.res) {
        x <- x^2
    }
    if (log.sqrd.res) {
        x <- log(x)
    }
    if (weighted) {
        if (type == "Monti") {
            METHOD <- "Weighted Monti test (Gamma Approximation)"
            cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE,
                       na.action = na.pass)
            obs <- cor$acf[1:lag]
        }
        else {
            cor <- acf(x, lag.max = lag, type = "correlation",
                       plot = FALSE, na.action = na.pass)
            obs <- cor$acf[2:(lag + 1)]
        }
        if (type == "Ljung-Box") {
            METHOD <- "Weighted Ljung-Box test (Gamma Approximation)"
        }
        n <- sum(!is.na(x))
        weights <- (lag - 1:lag + 1)/(lag)
        if (type == "Box-Pierce") {
            METHOD <- "Weighted Box-Pierce test (Gamma Approximation)"
            STATISTIC <- n * sum(weights * obs^2)
        }
        else {
            STATISTIC <- n * (n + 2) * sum(weights * (1/seq.int(n - 1, n - lag) * obs^2))
        }
        if (sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "Weighted X-squared on Squared Residuals for detecting nonlinear processes"
        }
        else if (log.sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "Weighted X-squared on Log-Squared Residuals for detecting nonlinear processes"
        }
        else if (abs.res) {
            fitdf <- 0
            names(STATISTIC) <- "Weighted X-squared on Absolute valued Residuals for detecting nonlinear processes"
        }
        else {
            names(STATISTIC) <- "Weighted X-squared on Residuals for fitted ARMA process"
        }
        shape <- (3/4) * (lag + 1)^2 * lag/(2 * lag^2 + 3 * lag + 1 - 6 * lag * fitdf)
        scale <- (2/3) * (2 * lag^2 + 3*lag + 1 - 6 * lag * fitdf)/(lag*(lag + 1))
        PARAMETER <- c(shape, scale)
        names(PARAMETER) <- c("Shape", "Scale")
        PVAL <- 1 - pgamma(STATISTIC, shape = shape, scale = scale)
        names(PVAL) <- "Approximate p-value"
    }
    else {
        if (type == "Monti") {
            METHOD <- "Monti test"
            cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE,
                       na.action = na.pass)
            obs <- cor$acf[1:lag]
        }
        else {
            cor <- acf(x, lag.max = lag, type = "correlation",
                       plot = FALSE, na.action = na.pass)
            obs <- cor$acf[2:(lag + 1)]
        }
        if (type == "Ljung-Box") {
            METHOD <- "Ljung-Box test"
        }
        n <- sum(!is.na(x))
        if (type == "Box-Pierce") {
            METHOD <- "Box-Pierce test"
            STATISTIC <- n * sum(obs^2)
        }
        else {
            STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1,
                                                      n - lag) * obs^2))
        }
        if (sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "X-squared on Squared Residuals for detecting nonlinear processes"
        }
        else if (log.sqrd.res) {
            fitdf <- 0
            names(STATISTIC) <- "X-squared on Log-Squared Residuals for detecting nonlinear processes"
        }
        else if (abs.res) {
            fitdf <- 0
            names(STATISTIC) <- "X-squared on Absolute valued Residuals for detecting nonlinear processes"
        }
        else {
            names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
        }
        mydf <- lag - fitdf
        PARAMETER <- c(mydf)
        names(PARAMETER) <- c("df")
        PVAL <- 1 - pchisq(STATISTIC, df = mydf)
        names(PVAL) <- "p-value"
    }
    structure(list(statistic = STATISTIC, parameter = PARAMETER,
                   p.value = PVAL, method = METHOD, data.name = DNAME),
              class = "htest")
}
