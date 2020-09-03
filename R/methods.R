fitted.tsissm.estimate <- function(object, ...)
{
    return(xts(object$model$fitted, object$spec$target$index))
}

residuals.tsissm.estimate <- function(object, raw = FALSE, h = 1, cores = 1, seed = NULL, trace = FALSE, index_start = 1, ...)
{
    if (h > 1) {
        out <- hresiduals.tsissm.estimate(object, h = h, cores = cores, seed = seed, trace = trace, raw = raw, index_start = 1, ...)
    } else {
        parameters <- NULL
        if (!raw) {
            out <- xts(object$model$residuals, object$spec$target$index)
        } else {
            f <- object$spec$transform$transform(object$model$fitted, lambda = object$parmatrix[parameters == "lambda"]$optimal)
            out <- xts(object$spec$target$y - f, object$spec$target$index)
        }
    }
    return(out)
}


summary.tsissm.estimate <- function(object, digits = 4, ...)
{
    estimate <- NULL
    cf <- object$parmatrix[estimate == 1, c("parameters","optimal","lower","upper")]
    cf$optimal <- round(cf$optimal, digits + 2)
    cf$lower <- round(cf$lower,2)
    cat("\n-------------------------------------")
    cat("\ntsissm:      Summary                 ")
    cat("\n-------------------------------------\n")
    print(cf, row.names = FALSE, digits = digits)
    cat("\n-------------------------------------")
    tsm <- tsmetrics(object)
    return(invisible(tsm))
}

.resid.variance <- function(object)
{
    parameters <- NULL
    td <- tsdecompose(object)
    snames <- colnames(td)
    n <- length(snames)
    v <- rep(0, ncol(td))
    r <- xts(matrix(0, ncol = ncol(td), nrow = nrow(td)), object$spec$target$index)
    A <- xts(object$spec$transform$transform(object$spec$target$y_orig, lambda = object$parmatrix[parameters == "lambda"]$optimal), object$spec$target$index)
    r[,1] = A - td[,1]
    v[1] <- var(r[,1])
    if (n > 1) {
        for (i in 2:n) {
            v[i] <- var(r[,i - 1] - td[,i])
        }
    }
    names(v) <- c("V[Actual - Level]", paste0("V[...-",snames[-1],"])"))
    attr(v, "states") <- snames
    return(v)
}
coef.tsissm.estimate <- function(object, ...)
{
    estimate <- NULL
    return(structure(object$parmatrix[estimate == 1]$optimal, names = object$parmatrix[estimate == 1]$parameters))
}

tsdecompose.tsissm.estimate <- function(object, ...)
{
    estimate <- NULL
    idx <- object$spec$idmatrix
    rnames <- rownames(idx)
    states <- rbind(object$model$xseed, object$model$states)
    indx <- object$spec$target$index
    pars <- object$parmatrix[estimate == 1]$optimal
    names(pars) <- object$parmatrix[estimate == 1]$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    mdim = object$spec$dims
    w <- matrix(S[list("w")]$values, nrow = mdim[1], ncol = 1)
    w_t <- t(w)
    if (idx["Level","n"] > 0) {
        Level <- xts(states[-nrow(states), idx["Level","start"]:idx["Level","end"]], indx)
    } else {
        Level <- NULL
    }
    if (idx["Slope","n"] > 0) {
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- xts(w_t[,nstart:nend] * states[-nrow(states), idx["Slope","start"]:idx["Slope","end"]], indx)
        colnames(Slope) <- "Slope"
    } else {
        Slope <- NULL
    }
    if (any(idx[grepl("Seasonal",rnames),"n"] > 0)) {
        ns <- length(idx[grepl("Seasonal",rnames),"n"])
        nstart <- idx[grepl("Seasonal",rnames),"start"]
        nend <- idx[grepl("Seasonal",rnames),"end"]
        frequency <- object$spec$seasonal$seasonal_frequency
        Seasonal <- matrix(0, ncol = ns, nrow = length(indx))
        colnames(Seasonal) <- names(nstart)
        if (object$spec$seasonal$seasonal_type == "trigonometric") {
            K <- object$spec$seasonal$seasonal_harmonics
            for (i in 1:ns) {
                Seasonal[,i] <- t(w_t[,nstart[i]:(nstart[i] + K[i] - 1)] %*% t(states[-nrow(states), nstart[i]:(nstart[i] + K[i] - 1)]))
            }
        } else {
            for (i in 1:ns) {
                Seasonal[,i] <- t(w_t[,nstart[i]:(nstart[i] + frequency[i] - 1)] %*% t(states[-nrow(states), nstart[i]:(nstart[i] + frequency[i] - 1)]))
            }
        }
        Seasonal <- xts(Seasonal, indx)
    } else {
        Seasonal <- NULL
    }
    if (idx["AR","n"] > 0) {
        AR <- t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[-nrow(states), idx["AR","start"]:idx["AR","end"]]))
        AR <- xts(AR, indx)
        colnames(AR) <- paste0("AR",idx["AR","n"])
    } else {
        AR <- NULL
    }
    if (idx["MA","n"] > 0) {
        MA <- t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[-nrow(states), idx["MA","start"]:idx["MA","end"]]))
        MA <- xts(MA, indx)
        colnames(MA) <- paste0("MA",idx["MA","n"])
    } else {
        MA <- NULL
    }
    if (idx["X","n"] > 0) {
        beta <- matrix(object$parmatrix[which(grepl("kappa", object$parmatrix$parameters))]$optimal, ncol = 1)
        X <- object$spec$xreg$xreg
        xreg <- xts(X %*% beta, indx)
        colnames(xreg) <- "X"
    } else {
        xreg <- NULL
    }
    S <- cbind(Level, Slope, Seasonal, AR, MA, xreg)
    return(S)
}


tsdecompose.tsissm.predict <- function(object, ...)
{
    estimate <- NULL
    idx <- object$spec$idmatrix
    rnames <- rownames(idx)
    states <- object$states
    indx <- object$spec$target$index
    pars <- object$spec$parmatrix[estimate == 1]$optimal
    names(pars) <- object$spec$parmatrix[estimate == 1]$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    mdim <- object$spec$dims
    nsim <- nrow(object$distribution)
    fdates <- colnames(object$distribution)
    date_class <- attr(object$distribution, "date_class")
    L <- list()
    k <- 1
    w <- matrix(S[list("w")]$values, nrow = mdim[1], ncol = 1)
    w_t <- t(w)
    if (idx["Level","n"] > 0) {
        Level <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(states[-1, idx["Level","start"]:idx["Level","end"], i],nrow = 1)
        }))
        colnames(Level) <- fdates
        class(Level) <- "tsmodel.distribution"
        attr(Level, "date_class") <- date_class
        Level <- list(original_series = object$decomp$Level, distribution = Level)
        class(Level) <- "tsmodel.predict"
        L[[1]] <- Level
        k <- k + 1
    }
    if (idx["Slope","n"] > 0) {
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(w_t[,nstart:nend] * states[-1, idx["Slope","start"]:idx["Slope","end"], i], nrow = 1)
        }))
        colnames(Slope) <- fdates
        class(Slope) <- "tsmodel.distribution"
        attr(Slope, "date_class") <- date_class
        Slope <- list(original_series = object$decomp$Slope, distribution = Slope)
        class(Slope) <- "tsmodel.predict"
        L[[k]] <- Slope
        k <- k + 1
    }

    if (any(idx[grepl("Seasonal",rnames),"n"] > 0)) {
        ns <- length(idx[grepl("Seasonal",rnames),"n"])
        nstart <- idx[grepl("Seasonal",rnames),"start"]
        nend <- idx[grepl("Seasonal",rnames),"end"]
        frequency <- object$spec$seasonal$seasonal_frequency
        if (object$spec$seasonal$seasonal_type == "trigonometric") {
            K <- object$spec$seasonal$seasonal_harmonics
            for (j in 1:ns) {
                tmp <- do.call(rbind, lapply(1:nsim, function(i){
                    matrix(t(w_t[,nstart[j]:(nstart[j] + K[j] - 1)] %*% t(states[-1, nstart[j]:(nstart[j] + K[j] - 1), i])), nrow = 1)
                }))
                colnames(tmp) <- fdates
                class(tmp) <- "tsmodel.distribution"
                attr(tmp, "date_class") <- date_class
                tmp <- list(original_series = object$decomp[,k], distribution = tmp)
                class(tmp) <- "tsmodel.predict"
                L[[k]] <- tmp
                k <- k + 1
            }
        } else {
            for (i in 1:ns) {
                tmp <- do.call(rbind, lapply(1:nsim, function(i){
                    matrix(t(w_t[,nstart[j]:(nstart[j] + frequency[j] - 1)] %*% t(states[-1, nstart[j]:(nstart[j] + frequency[j] - 1), i])), nrow = 1)
                }))
                colnames(tmp) <- fdates
                class(tmp) <- "tsmodel.distribution"
                attr(tmp, "date_class") <- date_class
                tmp <- list(original_series = object$decomp[,k], distribution = tmp)
                class(tmp) <- "tsmodel.predict"
                L[[k]] <- tmp
                k <- k + 1
            }
        }
    }
    if (idx["AR","n"] > 0) {
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[-1, idx["AR","start"]:idx["AR","end"], i])), nrow = 1)
        }))
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        tmp <- list(original_series = object$decomp[,k], distribution = tmp)
        class(tmp) <- "tsmodel.predict"
        L[[k]] <- tmp
        k <- k + 1
    }
    if (idx["MA","n"] > 0) {
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[-1, idx["MA","start"]:idx["MA","end"], i])), nrow = 1)
        }))
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        tmp <- list(original_series = object$decomp[,k], distribution = tmp)
        class(tmp) <- "tsmodel.predict"
        L[[k]] <- tmp
        k <- k + 1
    }
    names(L) <- colnames(object$decomp[,1:length(L)])
    return(L)
}

logLik.tsissm.estimate <- function(object, ...)
{
    parameters <- NULL
    estimate <- NULL
    np <- NROW(object$parmatrix[estimate == 1]) + ncol(object$model$xseed)
    r <- residuals(object, scaled = TRUE)
    v <- sum(r^2)
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    n <- NROW(object$spec$target$y_orig)
    loglik <- -0.5 * n * log(2 * pi * (v/n)) - (1/(2 * (v/n))) * v + (lambda - 1) * sum(log(object$spec$target$y_orig))
    return(structure(loglik, df = np + 1, class = "logLik"))
}

AIC.tsissm.estimate <- function(object, ..., k = 2)
{
    estimate <- NULL
    nr <- NROW(object$parmatrix[estimate == 1]) + NCOL(object$model$xseed)
    return(object$model$loglik + k * nr)
}

tsmetrics.tsissm.estimate <- function(object, ...)
{
    parameters <- NULL
    estimate <- NULL
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    nr <- NROW(object$parmatrix[estimate == 1]) + NCOL(object$model$xseed)
    AIC <- object$model$loglik + 2 * nr
    MAPE <- mape(object$spec$target$y_orig, fitted(object))
    BIAS <- bias(object$spec$target$y_orig, fitted(object))
    MSLRE <- mslre(object$spec$target$y_orig, fitted(object))
    yt <- object$spec$transform$transform(object$spec$target$y_orig, lambda = lambda)
    ft <- object$spec$transform$transform(as.numeric(fitted(object)), lambda = lambda)
    r <- yt - ft
    cat("\ntsissm: Performance Metrics")
    cat("\n----------------------------------\n")
    cat(paste0("AIC\t: ", round(AIC,2), " (n = ", nr,")"))
    cat("\n")
    cat(paste0("MAPE\t: ", round(MAPE,5)))
    cat("\n")
    cat(paste0("BIAS\t: ", round(BIAS,5)))
    cat("\n")
    cat(paste0("MSLRE\t: ", round(MSLRE, 5)))
    metrics = c(AIC, MAPE, BIAS, MSLRE)
    names(metrics) <- c("AIC","MAPE","BIAS","MSLRE")
    return(invisible(metrics))
}

tsmetrics.tsissm.predict = function(object, actual, alpha = 0.1, ...)
{
    n <- NCOL(object$distribution)
    if (NROW(actual) != n) stop("\nactual length not equal to forecast length")
    m_mape <- mape(actual, colMeans(object$distribution))
    m_bias <- bias(actual, colMeans(object$distribution))
    m_mslre <- mslre(actual, colMeans(object$distribution))
    m_mase <- mase(actual, colMeans(object$distribution), object$original_series, frequency = object$spec$target$frequency[1])
    if (!is.null(alpha)) {
        m_mis <- mis(actual, lower = apply(object$distribution, 2, quantile, alpha/2), upper = apply(object$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
    } else {
        m_mis <- as.numeric(NA)
    }
    return(data.frame("h" = n, "MAPE" = m_mape, "MASE" = m_mase, "MSLRE" = m_mslre, "BIAS" = m_bias, "MIS" = m_mis))
}
