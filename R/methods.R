fitted.tsissm.estimate <- function(object, ...)
{
    return(xts(object$model$fitted, object$spec$target$index))
}

residuals.tsissm.estimate <- function(object, raw = FALSE, h = 1, seed = NULL, trace = FALSE, index_start = 1, ...)
{
    if (h > 1) {
        out <- hresiduals.tsissm.estimate(object, h = h, seed = seed, trace = trace, raw = raw, index_start = 1, ...)
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



.make_standard_errors <- function(object)
{
    pmatrix <- object$parmatrix
    pars <- pmatrix[estimate == 1]$optimal
    H <- object$hessian
    se <- suppressWarnings(sqrt(diag(solve(H))))
    tvalues <- pars/se
    pvalues <- 2*(1 - pnorm(abs(tvalues)))
    return(data.frame("Std. Error" = se,"t value" = tvalues, "Pr(>|t|)" = pvalues, check.names = FALSE))
}

.make_model_description <- function(object) {
    model <- "A"
    if (object$spec$slope$include_slope) {
        model <- c(model,"A")
            if (object$spec$slope$include_damped) {
                model <- c(model,"d")                
            }
    } else {
        model = c(model, "N")
    }
    if (object$spec$seasonal$include_seasonal) {
        model <- c(model, "A")
        if (object$spec$seasonal$seasonal_type == "trigonometric") {
            s <- sapply(1:length(object$spec$seasonal$seasonal_frequency), function(i) {
                if (i < length(object$spec$seasonal$seasonal_frequency)) {
                    paste0(object$spec$seasonal$seasonal_frequency[i],"{",object$spec$seasonal$seasonal_harmonics[i],"}/")
                } else {
                    paste0(object$spec$seasonal$seasonal_frequency[i],"{",object$spec$seasonal$seasonal_harmonics[i],"}")
                }
            })
        } else {
            s <- sapply(1:length(object$spec$seasonal$seasonal_frequency), function(i) {
                paste0(object$spec$seasonal$seasonal_frequency[i],"{}")
            })
        }
        model <- c(model,"(",s,")")
    } else {
        model <- c(model,"N")
    }
    model <- c(model,"+ARMA(",object$spec$arma$order[1],",",object$spec$arma$order[2],")")
    if (object$spec$xreg$include_xreg) {
        model <- c(model,"+X(",NCOL(object$spec$xreg$xreg),")")
    }
    model <- paste0(model,collapse = "")
    return(model)
}

summary.tsissm.estimate <- function(object, digits = 4, ...)
{
    parameters <- NULL
    estimate <- NULL
    model <- .make_model_description(object)
    if (object$autodiff) {
        printout <- object$parmatrix[estimate == 1, c("parameters","optimal","lower","upper")]
        colnames(printout) <- c("Parameter","Est[Value]","Lower","Upper")
        printout <- as.data.frame(printout)
        S <- try(suppressWarnings(.make_standard_errors(object)), silent = TRUE)
        use_se <- FALSE
        if (!inherits(S,'try-error')) {
            printout <- cbind(printout, S)
            use_se <- TRUE
        }
        cat("ISSM Model:",model)
        if (object$spec$transform$name == "box-cox") {
            if (object$parmatrix[parameters == "lambda"]$estimate == 0) {
                lambda_df <- as.data.frame(object$parmatrix[parameters == "lambda",c("parameters","optimal","lower","upper")])
                colnames(lambda_df) <- c("Parameter","Est[Value]","Lower","Upper")
                if (use_se) {
                    lambda_df <- cbind(lambda_df, data.frame("Std. Error" = as.numeric(NaN),"t value" = as.numeric(NaN),"Pr(>|t|)" = as.numeric(NaN), check.names = FALSE))
                }
                printout <- rbind(printout, lambda_df)
            }
        }
        print(kable(printout, right = FALSE, digits = digits, row.names = FALSE, format = "simple"))
        tsm <- tsmetrics(object)
    } else {
        printout <- object$parmatrix[estimate == 1, c("parameters","optimal","lower","upper")]
        colnames(printout) <- c("Parameter","Est[Value]","Lower","Upper")
        printout <- as.data.frame(printout)
        if (object$spec$transform$name == "box-cox") {
            if (object$parmatrix[parameters == "lambda"]$estimate == 0) {
                lambda_df <- as.data.frame(object$parmatrix[parameters == "lambda",c("parameters","optimal","lower","upper")])
                colnames(lambda_df) <- c("Parameter","Est[Value]","Lower","Upper")
                printout <- rbind(printout, lambda_df)
            }
        }
        cat("ISSM Model:",model)
        print(kable(printout, right = FALSE, digits = digits, row.names = FALSE, format = "simple"))
        tsm <- tsmetrics(object)
    }
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

tsdecompose.tsissm.estimate <- function(object, simplify = FALSE, ...)
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
    if (simplify) {
        zero_matrix <- xts(rep(0, length(indx)), indx)
        Trend <- S <- ARMA <- Irregular <- zero_matrix
        colnames(Trend) <- "Trend"
        colnames(S) <- "Seasonal"
        colnames(ARMA) <- "ARMA"
        colnames(Irregular) <- "Irregular"
    }
    # states includes seed (t=0) state
    states <- states[-nrow(states), ]
    if (idx["Level","n"] > 0) {
        Level <- xts(states[, idx["Level","start"]:idx["Level","end"]], indx)
        if (simplify) Trend <- Trend + Level
    } else {
        Level <- NULL
    }
    if (idx["Slope","n"] > 0) {
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- xts(w_t[,nstart:nend] * states[, idx["Slope","start"]:idx["Slope","end"]], indx)
        colnames(Slope) <- "Slope"
        if (simplify) Trend <- Trend + Slope
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
                Seasonal[,i] <- t(w_t[,nstart[i]:(nstart[i] + K[i] - 1)] %*% t(states[, nstart[i]:(nstart[i] + K[i] - 1)]))
            }
        } else {
            for (i in 1:ns) {
                Seasonal[,i] <- t(w_t[,nstart[i]:(nstart[i] + frequency[i] - 1)] %*% t(states[, nstart[i]:(nstart[i] + frequency[i] - 1)]))
            }
        }
        Seasonal <- xts(Seasonal, indx)
        if (simplify) S <- S + xts(rowSums(Seasonal), indx)
    } else {
        Seasonal <- S <- NULL
    }
    if (idx["AR","n"] > 0) {
        AR <- t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[, idx["AR","start"]:idx["AR","end"]]))
        AR <- xts(AR, indx)
        colnames(AR) <- paste0("AR",idx["AR","n"])
        if (simplify) ARMA <- ARMA + AR
    } else {
        AR <- NULL
    }
    if (idx["MA","n"] > 0) {
        MA <- t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[, idx["MA","start"]:idx["MA","end"]]))
        MA <- xts(MA, indx)
        colnames(MA) <- paste0("MA",idx["MA","n"])
        if (simplify) ARMA <- ARMA + MA
    } else {
        MA <- NULL
    }
    if (idx["AR","n"] == 0 & idx["MA","n"] == 0) {
        ARMA <- NULL
    }
    if (idx["X","n"] > 0) {
        beta <- matrix(object$parmatrix[which(grepl("kappa", object$parmatrix$parameters))]$optimal, ncol = 1)
        X <- object$spec$xreg$xreg
        xreg <- xts(X %*% beta, indx)
        colnames(xreg) <- "X"
    } else {
        xreg <- NULL
    }
    Irregular <- residuals(object, raw = TRUE)
    colnames(Irregular) <- "Irregular"
    if (simplify) {
        decomposition <- cbind(Trend, S, ARMA, xreg, Irregular)
    } else {
        decomposition <- cbind(Level, Slope, Seasonal, AR, MA, xreg, Irregular)
    }
    return(decomposition)
}


tsdecompose.tsissm.predict <- function(object, simplify = FALSE, ...)
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
    if (simplify) {
        zero_matrix <- matrix(0, nrow = nsim, ncol = object$h)
        empty_list <- list(distribution = zero_matrix, original_series = xts(rep(0, NROW(object$original_series)), index(object$original_series)))
        Trend <- Seasonal <- ARMA <- Irregular <- empty_list
    }
    if (idx["Level","n"] > 0) {
        Level <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(states[, idx["Level","start"]:idx["Level","end"], i],nrow = 1)
        }))
        # start at time zero
        Level <- Level[,1:(ncol(Level) - 1)]
        colnames(Level) <- fdates
        class(Level) <- "tsmodel.distribution"
        attr(Level, "date_class") <- date_class
        Level <- list(original_series = object$decomp$Level, distribution = Level)
        class(Level) <- "tsmodel.predict"
        L[[1]] <- Level
        k <- k + 1
        if (simplify) {
            Trend <- Level
        }
    }
    if (idx["Slope","n"] > 0) {
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(w_t[,nstart:nend] * states[, idx["Slope","start"]:idx["Slope","end"], i], nrow = 1)
        }))
        Slope <- Slope[,1:(ncol(Slope) - 1)]
        colnames(Slope) <- fdates
        class(Slope) <- "tsmodel.distribution"
        attr(Slope, "date_class") <- date_class
        Slope <- list(original_series = object$decomp$Slope, distribution = Slope)
        class(Slope) <- "tsmodel.predict"
        L[[k]] <- Slope
        k <- k + 1
        if (simplify) {
            Trend$distribution <- Trend$distribution + Slope$distribution
            Trend$original_series <- Trend$original_series + Slope$original_series
        }
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
                    matrix(t(w_t[,nstart[j]:(nstart[j] + K[j] - 1)] %*% t(states[, nstart[j]:(nstart[j] + K[j] - 1), i])), nrow = 1)
                }))
                tmp <- tmp[,1:(ncol(tmp) - 1)]
                colnames(tmp) <- fdates
                class(tmp) <- "tsmodel.distribution"
                attr(tmp, "date_class") <- date_class
                tmp <- list(original_series = object$decomp[,k], distribution = tmp)
                class(tmp) <- "tsmodel.predict"
                L[[k]] <- tmp
                if (simplify) {
                    Seasonal$distribution <- Seasonal$distribution + tmp$distribution
                    Seasonal$original_series <- Seasonal$original_series + tmp$original_series
                }
                k <- k + 1
            }
        } else {
            for (i in 1:ns) {
                tmp <- do.call(rbind, lapply(1:nsim, function(i){
                    matrix(t(w_t[,nstart[j]:(nstart[j] + frequency[j] - 1)] %*% t(states[, nstart[j]:(nstart[j] + frequency[j] - 1), i])), nrow = 1)
                }))
                colnames(tmp) <- fdates
                tmp <- tmp[,1:(ncol(tmp) - 1)]
                class(tmp) <- "tsmodel.distribution"
                attr(tmp, "date_class") <- date_class
                tmp <- list(original_series = object$decomp[,k], distribution = tmp)
                class(tmp) <- "tsmodel.predict"
                L[[k]] <- tmp
                if (simplify) {
                    Seasonal$distribution <- Seasonal$distribution + tmp$distribution
                    Seasonal$original_series <- Seasonal$original_series + tmp$original_series
                }
                k <- k + 1
            }
        }
    } else {
        Seasonal <- NULL
    }
    if (idx["AR","n"] > 0) {
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[, idx["AR","start"]:idx["AR","end"], i])), nrow = 1)
        }))
        tmp <- tmp[,1:(ncol(tmp) - 1)]
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        tmp <- list(original_series = object$decomp[,k], distribution = tmp)
        class(tmp) <- "tsmodel.predict"
        L[[k]] <- tmp
        k <- k + 1
        if (simplify) {
            ARMA$distribution <- ARMA$distribution + tmp$distribution
            ARMA$original_series <- ARMA$original_series + tmp$original_series
        }
    }
    if (idx["MA","n"] > 0) {
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[, idx["MA","start"]:idx["MA","end"], i])), nrow = 1)
        }))
        tmp <- tmp[,1:(ncol(tmp) - 1)]
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        tmp <- list(original_series = object$decomp[,k], distribution = tmp)
        class(tmp) <- "tsmodel.predict"
        L[[k]] <- tmp
        k <- k + 1
        if (simplify) {
            ARMA$distribution <- ARMA$distribution + tmp$distribution
            ARMA$original_series <- ARMA$original_series + tmp$original_series
        }
    }
    if (idx["AR","n"] == 0 & idx["MA","n"] == 0) {
        ARMA <- NULL
    }
    # Innovations
    Irregular <- object$innov
    L[[k]] <- Irregular
    names(L) <- c(colnames(object$decomp[,1:(length(L) - 1)]), "Irregular")
    if (simplify) {
        if (!is.null(Seasonal)) class(Seasonal) <- "tsmodel.predict"
        if (!is.null(ARMA)) class(ARMA) <- "tsmodel.predict"
        class(Irregular) <- "tsmodel.predict"
        L <- list()
        L$Trend <- Trend
        if (!is.null(Seasonal)) L$Seasonal <- Seasonal
        if (!is.null(ARMA)) L$ARMA <- ARMA
        L$Irregular <- Irregular
    }
    return(L)
}

logLik.tsissm.estimate <- function(object, ...)
{
    parameters <- NULL
    estimate <- NULL
    np <- NROW(object$parmatrix[estimate == 1]) + ncol(object$model$xseed)
    r <- na.omit(residuals(object, raw = TRUE))
    v <- sum(r^2)
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    n <- NROW(object$spec$target$y_orig[which(object$spec$good == 1)])
    if (object$spec$transform$name == "box-cox") {
        loglik <- -0.5 * n * log(2 * pi * (v/n)) - (1/(2 * (v/n))) * v + (lambda - 1) * sum(log(object$spec$target$y_orig[which(object$spec$good == 1)]))
    } else {
        loglik <- -0.5 * n * log(2 * pi * (v/n)) - (1/(2 * (v/n))) * v
    }
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
    # yt <- object$spec$transform$transform(object$spec$target$y_orig, lambda = lambda)
    # ft <- object$spec$transform$transform(as.numeric(fitted(object)), lambda = lambda)
    # r <- yt - ft
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
    m_mase <- mase(actual, colMeans(object$distribution), object$original_series, frequency = object$frequency[1])
    if (!is.null(alpha)) {
        m_mis <- mis(actual, lower = apply(object$distribution, 2, quantile, alpha/2), upper = apply(object$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
    } else {
        m_mis <- as.numeric(NA)
    }
    m_crps <- crps(actual, object$distribution)
    return(data.frame("h" = n, "MAPE" = m_mape, "MASE" = m_mase, "MSLRE" = m_mslre, "BIAS" = m_bias, "MIS" = m_mis,"CRPS" = m_crps))
}

