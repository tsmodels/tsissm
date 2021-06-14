issm_modelspec <- function(y, slope = TRUE, slope_damped = FALSE, seasonal = FALSE, seasonal_frequency = 1,
                           seasonal_type = c("trigonometric","regular"), seasonal_harmonics = NULL,
                           ar = 0, ma = 0, xreg = NULL, lambda = 1, sampling = NULL, ...)
{
    # checks
    if (!is.xts(y)) {
        stop("\ny is not an xts object...refusing to continue. Please fix and resubmit.")
    }
    n <- NROW(y)
    good <- rep(1, n)
    if (any(is.na(y))) {
        ix <- which(is.na(y))
        good[ix] <- 0
    }
    if (is.null(sampling)) {
        sampling <- sampling_frequency(index(y))
    }
    if (!is.null(xreg)) {
        if (nrow(xreg) != n) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (n/2)) warning("\nnumber of regressors greater than half the length of y!")
        if (is.null(colnames(xreg))) colnames(xreg) <- paste0("reg.",1:ncol(xreg))
        include_xreg <- TRUE
        xreg <- coredata(xreg)
    } else {
        # fix to zero valued matrix
        xreg <- matrix(0, ncol = 1, nrow = n)
        colnames(xreg) <- paste0("reg.",1)
        include_xreg <- FALSE
    }
    spec <- list()
    spec$target$y_orig <- as.numeric(y)
    spec$target$index <- index(y)
    spec$target$sampling <- sampling
    spec$slope$include_slope <- as.logical(slope)
    if (slope) {
        if (slope_damped) {
            include_damped <- TRUE
        } else {
            include_damped <- FALSE
        }
        spec$slope$include_damped <- include_damped
    } else {
        spec$slope$include_damped <- FALSE
    }
    spec$seasonal$include_seasonal <- as.logical(seasonal)
    if (seasonal) {
        if (any(is.null(seasonal_frequency))) {
            stop("\nseasonal_frequency cannot be NULL if seasonal is set to TRUE")
        }
        if (any(seasonal_frequency < 2)) {
            stop("\nseasonal_frequency must be 2 or greater")
        }
        if (any((seasonal_frequency) %% 1 != 0)) {
            if (seasonal_type[1] == "regular") {
                stop("\nseasonal_type must be trigonometric when seasonal.frequency is not a whole number.")
            }
        }
        if (seasonal_type == "regular" & length(seasonal_frequency) > 1) {
            stop("\nregular seasonality currently only supports single frequency models.")
        }
        if (any(seasonal_frequency > (n/2))) {
            warning("\nseasonal_frequency contains a value(s) which are greater than half the length of y!")
        }
        # length of frequency vs harmonics for trig
        if (length(seasonal_frequency) > 2) {
            seasonal_frequency <- sort(seasonal_frequency)
            if (any(diff(seasonal_frequency) == 0)) {
                stop("\nfound duplicate values in seasonal_frequency.")
            }
        }
    } else {
        seasonal_harmonics <- NULL
        seasonal_frequency <- 1
    }
    spec$seasonal$seasonal_harmonics <- seasonal_harmonics[1:length(seasonal_frequency)]
    spec$seasonal$seasonal_frequency <- seasonal_frequency
    spec$seasonal$seasonal_type <- match.arg(seasonal_type[1], c("trigonometric","regular"))
    spec$xreg$include_xreg <- include_xreg
    spec$xreg$xreg <- xreg
    
    if (is.null(lambda)) lambda <- 1
    if (is.na(lambda)) {
        include_lambda <- TRUE
        transform <- box_cox(lambda = lambda, lower = 0, upper = 1.5)
        tmp <- transform$transform(y = y, frequency = seasonal_frequency[1])
        transform$lambda <- attr(tmp, "lambda")
        transform$include_lambda <- TRUE
        lambda <- transform$lambda
    } else {
        include_lambda <- FALSE
        transform <- box_cox(lambda = lambda, lower = 0, upper = 1.5)
        transform$include_lambda <- FALSE
        transform$lambda <- lambda
        tmp <- transform$transform(y = y, frequency = seasonal_frequency[1])
    }
    spec$transform <- transform
    spec$arma$order <- c(ar, ma)
    M <- ss_matrices(y = y, slope = spec$slope$include_slope, damped_slope = spec$slope$include_damped,
                     frequency = spec$seasonal$seasonal_frequency, type = spec$seasonal$seasonal_type,
                     K = spec$seasonal$seasonal_harmonics, ar = spec$arma$order[1], ma = spec$arma$order[2],
                     include_xreg = include_xreg, xreg = xreg)
    S <- M$setup
    ipars <- issm_init_pars(lambda = lambda, frequency = seasonal_frequency)
    S <- rbind(S, data.table(matrix = "xreg", values = as.vector(coredata(xreg)), pars = NA))
    pars <- M$pars
    est <- rep(1, length(pars))
    dims <- M$dims
    lower <- M$lower
    upper <- M$upper
    if (slope) {
        if (!spec$slope$include_damped) {
            ix <- which(grepl("phi",pars))
            est[ix] <- 0
            lower[ix] <- 1
            upper[ix] <- 1
        }
    }
    # else already included
    if (include_lambda) {
        pars <- c(pars, "lambda")
        lower <- c(lower, 1e-12)
        upper <- c(upper, 1.5)
        est <- c(est, 1)
    } else {
        pars <- c(pars, "lambda")
        lower <- c(lower, lambda)
        upper <- c(upper, lambda)
        est <- c(est, 0)
    }
    if (include_xreg) {
        S <- rbind(S, data.table(matrix = "kappa", values = rep(as.numeric(NA), ncol(xreg)), pars = paste0("kappa.",1:ncol(xreg))))
        pars <- c(pars, paste0("kappa.",1:ncol(xreg)))
        lower <- c(lower, rep(-Inf, ncol(xreg)))
        upper <- c(upper, rep( Inf, ncol(xreg)))
        est <- c(est, rep(1, ncol(xreg)))
    } else {
        S <- rbind(S, data.table(matrix = "kappa", values = rep(as.numeric(0), ncol(xreg)), pars = paste0("kappa.",1:ncol(xreg))))
        pars <- c(pars, paste0("kappa.",1:ncol(xreg)))
        lower <- c(lower, rep(0, ncol(xreg)))
        upper <- c(upper, rep(0, ncol(xreg)))
        est <- c(est, rep(0, ncol(xreg)))
    }
    init <- ipars[sapply(1:length(pars), function(i) which(grepl(gsub("[^a-zA-A]","",pars[i]), ipars$variable)))]
    iscale <- init$scale
    init <- init$value
    names(init) <- pars
    parmatrix <- data.table(parameters = pars, initial = init, lower = lower, upper = upper, estimate = est, scale = iscale)
    setkeyv(S, "matrix")
    spec$S <- S
    spec$dims <- dims
    spec$parmatrix <- parmatrix
    spec$idmatrix <- M$idmatrix
    class(spec) <- c("tsissm.spec", "tsmodel.spec")
    return(spec)
}

issm_init_pars <- function(lambda, frequency)
{
    rbind(
        data.table(variable = "alpha", value = 0.09, scale = 0.01),
        data.table(variable = "beta", value = 0.01, scale = 0.01),
        data.table(variable = "gamma", value = 0.0, scale = 1e-05),
        data.table(variable = "phi", value = 0.98, scale = 0.01),
        data.table(variable = "kappa", value = 0.0, scale = 1),
        data.table(variable = "theta", value = 0.01, scale = 0.1),
        data.table(variable = "psi", value = 0.01, scale = 0.1),
        data.table(variable = "lambda", value = lambda, scale = 0.001))
}

tsspec.tsissm.estimate <- function(object, y = NULL, lambda = NULL, xreg = NULL, ...)
{
    if (is.null(y)) {
        y <- xts(object$spec$target$y_orig, object$spec$target$index)
    }
    if (!is.null(xreg)) {
        xreg <- coredata(xreg)
        if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
    } else {
        if (object$spec$xreg$include_xreg) {
            xreg <- object$spec$xreg$xreg
            if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
            if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
        } else {
            xreg <- NULL
        }
    }
    if (is.null(lambda)) {
        lambda <- object$spec$transform$lambda
    }
    issm_modelspec(y, slope = object$spec$slope$include_slope,
                   slope_damped = object$spec$slope$include_damped,
                   seasonal = object$spec$seasonal$include_seasonal,
                   seasonal_frequency = object$spec$seasonal$seasonal_frequency,
                   seasonal_type = object$spec$seasonal$seasonal_type,
                   seasonal_harmonics = object$spec$seasonal$seasonal_harmonics,
                   ar = object$spec$arma$order[1],
                   ma = object$spec$arma$order[2],
                   xreg = xreg, lambda = lambda,
                   sampling = object$spec$target$sampling)
}

get_pars.tsissm.spec <- function(object)
{
    return(object$parmatrix)
}

set_pars.tsissm.spec <- function(object, value)
{

}
