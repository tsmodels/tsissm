#' Model Specification
#'
#' @description Specifies an ISSM model prior to estimation.
#' @details The specification object holds the information and data which is
#' then passed to the maximum likelihood estimation routines.
#' @param y an xts vector.
#' @param slope (Logical) slope component.
#' @param slope_damped (Logical) slope dampening component.
#' @param seasonal (Logical) seasonal component(s).
#' @param seasonal_frequency vector of numeric seasonal frequencies. For 
#' trigonometric this can be fractional, but must be integer for regular 
#' seasonality.
#' @param seasonal_type either trigonometric or regular. The latter currently 
#' does not allow multiple seasonality.
#' @param seasonal_harmonics the number of harmonics per seasonal frequency for 
#' the trigonometric type.
#' @param ar AR order.
#' @param ma MA order.
#' @param xreg an xts matrix of external regressors.
#' @param transformation a valid transformation for y from the \dQuote{tstransform} 
#' function in the \dQuote{tsaux} package (currently box-cox or logit are available).
#' @param lambda the Box Cox lambda. If not NULL, then either a numeric value or NA 
#' denoting automatic calculation.
#' @param lower lower bound for the transformation.
#' @param upper upper bound for the transformation.
#' @param sampling (optional) sampling frequency of the dataset. If NULL, will 
#' try to identify from the timestamps of y. This is useful for plotting and 
#' extending the timestamps in the prediction horizon.
#' @param ... not used.
#' @details The specification performs some sanity checks on the arguments provided 
#' and sets up the required state space matrices and parameters which are used in 
#' the estimation stage.
#' @return An object of class \dQuote{tsissm.spec} with the following slots:\cr
#' \item{target}{A list with original data series, the data series index and the 
#' sampling frequency}
#' \item{slope}{A list with details about the slope state}
#' \item{seasonal}{A list with details about the seasonal state}
#' \item{xreg}{A list with details on the external regressors}
#' \item{transform}{A list with details on the transformation}
#' \item{arma}{A list with details on the ARMA state}
#' \item{S}{A data.table with the vectorized state matrices}
#' \item{dims}{A vector with dimensions and flags used in the estimation code}
#' \item{parmatrix}{A data.table of the model parameters}
#' \item{idmatrix}{A matrix with index information on the parameters}
#' @references De Livera, Alysha M and Hyndman, Rob J and Snyder, Ralph D, 2011, 
#' Forecasting time series with complex seasonal patterns using exponential smoothing, 
#' \emph{Journal of the American Statistical Association}, \bold{106(496)}, 1513--1527.
#' @aliases issm_modelspec
#' @rdname issm_modelspec
#' @export
#'
#'
#'
#'
issm_modelspec <- function(y, slope = TRUE, slope_damped = FALSE, seasonal = FALSE, 
                           seasonal_frequency = 1, seasonal_type = c("trigonometric","regular"), 
                           seasonal_harmonics = NULL, ar = 0, ma = 0, xreg = NULL, 
                           transformation = "box-cox", lambda = 1, lower = 0, upper = 1, 
                           sampling = NULL, ...)
{
    if (!is.xts(y)) {
        stop("\ny is not an xts object...refusing to continue. Please fix and resubmit.")
    }
    if (any(y <= 0, na.rm = TRUE)) 
        stop("\nThe issm model only supports strictly positive data.")
    n <- NROW(y)
    if (is.null(sampling)) {
        sampling <- sampling_frequency(index(y))
    }
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
    }
    if (!is.null(xreg)) {
        if (nrow(xreg) != n) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (n/2)) warning("\nnumber of regressors greater than half the length of y!")
        if (is.null(colnames(xreg))) 
            colnames(xreg) <- paste0("reg.", 1:ncol(xreg))
        include_xreg <- TRUE
        xreg <- coredata(xreg)
    } else {
        xreg <- matrix(0, ncol = 1, nrow = n)
        colnames(xreg) <- paste0("reg.", 1)
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
        seasonal_type <- match.arg(seasonal_type[1], c("trigonometric","regular"))
        if (seasonal_type == "regular" & length(seasonal_frequency) > 1) {
            stop("\nregular seasonality currently only supports single frequency models.")
        }
        if (any(seasonal_frequency > (n/2))) {
            warning("\nseasonal_frequency contains a value(s) which are greater than half the length of y!")
        }
        if (length(seasonal_frequency) > 2) {
            seasonal_frequency <- sort(seasonal_frequency)
            if (any(diff(seasonal_frequency) == 0)) {
                stop("\nfound duplicate values in seasonal_frequency.")
            }
        }
    } else {
        seasonal_type <- "regular"
        seasonal_harmonics <- NULL
        seasonal_frequency <- 1
    }
    spec$seasonal$seasonal_harmonics <- seasonal_harmonics[1:length(seasonal_frequency)]
    spec$seasonal$seasonal_frequency <- seasonal_frequency
    spec$seasonal$seasonal_type <- match.arg(seasonal_type[1], c("trigonometric", "regular"))
    spec$xreg$include_xreg <- include_xreg
    spec$xreg$xreg <- xreg
    # 
    transformation <- match.arg(transformation[1], c("box-cox","logit"))
    if (transformation == "logit") {
        lambda <- 1
        transform <- tstransform(method = transformation, lower = lower, upper = upper)
        transform$lambda <- 1
        transform$name <- "logit"
        include_lambda <- FALSE
        transform$include_lambda <- FALSE
        transform$lower <- lower
        transform$upper <- upper
        y <- transform$transform(y)
        spec$target$y <- as.numeric(y)
    } else {
        if (is.null(lambda)) lambda <- 1
        if (is.na(lambda)) {
            include_lambda <- TRUE
            transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
            tmp <- transform$transform(y = y, frequency = seasonal_frequency[1])
            transform$lambda <- attr(tmp, "lambda")
            transform$include_lambda <- TRUE
            lambda <- transform$lambda
            transform$lower <- lower
            transform$upper <- upper
        } else {
            include_lambda <- FALSE
            transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
            transform$include_lambda <- FALSE
            transform$lambda <- lambda
            transform$lower <- lower
            transform$upper <- upper
        }
        transform$name <- "box-cox"
        spec$target$y <- as.numeric(y)
    }
    spec$transform <- transform
    spec$arma$order <- c(ar, ma)
    M <- ss_matrices(y = y, slope = spec$slope$include_slope, 
                     damped_slope = spec$slope$include_damped, frequency = spec$seasonal$seasonal_frequency, 
                     type = spec$seasonal$seasonal_type, K = spec$seasonal$seasonal_harmonics, 
                     ar = spec$arma$order[1], ma = spec$arma$order[2], include_xreg = include_xreg, 
                     xreg = xreg)
    S <- M$setup
    ipars <- issm_init_pars(lambda = lambda, frequency = seasonal_frequency)
    S <- rbind(S, data.table(matrix = "xreg", values = as.vector(coredata(xreg)), pars = NA))
    pars <- M$pars
    est <- rep(1, length(pars))
    dims <- M$dims
    # tells the C++ estimation codes to not apply a box cox transform
    if (transformation == "logit") {
        dims <- c(dims, 1)
    } else {
        dims <- c(dims, 0)
    }
    lower <- M$lower
    upper <- M$upper
    if (slope) {
        if (!spec$slope$include_damped) {
            ix <- which(grepl("phi", pars))
            est[ix] <- 0
            lower[ix] <- 1
            upper[ix] <- 1
        }
    }
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
        lower_x <- -1e8
        upper_x <-  1e8
        S <- rbind(S, data.table(matrix = "kappa", values = rep(as.numeric(NA), ncol(xreg)), pars = paste0("kappa.",1:ncol(xreg))))
        pars <- c(pars, paste0("kappa.",1:ncol(xreg)))
        lower <- c(lower, rep(lower_x, ncol(xreg)))
        upper <- c(upper, rep(upper_x, ncol(xreg)))
        est <- c(est, rep(1, ncol(xreg)))
    } else {
        S <- rbind(S, data.table(matrix = "kappa", values = rep(as.numeric(0), ncol(xreg)), pars = paste0("kappa.", 1:ncol(xreg))))
        pars <- c(pars, paste0("kappa.", 1:ncol(xreg)))
        lower <- c(lower, rep(0, ncol(xreg)))
        upper <- c(upper, rep(0, ncol(xreg)))
        est <- c(est, rep(0, ncol(xreg)))
    }
    init <- ipars[sapply(1:length(pars), function(i) which(grepl(gsub("[^a-zA-A]", "", pars[i]), ipars$variable)))]
    iscale <- init$scale
    init <- init$value
    names(init) <- pars
    parmatrix <- data.table(parameters = pars, initial = init, lower = lower, upper = upper, estimate = est, scale = iscale)
    setkeyv(S, "matrix")
    spec$S <- S
    spec$dims <- dims
    spec$parmatrix <- parmatrix
    spec$idmatrix <- M$idmatrix
    spec$good <- good
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

#' Model Specification Extractor
#'
#' @description Extracts a model specification (class \dQuote{tsissm.spec}) from
#' an object of class \dQuote{tsissm.estimate}.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param y an optional new xts vector.
#' @param lambda an optional lambda parameter for the Box Cox transformation (if
#' previously used).
#' @param xreg an optional matrix of regressors.
#' @param ... not currently used.
#' @note This function is used by other functions in the package such as the
#' backtest which requires rebuilding the specification for each re-estimation
#' step with updated data but keeping all else equal.
#' @return An object of class \dQuote{tsissm.spec}.
#' @aliases tsspec
#' @method tsspec tsissm.estimate
#' @rdname tsspec
#' @export
#'
#'
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
                   transformation = object$spec$transform$name,
                   lower = object$spec$transform$lower,
                   upper = object$spec$transform$upper,
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
