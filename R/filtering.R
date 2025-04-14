#' Online Model Filtering
#'
#' @description Online filter which updates the states and fitted values using
#' new data.
#' @param object an object of class \dQuote{tsissm.estimate} or \dQuote{tsissm.selection}.
#' @param y an xts vector of new information related to y. The function checks
#' whether y contains indices (dates) which are not in the passed object and
#' only filters new information.
#' @param newxreg An xts matrix of new information related to external regressors
#' (if those were used in the original model estimated).
#' @param ... not currently used.
#' @details The new data is filtered (1 step ahead) using the last state of the
#' object. Once this is complete, the object is updated with the new states
#' and information so that the process can be continued on the same object as
#' new information arrives.
#' @return An object of class \dQuote{tsissm.estimate} or  \dQuote{tsissm.selection}.
#' @aliases tsfilter
#' @method tsfilter tsissm.estimate
#' @rdname tsfilter
#' @export
#'
#'
tsfilter.tsissm.estimate <- function(object, y = NULL, newxreg = NULL, ...)
{
    if (!is.null(y)) {
        if (!is.xts(y)) stop("\ny must be an xts object")
    }
    out <- switch(object$spec$variance$type,
                  "constant" = .filter_constant(object, y, newxreg),
                  "dynamic" = .filter_dynamic(object, y, newxreg))
    return(out)
}


#' @method tsfilter tsissm.selection
#' @rdname tsfilter
#' @export
#'
#'
tsfilter.tsissm.selection <- function(object, y = NULL, newxreg = NULL, ...)
{
    if (!is.null(y)) {
        if (!is.xts(y)) stop("\ny must be an xts object")
    }
    out <- future_lapply(1:length(object$models), function(i){
        tsfilter(object$models[[i]], y = y, newxreg = newxreg)
    }, future.seed = TRUE, future.packages = "tsissm")
    out <- eval(out)
    object$models <- out
    return(object)
}

.filter_constant <- function(object, y = NULL, newxreg = NULL)
{
    parameters <- NULL
    yold <- xts(object$spec$target$y_orig, object$spec$target$index)
    ynew <- y
    exc <- which(index(ynew) %in% index(yold))
    if (length(exc) == 0) {
        y <- ynew
    } else{
        y <- ynew[-exc]
        if (NROW(y) == 0) {
            warning("\nno new data in y which is not already in object!")
            return(object)
        }
    }
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
    }
    
    if (object$spec$xreg$include_xreg) {
        nx <- NCOL(object$spec$xreg$xreg)
        if (!is.null(newxreg)) {
            if (ncol(newxreg) != nx) stop(paste0("\nExpected ", nx, " columns in newxreg but got ", NCOL(newxreg)))
            newcnames <- colnames(newxreg)
            oldcnames <- colnames(object$spec$xreg$xreg)
            if (!is.null(newcnames) & !is.null(oldcnames)) {
                if (!all(sort(oldcnames) %in% sort(newcnames))) {
                    stop("\ncolnames of newxreg do not match those of original xreg")
                }
                newxreg <- newxreg[, oldcnames]
            }
            X <- coredata(newxreg)
            if (length(exc) > 0) {
                X <- X[-exc,,drop = FALSE]
            }
        } else {
            X <- matrix(0, ncol = nx, nrow = nrow(y))
        }
        X <- rbind(matrix(0, ncol = nx, nrow = 1), X)
    } else {
        X <- matrix(0, ncol = 1, nrow = NROW(y) + 1)
    }
    newindex <- index(y)
    yneworig <- y
    xseed <- tail(object$model$states, 1)
    pars <- object$parmatrix$value
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    parnames <- names(pars)
    mdim <- object$spec$dims
    mdim[2] <- NROW(y)
    n_states <- mdim[1]
    timesteps <- mdim[2]
    f0 <- matrix(S[list("F0")]$values, n_states, n_states)
    f1 <- matrix( S[list("F1")]$values, n_states, n_states)
    f2 <- matrix( S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    y <- as.numeric(y)
    kappa <- as.numeric(S[list("kappa")]$values)
    xseed <- as.vector(xseed)
    f <- .issfilter_constant(F0 = f0, F1 = f1, F2 = f2, w = w, g = g, y = c(1.0, y), X = X,  kappa = kappa, 
                             xseed = xseed, good = c(0, good), mdim = mdim, lambda = as.numeric(pars["lambda"]))
    
    object$spec$target$y_orig <- c(object$spec$target$y_orig, as.numeric(yneworig))
    object$spec$target$index <- c(object$spec$target$index, newindex)
    object$spec$good <- c(object$spec$good, good)
    object$spec$mdim[2] <- length(object$spec$target$y_orig)
    if (!is.null(newxreg)) {
        object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, X[-1, , drop = FALSE])
    } else {
        object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = NROW(yneworig)))
    }
    filt <- object$spec$transform$inverse(f$fitted, lambda = as.numeric(pars["lambda"]))
    object$spec$target$y <- c(object$spec$target$y, object$spec$transform$transform(as.numeric(y), lambda = as.numeric(pars["lambda"])))
    filt <- filt[-1]
    err <- as.numeric(y) - filt
    object$model$fitted <- c(object$model$fitted, filt)
    object$model$states <- rbind(object$model$states, f$states[-1,,drop = FALSE])
    object$model$residuals <- c(object$model$residuals, err)
    object$model$error <- c(object$model$error, f$error[-1])
    
    good_index <- which(object$spec$good == 1)
    sig <- object$model$sigma
    n <- length(object$spec$target$y_orig)
    ngood <- sum(object$spec$good)
    if (object$spec$distribution$type == "norm") {
        nll <- -log(dnorm(object$model$error[good_index]/sig, 0, 1)/sig) - (pars["lambda"] - 1.0) * log(as.numeric(object$spec$target$y_orig[good_index]))
    } else if (object$spec$distribution$type == "std") {
        shape <- object$parmatrix[parameters == "shape"]$value
        nll <- -log(ddist("std", object$model$error[good_index]/sig, 0, 1, shape = shape)/sig) - (pars["lambda"] - 1.0) * log(as.numeric(object$spec$target$y_orig[good_index]))
    } else if (object$spec$distribution$type == "jsu") {
        skew <- object$parmatrix[parameters == "skew"]$value
        shape <- object$parmatrix[parameters == "shape"]$value
        nll <- -log(ddist("jsu", object$model$error[good_index]/sig, 0, 1, skew = skew, shape = shape)/sig) - (pars["lambda"] - 1.0) * log(as.numeric(object$spec$target$y_orig[good_index]))
    }
    loglik <- ngood/n * sum(nll)
    object$model$nll <- nll
    object$model$loglik <- loglik
    
    
    return(object)
}

.filter_dynamic <- function(object, y = NULL, newxreg = NULL)
{
    group <- parameters <- NULL
    yold <- xts(object$spec$target$y_orig, object$spec$target$index)
    ynew <- y
    exc <- which(index(ynew) %in% index(yold))
    if (length(exc) == 0) {
        y <- ynew
    } else{
        y <- ynew[-exc]
        if (NROW(y) == 0) {
            warning("\nno new data in y which is not already in object!")
            return(object)
        }
    }
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
    }
    
    if (object$spec$xreg$include_xreg) {
        nx <- NCOL(object$spec$xreg$xreg)
        if (!is.null(newxreg)) {
            if (ncol(newxreg) != nx) stop(paste0("\nExpected ", nx, " columns in newxreg but got ", NCOL(newxreg)))
            newcnames <- colnames(newxreg)
            oldcnames <- colnames(object$spec$xreg$xreg)
            if (!is.null(newcnames) & !is.null(oldcnames)) {
                if (!all(sort(oldcnames) %in% sort(newcnames))) {
                    stop("\ncolnames of newxreg do not match those of original xreg")
                }
                newxreg <- newxreg[, oldcnames]
            }
            X <- coredata(newxreg)
            if (length(exc) > 0) {
                X <- X[-exc,,drop = FALSE]
            }
        } else {
            X <- matrix(0, ncol = nx, nrow = nrow(y))
        }
        X <- rbind(matrix(0, ncol = nx, nrow = 1), X)
    } else {
        X <- matrix(0, ncol = 1, nrow = NROW(y) + 1)
    }
    newindex <- index(y)
    yneworig <- y
    xseed <- tail(object$model$states, 1)
    pars <- object$parmatrix$value
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    parnames <- names(pars)
    mdim <- object$spec$dims
    mdim[2] <- NROW(y)
    n_states <- mdim[1]
    timesteps <- mdim[2]
    f0 <- matrix(S[list("F0")]$values, n_states, n_states)
    f1 <- matrix( S[list("F1")]$values, n_states, n_states)
    f2 <- matrix( S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    y <- as.numeric(y)
    kappa <- as.numeric(S[list("kappa")]$values)
    eta <- object$parmatrix[group == "eta"]$value
    delta <- object$parmatrix[group == "delta"]$value
    omega <- object$model$target_omega
    vmodel <- object$spec$variance$vmodel
    maxpq <- vmodel[1]
    initial_arch <- tail(object$model$error^2, maxpq)
    initial_variance <- tail(object$model$sigma^2, maxpq)
    xseed <- as.vector(xseed)
    f <- .issfilter_dynamic(F0 = f0, F1 = f1, F2 = f2, w = w, g = g, y = c(1.0, y), X = X,  kappa = kappa, 
                            xseed = xseed, good = c(0, good), eta = eta, delta = delta, 
                            initial_arch = initial_arch, initial_variance = initial_variance, 
                            mdim = mdim, omega = omega, vmodel = as.integer(vmodel), 
                            lambda = as.numeric(pars["lambda"]))
    
    object$spec$target$y_orig <- c(object$spec$target$y_orig, as.numeric(yneworig))
    object$spec$target$index <- c(object$spec$target$index, newindex)
    object$spec$good <- c(object$spec$good, good)
    object$spec$mdim[2] <- length(object$spec$target$y_orig)
    if (!is.null(newxreg)) {
        object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, X[-1, , drop = FALSE])
    } else {
        object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = NROW(yneworig)))
    }
    filt <- object$spec$transform$inverse(f$fitted, lambda = as.numeric(pars["lambda"]))
    object$spec$target$y <- c(object$spec$target$y, object$spec$transform$transform(as.numeric(y), lambda = as.numeric(pars["lambda"])))
    filt <- filt[-1]
    err <- as.numeric(y) - filt
    object$model$fitted <- c(object$model$fitted, filt)
    object$model$sigma <- c(object$model$sigma, f$sigma[-seq(1,vmodel[1])])
    object$model$states <- rbind(object$model$states, f$states[-1,,drop = FALSE])
    object$model$residuals <- c(object$model$residuals, err)
    object$model$error <- c(object$model$error, f$error[-1])
    
    good_index <- which(object$spec$good == 1)
    sig <- object$model$sigma[good_index]
    n <- length(object$spec$target$y_orig)
    ngood <- sum(object$spec$good)
    z_std <- object$model$error[good_index]/sig
    if (object$spec$distribution$type == "norm") {
        nll <- -log(dnorm(z_std, 0, 1)/sig) - (pars["lambda"] - 1.0) * log(as.numeric(object$spec$target$y_orig[good_index]))
    } else if (object$spec$distribution$type == "std") {
        shape <- object$parmatrix[parameters == "shape"]$initial
        nll <- -log(ddist("std", z_std, 0, 1, shape = shape)/sig) - (pars["lambda"] - 1.0) * log(as.numeric(object$spec$target$y_orig[good_index]))
    } else if (object$spec$distribution$type == "jsu") {
        skew <- object$parmatrix[parameters == "skew"]$value
        shape <- object$parmatrix[parameters == "shape"]$value
        nll <- -log(ddist("jsu", z_std, 0, 1, skew = skew, shape = shape)/sig) - (pars["lambda"] - 1.0) * log(as.numeric(object$spec$target$y_orig[good_index]))
    }
    loglik <- ngood/n * sum(nll)
    object$model$nll <- nll
    object$model$loglik <- loglik
    
    return(object)
}
