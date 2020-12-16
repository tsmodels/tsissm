predict.tsissm.estimate <- function(object, h = 12, newxreg = NULL, nsim = 1000, forc_dates = NULL, innov = NULL, init_states = NULL, ...)
{
    parameters <- NULL
    if (!is.null(forc_dates)) {
        if (h != length(forc_dates))
            stop("\nforc_dates must have length equal to h")
    }
    if (!object$spec$xreg$include_xreg) {
        newxreg = matrix(0, ncol = 1, nrow = h)
        if (is.null(forc_dates)) {
            forc_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
        }
    } else {
        if (!is.null(newxreg)) {
            forc_dates = index(newxreg)
        }
        else {
            if (is.null(forc_dates)) {
                forc_dates = future_dates(tail(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
            newxreg <- xts(matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h), forc_dates)
            colnames(newxreg) <- colnames(object$spec$xreg$xreg)
        }
    }
    act <- object$spec$transform$transform(object$spec$target$y_orig, lambda = object$parmatrix[parameters == "lambda"]$optimal)
    fit <- object$spec$transform$transform(object$model$fitted, lambda = object$parmatrix[parameters == "lambda"]$optimal)
    res <- act - fit
    sigma.res <- sd(res)
    
    if (is.null(innov)) {
        E <- matrix(rnorm(h * nsim, 0, sigma.res), ncol = h, nrow = nsim)
    } else {
        if (length(innov) != (h * nsim)) {
            stop("\nlength innov must be nsim x h")
        }
        if (any(innov == 0)) {
            innov[which(innov == 0)] <- 1e-12
        }
        if ( any(innov == 1)) {
            innov[which(innov == 1)] <- 1 - 1e-12
        }
        innov <- matrix(innov, h, nsim)
        E <- qnorm(innov, sd = sigma.res)
    }
    xseed <- tail(object$model$states, 1)
    if (!is.null(init_states)) {
        if (length(as.vector(init_states) != ncol(xseed))) {
            stop(paste0("\ninit_states must be a vector of length ", ncol(xseed)))
        } else {
            xseed <- matrix(as.numeric(init_states), nrow = 1, ncol = ncol(xseed))
        }
    }
    pars <- object$parmatrix$optimal
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    mdim = c(object$spec$dims[1], nsim, h, ncol(object$spec$xreg$xreg))
    f <- isspredict(f0_ = S[list("F0")]$values,
                    f1_ = S[list("F1")]$values,
                    f2_ = S[list("F2")]$values,
                    w_ = as.numeric(S[list("w")]$values),
                    g_ = as.numeric(S[list("g")]$values),
                    error_ = as.vector(E),
                    xseed_ = xseed,
                    xreg_ = as.numeric(newxreg),
                    kappa_ = S[list("kappa")]$values,
                    mdim = mdim)
    date_class <- attr(object$spec$target$sampling, "date_class")
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    ysim <- object$spec$transform$inverse(f$simulated[,-1], lambda = lambda)
    if (NCOL(ysim) == 1) ysim <- matrix(ysim, ncol = 1)
    colnames(ysim) <- as.character(forc_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    states <- f$states
    spec <- object$spec
    spec$parmatrix <- object$parmatrix
    zList <- list(original_series = xts(object$spec$target$y_orig, object$spec$target$index),
                  distribution = ysim, mean = zoo(colMeans(ysim), forc_dates), h = h,
                  spec = spec, states = states, decomp = tsdecompose(object))
    class(zList) <- c("tsissm.predict","tsmodel.predict")
    return(zList)
}
