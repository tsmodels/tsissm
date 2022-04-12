predict.tsissm.estimate <- function(object, h = 12, newxreg = NULL, nsim = 1000, forc_dates = NULL, innov = NULL, init_states = NULL, 
                                    exact_moments = TRUE, innov_type = "q", sigma_scale = NULL, ...)
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
    if (!is.null(sigma_scale)) {
        sigma_scale <- as.numeric(sigma_scale)
        if (any(sigma_scale <= 0)) stop("\nsigma_scale must be strictly positive")
        if (length(sigma_scale) == 1) sigma_scale <- rep(sigma_scale, h)
        if (length(sigma_scale) != h) stop("\nsigma_scale must be of length h or 1 (recycled to h)")
    }
    act <- object$spec$transform$transform(object$spec$target$y_orig, lambda = object$parmatrix[parameters == "lambda"]$optimal)
    fit <- object$spec$transform$transform(object$model$fitted, lambda = object$parmatrix[parameters == "lambda"]$optimal)
    res <- act - fit
    sigma.res <- sd(res, na.rm = TRUE)
    
    if (!is.null(innov)) {
        requiredn <- h * nsim
        if (length(innov) != requiredn) {
            stop("\nlength of innov must be nsim x h")
        }
        # check that the innovations are uniform samples (from a copula)
        if (innov_type == "q") {
            if (any(innov < 0 | innov > 1 )) {
                stop("\ninnov must be >0 and <1 (uniform samples) for innov_type = 'q'")
            }
            if (any(innov == 0)) innov[which(innov == 0)] <- 1e-12
            if (any(innov == 1)) innov[which(innov == 1)] <- (1 - 1e-12)
            innov <- matrix(innov, nrow = nsim, ncol = h)
            E <- qnorm(innov, sd = sigma.res)
        } else {
            E <- matrix(innov * sigma.res, nrow = nsim, ncol = h)
        }
    } else {
        E <- matrix(rnorm(h * nsim, 0, sigma.res), ncol = h, nrow = nsim)
    }
    xseed <- tail(object$model$states, 1)
    if (!is.null(init_states)) {
        if (length(as.vector(init_states)) != ncol(xseed)) {
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
    if (!is.null(sigma_scale)) {
        ysim <- f$simulated[,-1,drop = FALSE]
        if (exact_moments) {
            pmoments <- tsmoments.tsissm.estimate(object, h = h, newxreg = newxreg, init_states = xseed)
            for (i in 1:ncol(ysim)) {
                ysim[,i] <- scale(ysim[,i])
                ysim[,i] <- (ysim[,i]*sqrt(pmoments$variance[i])) + pmoments$mean[i]
            }
        }
        mu <- colMeans(ysim)
        ysim <- sweep(ysim, 2, mu, "-")
        ysim <- sweep(ysim, 2, sigma_scale, "*")
        ysim <- sweep(ysim, 2, mu, "+")
    } else {
        ysim <- f$simulated[,-1,drop = FALSE]
        if (exact_moments) {
            pmoments <- tsmoments.tsissm.estimate(object, h = h, newxreg = newxreg, init_states = xseed)
            for (i in 1:ncol(ysim)) {
                ysim[,i] <- scale(ysim[,i])
                ysim[,i] <- (ysim[,i]*sqrt(pmoments$variance[i])) + pmoments$mean[i]
            }
        }
    }
    ysim <- object$spec$transform$inverse(ysim, lambda = lambda)
    analytic_mean <- tsmoments.tsissm.estimate(object, h = h, newxreg = newxreg, init_states = xseed, transform = TRUE)$mean

    if (NCOL(ysim) == 1) ysim <- matrix(ysim, ncol = 1)
    colnames(ysim) <- as.character(forc_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    states <- f$states
    spec <- object$spec
    spec$parmatrix <- object$parmatrix
    
    colnames(E) <- as.character(forc_dates)
    attr(E, "date_class") <- date_class
    class(E) <- "tsmodel.distribution"
    error <- list(distribution = E, original_series = residuals(object, raw = TRUE))
    class(error) <- "tsmodel.predict"

    zList <- list(original_series = xts(object$spec$target$y_orig, object$spec$target$index),
                  distribution = ysim, mean = zoo(colMeans(ysim), forc_dates), 
                  analytic_mean = zoo(analytic_mean, forc_dates), h = h,
                  spec = spec, states = states, 
                  decomp = tsdecompose(object),
                  xseed = xseed,
                  innov = error,
                  sigma = sigma.res,
                  frequency = ifelse(is.null(object$spec$seasonal$seasonal_frequency[1]),1,object$spec$seasonal$seasonal_frequency[1]))
    class(zList) <- c("tsissm.predict","tsmodel.predict")
    return(zList)
}


tsmoments.tsissm.estimate <- function(object, h, newxreg = NULL, init_states = NULL, transform = FALSE, ...)
{
    mu <- rep(0, h)
    sig2 <- rep(0, h)
    sigma_r <- sd(residuals(object, raw = TRUE), na.rm = T)
    w <- t(object$model$w)
    Fmat <- object$model$F
    g <- object$model$g
    if (!is.null(init_states)) {
        x <- t(init_states)
    } else {
        x <- t(tail(object$model$states,1))
    }
    cj <- 0
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    if (object$spec$xreg$include_xreg) {
        beta <- object$parmatrix[grepl("kappa",parameters)]$optimal
        X <- as.numeric(newxreg %*% beta)
    } else {
        X <- rep(0, h)
    }
    for (i in 1:h) {
        Fmat_power <- matpower(Fmat,i - 1)
        mu[i] <- w %*% (Fmat_power %*% x) + X[i]
        if (i == 1) {
            sig2[i] <- sigma_r^2
        } else {
            Fmat_power <- matpower(Fmat,i - 2)
            cj <- cj + (w %*% (Fmat_power %*% g))^2
            sig2[i] <- (sigma_r^2)*(1 + cj)
        }
        if (transform) {
            if (object$spec$transform$name == "box-cox") {
                if (lambda == 0) {
                    mu[i] <- exp(mu[i]) * (1 + sig2[i]/2)
                } else {
                    mu[i] <- (lambda * mu[i] + 1)^(1/lambda) * (1 + (sig2[i] * (1 - lambda))/(2 * (lambda * mu[i] + 1)^2))
                }
            } else {
                mu[i] <- object$spec$transform$inverse(mu[i])
            }
        }
    }
    return(list(mean = mu, variance = sig2))
}
