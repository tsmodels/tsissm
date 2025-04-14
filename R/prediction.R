#' Model Prediction
#'
#' @description Prediction function for class \dQuote{tsissm.estimate} or \dQuote{tsissm.selection}.
#' @param object an object of class \dQuote{tsissm.estimate} or \dQuote{tsissm.selection}.
#' @param h the forecast horizon.
#' @param newxreg a matrix of external regressors in the forecast horizon.
#' @param nsim the number of simulations to use for generating the simulated
#' predictive distribution.
#' @param forc_dates an optional vector of forecast dates equal to h. If NULL will
#' use the implied periodicity of the data to generate a regular sequence of
#' dates after the last available date in the data.
#' @param innov an optional vector of innovations (see innov_type). 
#' The length of this vector should be equal to nsim x horizon.
#' @param innov_type if \sQuote{innov} is not NULL, then this denotes the type of values
#' passed, with \dQuote{q} denoting quantile probabilities (default and
#' backwards compatible) and \dQuote{z} for standardized errors.
#' @param init_states an optional vector of states to initialize the forecast
#' and override the initial state vector. If NULL, will use the last available 
#' state from the estimated model.
#' @param seed an object specifying if and how the random number generator should be
#' initialized (\sQuote{seeded}). Either NULL or an integer that will be used in
#' a call to set.seed before simulating the response vectors.
#' @param ... not currently used.
#' @details Like all models in the tsmodels framework, prediction is done by
#' simulating h-steps ahead in order to build a predictive distribution.
#' @return An object of class \dQuote{tsissm.predict} which also inherits
#' \dQuote{tsmodel.predict}, with slots for the simulated prediction distribution,
#' the original series (as a zoo object), the original specification object and
#' the mean forecast. The predictive distribution is back transformed if lambda was
#' not set to NULL in the specification. If the input class is \dQuote{tsissm.selection}
#' then the returned class with be \dQuote{tsissm.selection_predict} which hold the
#' list of predicted objects. For this dispatch method, custom innovations are not allowed
#' since correlated innovations are passed to each model predicted (using a Normal Copula)
#' in order to enable ensembling of the simulated predictive distributions.
#' @aliases predict
#' @method predict tsissm.estimate
#' @rdname predict
#' @export
#'
#'
predict.tsissm.estimate <- function(object, h = 12, seed = NULL, newxreg = NULL, nsim = 1000, 
                                    forc_dates = NULL, innov = NULL, innov_type = "q", 
                                    init_states = NULL,
                                    ...)
{
    switch(object$spec$variance$type, 
           "constant" = .predict_constant(object = object, h = h, seed = seed, newxreg = newxreg, nsim = nsim, 
                                                      forc_dates = forc_dates, innov = innov, innov_type = innov_type, 
                                                      init_states = init_states, ...),
           "dynamic" = .predict_dynamic(object = object, h = h, seed = seed, newxreg = newxreg, nsim = nsim, 
                                          forc_dates = forc_dates, innov = innov, innov_type = innov_type, 
                                          init_states = init_states, ...)
           )
    
}

.predict_constant <- function(object, h = 12, seed = NULL, newxreg = NULL, nsim = 1000, 
                              forc_dates = NULL, innov = NULL, innov_type = "q", 
                              init_states = NULL, ...)
{
    if (nsim < 100) stop("\nnsim must be at least 100.")
    if (!is.null(seed)) set.seed(seed)
    group <- parameters <- NULL
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
    newxreg <- rbind(matrix(0, ncol = ncol(newxreg), nrow = 1), coredata(newxreg))
    
    sigma_res <- sigma(object)
    skew <- object$parmatrix[parameters == "skew"]$value
    shape <- object$parmatrix[parameters == "shape"]$value
    
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
            # norm, std, jsu
            E <- matrix(qdist(object$spec$distribution$type, innov, 0, sigma = sigma_res, skew = skew, shape = shape), nrow = nsim, ncol = h)
        } else {
            # allows to ovveride the standardized distribution
            E <- matrix(innov * sigma_res, nrow = nsim, ncol = h)
        }
    } else {
        # norm, std, jsu
        r <- rdist(object$spec$distribution$type, nsim * h, 0, sigma_res, skew = skew, shape = shape)
        E <- matrix(r, ncol = h, nrow = nsim)
    }
    E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
    
    xseed <- tail(object$model$states, 1)
    if (!is.null(init_states)) {
        if (length(as.vector(init_states)) != ncol(xseed)) {
            stop(paste0("\ninit_states must be a vector of length ", ncol(xseed)))
        } else {
            xseed <- matrix(as.numeric(init_states), nrow = 1, ncol = ncol(xseed))
        }
    }
    pars <- object$parmatrix$value
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    n_states <- object$spec$dims[1]
    mdim <- as.integer(c(n_states, nsim, h))
    F0 <- matrix(S[list("F0")]$values, n_states, n_states)
    F1 <- matrix(S[list("F1")]$values, n_states, n_states)
    F2 <- matrix(S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    kappa <- as.numeric(S[list("kappa")]$values)
    f <- .isspredict_constant(F0 = F0, F1 = F1, F2 = F2, 
                    w = w, g = g, E = E, xseed = xseed,
                    X = newxreg, kappa = kappa, mdim = mdim)
    date_class <- attr(object$spec$target$sampling, "date_class")
    
    lambda <- object$parmatrix[parameters == "lambda"]$value
    ysim <- f$simulated[,-1,drop = FALSE]
    ysim <- object$spec$transform$inverse(ysim, lambda = lambda)
    analytic_moments <- tsmoments.tsissm.estimate(object, h = h, newxreg = newxreg[-1,,drop = FALSE], init_states = xseed, transform = TRUE)
    
    if (NCOL(ysim) == 1) ysim <- matrix(ysim, ncol = 1)
    colnames(ysim) <- as.character(forc_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    states <- array(0, dim = c(h + 1, n_states, nsim))
    for (i in 1:nsim) {
        states[,,i] <- f$states[[i]]
    }
    xspec <- object$spec
    xspec$parmatrix <- object$parmatrix
    E <- E[,-1, drop = FALSE]
    colnames(E) <- as.character(forc_dates)
    attr(E, "date_class") <- date_class
    class(E) <- "tsmodel.distribution"
    error <- list(distribution = E, original_series = residuals(object, transformed = TRUE))
    class(error) <- "tsmodel.predict"
    
    original_series <- xts(object$spec$target$y_orig, object$spec$target$index)
    original_states <- tsdecompose(object)
    fr <- ifelse(is.null(object$spec$seasonal$seasonal_frequency[1]),1,object$spec$seasonal$seasonal_frequency[1])
    
    out <- list(original_series = original_series,
                distribution = ysim, mean = zoo(colMeans(ysim), forc_dates), 
                analytic_mean = zoo(analytic_moments$mean, forc_dates), 
                approximate_variance = zoo(analytic_moments$variance, forc_dates), 
                h = h, spec = xspec, original_states = original_states,
                newxreg = newxreg[-1,,drop = FALSE],
                decomp = tsdecompose(object),
                states = states, xseed = xseed, innov = error,
                sigma = sigma_res, frequency = fr)
    class(out) <- c("tsissm.predict","tsmodel.predict")
    return(out)
}


.predict_dynamic <- function(object, h = 12, seed = NULL, newxreg = NULL, nsim = 1000, 
                              forc_dates = NULL, innov = NULL, innov_type = "q", 
                              init_states = NULL, ...)
{
    if (nsim < 100) stop("\nnsim must be at least 100.")
    group <- parameters <- NULL
    if (!is.null(seed)) set.seed(seed)
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
    newxreg <- rbind(matrix(0, ncol = ncol(newxreg), nrow = 1),  coredata(newxreg))
    
    skew <- object$parmatrix[parameters == "skew"]$value
    shape <- object$parmatrix[parameters == "shape"]$value
    maxpq <- max(object$spec$variance$order)
    
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
            # norm, std, jsu
            Z <- matrix(qdist(object$spec$distribution$type, innov, 0, sigma = 1, skew = skew, shape = shape), nrow = nsim, ncol = h)
        } else {
            # allows to ovveride the standardized distribution
            Z <- matrix(innov, nrow = nsim, ncol = h)
        }
    } else {
        # norm, std, jsu
        r <- rdist(object$spec$distribution$type, nsim * h, 0, 1, skew = skew, shape = shape)
        Z <- matrix(r, ncol = h, nrow = nsim)
    }
    Z <- cbind(matrix(0, ncol = maxpq, nrow = nsim), Z)
    # need the pre-residuals
    E <- rep(0, h + maxpq)
    E[1:maxpq] <- tail(object$model$error, maxpq)
    variance_intercept <- object$model$target_omega
    init_arch <- tail(object$model$error^2, maxpq)
    init_garch <- tail(object$model$sigma^2, maxpq)
    xseed <- tail(object$model$states, 1)
    if (!is.null(init_states)) {
        if (length(as.vector(init_states)) != ncol(xseed)) {
            stop(paste0("\ninit_states must be a vector of length ", ncol(xseed)))
        } else {
            xseed <- matrix(as.numeric(init_states), nrow = 1, ncol = ncol(xseed))
        }
    }
    pars <- object$parmatrix$value
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    n_states <- object$spec$dims[1]
    mdim <- as.integer(c(n_states, nsim, h))
    F0 <- matrix(S[list("F0")]$values, n_states, n_states)
    F1 <- matrix(S[list("F1")]$values, n_states, n_states)
    F2 <- matrix(S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    kappa <- as.numeric(S[list("kappa")]$values)
    alpha <- object$parmatrix[group == "eta"]$value
    beta <- object$parmatrix[group == "delta"]$value
    f <- .isspredict_dynamic(F0 = F0, F1 = F1, F2 = F2, w = w, g = g, e = E, 
                             Z = Z, xseed = xseed, X = newxreg, kappa = kappa, init_arch = init_arch, init_garch = init_garch,
                             alpha = alpha, beta = beta, 
                             variance_intercept = variance_intercept,
                             mdim = mdim, order = as.integer(object$spec$variance$order))
    
    date_class <- attr(object$spec$target$sampling, "date_class")
    
    lambda <- object$parmatrix[parameters == "lambda"]$value
    ysim <- f$simulated[,-1,drop = FALSE]
    ysim <- object$spec$transform$inverse(ysim, lambda = lambda)
    analytic_moments <- .tsmoments_dynamic(object, h = h, newxreg = newxreg[-1, , drop = FALSE], init_states = xseed, transform = TRUE)
    if (NCOL(ysim) == 1) ysim <- matrix(ysim, ncol = 1)
    colnames(ysim) <- as.character(forc_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    
    states <- array(0, dim = c(h + 1, n_states, nsim))
    for (i in 1:nsim) {
        states[,,i] <- f$states[[i]]
    }
    xspec <- object$spec
    xspec$parmatrix <- object$parmatrix
    E <- f$Error
    E <- E[,-1, drop = FALSE]
    colnames(E) <- as.character(forc_dates)
    attr(E, "date_class") <- date_class
    class(E) <- "tsmodel.distribution"
    error <- list(distribution = E, original_series = residuals(object, transformed = TRUE))
    class(error) <- "tsmodel.predict"
    
    sigma <- f$sigma[,-1, drop = FALSE]
    colnames(sigma) <- as.character(forc_dates)
    attr(sigma, "date_class") <- date_class
    class(sigma) <- "tsmodel.distribution"
    osig <- xts(object$model$sigma, object$spec$target$index)
    colnames(osig) <- "sigma"
    sigma <- list(distribution = sigma, original_series = osig)
    class(sigma) <- "tsmodel.predict"
    
    original_series <- xts(object$spec$target$y_orig, object$spec$target$index)
    original_states <- tsdecompose(object)
    fr <- ifelse(is.null(object$spec$seasonal$seasonal_frequency[1]),1,object$spec$seasonal$seasonal_frequency[1])
    
    out <- list(original_series = original_series,
                distribution = ysim, mean = zoo(colMeans(ysim), forc_dates), 
                analytic_mean = zoo(analytic_moments$mean, forc_dates), 
                approximate_variance = zoo(analytic_moments$variance, forc_dates),
                analytic_garch = zoo(analytic_moments$garch_variance, forc_dates),
                h = h, spec = xspec, original_states = original_states,
                decomp = tsdecompose(object), xreg = newxreg[-1,,drop = FALSE],
                states = states, xseed = xseed, innov = error,
                sigma = sigma, frequency = fr)
    class(out) <- c("tsissm.predict","tsmodel.predict")
    return(out)
}


#' @method predict tsissm.selection
#' @rdname predict
#' @export
#'
#'
predict.tsissm.selection <- function(object, h = 12, newxreg = NULL, nsim = 1000, 
                                    forc_dates = NULL, init_states = NULL, ...)
{
    C <- object$correlation
    # sample correlated quantiles from a copula
    cop <- normalCopula(P2p(C), dim = NCOL(C), dispstr = "un")
    U <- rCopula(h * nsim, cop)
    out <- future_lapply(1:length(object$models), function(i){
        r <- matrix(U[,i], nrow = nsim, ncol = h)
        predict(object$models[[i]], h = h, newxreg = newxreg, nsim = nsim, forc_dates = forc_dates, innov = r, innov_type = "q",
                init_states = init_states, ...)
    }, future.seed = TRUE, future.packages = c("tsissm","xts"))
    out <- eval(out)
    L <- list(models = out, table = object$table, correlation = object$correlation)
    class(L) <- c("tsissm.selection_predict")
    return(L)
}

#' Analytic Forecast Moments
#'
#' @description Prediction function for class \dQuote{tsissm.estimate}.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param h the forecast horizon.
#' @param newxreg a matrix of external regressors in the forecast horizon.
#' @param init_states optional vector of states to initialize the forecast.
#' If NULL, will use the last available state from the estimated model.
#' @param transform whether to back transform the mean forecast. For the Box-Cox 
#' transformation this uses a Taylor Series expansion to adjust for the bias.
#' @param ... not currently used.
#' @returns a list with a slot for the analytic mean and one for the process variance. In the 
#' case of dynamic variance, it returns an additional slot for the GARCH variance.
#' @details
#' For the constant variance model the conditional moments are given by:
#' \deqn{E\left(y_{n+h|n}^{(\omega)}\right) = \mathbf{w}' \mathbf{F}^{h-1}\mathbf{x}_n}
#' \deqn{
#' V\left(y_{n+h|n}^{(\omega)}\right) = 
#' \begin{cases}
#' \sigma^2, & \text{if } h = 1, \\[6pt]
#' \sigma^2 \left[1 + \sum_{j=1}^{h-1} c_j^2 \right], & \text{if } h \geq 2.
#' \end{cases}}
#' The backtransformed variance should be used with caution as it uses a 
#' higher-order expansion correction which is known to be inaccurate for h > 5. Instead,
#' use the simulated distribution to extract this information.
#' @aliases tsmoments
#' @method tsmoments tsissm.estimate
#' @rdname tsmoments
#' @export
#'
#'
tsmoments.tsissm.estimate <- function(object, h = 1, newxreg = NULL, init_states = NULL, transform = FALSE, ...)
{
    out <- switch(object$spec$variance$type,
                  "constant" = .tsmoments_constant(object, h = h, newxreg = newxreg, init_states = init_states, transform = transform, ...),
                  "dynamic" = .tsmoments_dynamic(object, h = h, newxreg = newxreg, init_states = init_states, transform = transform, ...))
    return(out)
}


.tsmoments_constant <- function(object, h = 1, newxreg = NULL, init_states = NULL, transform = FALSE, ...)
{
    parameters <- NULL
    mu <- rep(0, h)
    sigma_sqr_process <- rep(0, h)
    sigma_r <- sigma(object)
    w <- t(object$model$w)
    Fmat <- object$model$F
    g <- object$model$g
    if (!is.null(init_states)) {
        x <- t(init_states)
    } else {
        x <- t(tail(object$model$states,1))
    }
    cj <- 0
    lambda <- object$parmatrix[parameters == "lambda"]$value
    if (object$spec$xreg$include_xreg) {
        beta <- object$parmatrix[grepl("kappa",parameters)]$value
        X <- as.numeric(newxreg %*% beta)
    } else {
        X <- rep(0, h)
    }
    for (i in 1:h) {
        Fmat_power <- matpower(Fmat,i - 1)
        mu[i] <- w %*% (Fmat_power %*% x) + X[i]
        if (i == 1) {
            sigma_sqr_process[i] <- sigma_r^2
        } else {
            Fmat_power <- matpower(Fmat,i - 2)
            cj <- cj + (w %*% (Fmat_power %*% g))^2
            sigma_sqr_process[i] <- (sigma_r^2)*(1 + cj)
        }
    }
    if (transform) {
        M <- box_cox_bc(mu, sigma_sqr_process, lambda)
        mu <- M$mean
        variance <- M$variance
    } else {
        variance <- sigma_sqr_process
    }
    return(list(mean = mu, variance = variance))
}

.tsmoments_dynamic <- function(object, h = 1, newxreg = NULL, 
                              init_states = NULL, transform = FALSE, ...)
{
    parameters <- NULL
    mu <- rep(0, h)
    sigma_sqr_process <- rep(0, h)
    sigma_garch <- .forecast_garch_sigma(object, h)
    w <- t(object$model$w)
    Fmat <- object$model$F
    g <- object$model$g
    if (!is.null(init_states)) {
        x <- t(init_states)
    } else {
        x <- t(tail(object$model$states,1))
    }
    cj <- 0
    lambda <- object$parmatrix[parameters == "lambda"]$value
    if (object$spec$xreg$include_xreg) {
        beta <- object$parmatrix[grepl("kappa",parameters)]$value
        X <- as.numeric(newxreg %*% beta)
    } else {
        X <- rep(0, h)
    }
    for (i in 1:h) {
        Fmat_power <- matpower(Fmat,i - 1)
        mu[i] <- w %*% (Fmat_power %*% x) + X[i]
        if (i == 1) {
            sigma_sqr_process[i] <- sigma_garch[i]^2
        } else {
            Fmat_power <- matpower(Fmat,i - 2)
            cj <- cj + (w %*% (Fmat_power %*% g))^2
            sigma_sqr_process[i] <- (sigma_garch[i]^2)*(1 + cj)
        }
    }
    if (transform) {
        unc_garch <- object$model$target_omega/(1 - object$model$persistence)
        M <- box_cox_garch_bc(mu, sigma_sqr_process, lambda, sigma_garch^2, unc_garch)
        mu <- M$mean
        variance <- M$variance
    } else {
        variance <- sigma_sqr_process
    }
    return(list(mean = mu, variance = variance, garch_variance = sigma_garch^2))
}

.forecast_garch_sigma <- function(object, h = 1, ...) {
    group <- NULL
    omega <- object$model$target_omega
    alpha <- object$parmatrix[group == "eta"]$value
    beta <- object$parmatrix[group == "delta"]$value
    order <- as.integer(object$spec$variance$order)
    maxpq <- max(order)
    init_arch <- tail(object$model$error^2, maxpq)
    init_garch <- tail(object$model$sigma^2, maxpq)
    epsilon_sqr <- sigma_sqr <- rep(0, h + maxpq)
    sigma_sqr[1:maxpq] <- init_garch
    epsilon_sqr[1:maxpq] <- init_arch
    for (i in (maxpq + 1):(h + maxpq)) {
        sigma_sqr[i] <- omega
        if (order[1] > 0) {
            for (j in 1:order[2]) {
                sigma_sqr[i] <- sigma_sqr[i] + beta[j] * sigma_sqr[i - j]
            }
        }
        if (order[1] > 0) {
            for (j in 1:order[1]) {
                if ((i - maxpq) > j) {
                    s <- sigma_sqr[(i - j)]
                } else {
                    s <- epsilon_sqr[(i - j)]
                }
                sigma_sqr[i] <- sigma_sqr[i] + alpha[j] * s
            }
        }
    }
    sig <- sqrt(sigma_sqr)
    if (maxpq > 0) sig <- sig[-c(1:maxpq)]
    return(sig)
}

box_cox_bc <- function(mu, variance, lambda) {
    if (lambda == 0) {
        # Log-normal case (lambda = 0)
        mu_bc <- exp(mu) * (1 + variance / 2)
        var_bc <- exp(2 * mu) * (exp(variance) - 1)
    } else {
        # Bias-adjusted mean transformation
        mu_bc <- (lambda * mu + 1)^(1/lambda) * (1 + (variance * (1 - lambda)) / (2 * (lambda * mu + 1)^2))
        
        # Base variance transformation (normal assumption)
        var_bc <- variance * (lambda * mu + 1)^(2 * (1 - lambda)/lambda) 
        
        # Higher-order correction (normal assumption)
        var_bc <- var_bc + (variance^2 * (1 - lambda)^2 / (2 * (lambda * mu + 1)^4)) * 
            (lambda * mu + 1)^(4 * (1 - lambda)/lambda - 4)
    }
    
    return(list(mean = mu_bc, variance = var_bc))
}

box_cox_garch_bc <- function(mu, variance, lambda, garch_variance, unc_garch) {
    if (lambda == 0) {
        # Log-normal case (lambda = 0)
        mu_bc <- exp(mu) * (1 + variance / 2)
        var_bc <- exp(2 * mu) * (exp(variance) - 1)
    } else {
        # Bias-adjusted mean transformation
        mu_bc <- (lambda * mu + 1)^(1/lambda) * (1 + (variance * (1 - lambda)) / (2 * (lambda * mu + 1)^2))
        
        # Adjust variance transformation for time-varying variance
        scale_factor <- (garch_variance/unc_garch)
        
        # Base variance transformation
        var_bc <- variance * (lambda * mu + 1)^(2 * (1 - lambda)/lambda) 
        
        # Higher-order correction with GARCH variance adjustment
        var_bc <- var_bc + (variance^2 * (1 - lambda)^2 / (2 * (lambda * mu + 1)^4)) * 
            (lambda * mu + 1)^(4 * (1 - lambda)/lambda - 4)
        
        # Apply scaling correction for GARCH
        var_bc <- var_bc * scale_factor
    }
    
    return(list(mean = mu_bc, variance = var_bc))
}

