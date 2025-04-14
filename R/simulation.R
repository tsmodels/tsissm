#' Model Simulation
#'
#' @description Simulation function for class \dQuote{tsissm.estimate}.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param nsim the number of paths per complete set of time steps (h).
#' @param seed a value specifying if and how the random number generator should be
#' initialized (\sQuote{seeded}). Either NULL or an integer that will be used in
#' a call to set.seed before simulating the response vectors.
#' @param h the number of time steps to simulate paths for. If this is NULL,
#' it will use the same number of periods as in the original series.
#' @param newxreg an optional matrix of regressors to use for the simulation
#' if xreg was used in the estimation. If NULL and the estimated object had
#' regressors, and h was also set to NULL, then the original regressors will
#' be used.
#' @param sim_dates an optional vector of simulation dates equal to h. If NULL
#' will use the implied periodicity of the data to generate a regular sequence
#' of dates after the first available date in the data.
#' @param bootstrap whether to bootstrap the innovations from the estimated
#' object by re-sampling from the empirical distribution.
#' @param init_states An optional vector of states to initialize the forecast. 
#' If NULL, will use the first available states from the estimated model.
#' @param init_res For a dynamic variance model, the initialization for the
#' ARCH recursion of length equal to max(p,q).
#' @param init_sigma For a dynamic variance model, the standard deviation initialization 
#' for the GARCH recursion of length equal to max(p,q).
#' @param innov an optional vector of innovations (see innov_type). 
#' The length of this vector should be equal to nsim x horizon.
#' @param innov_type if \sQuote{innov} is not NULL, then this denotes the type of values
#' passed, with \dQuote{q} denoting quantile probabilities (default and
#' backwards compatible) and \dQuote{z} for standardized errors.
#' @param pars an optional named vector of model coefficients which override the
#' estimated coefficients. No checking is currently performed on the adequacy of
#' these coefficients.
#' @param ... not currently used.
#' @return An object of class \dQuote{tsissm.simulate} with slots for the simulated
#' series and states.
#' @aliases simulate
#' @method simulate tsissm.estimate
#' @rdname simulate
#' @export
#'
#'
simulate.tsissm.estimate <- function(object, nsim = 1, seed = NULL, h = 1, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, innov = NULL,
                                     innov_type = "q", pars = coef(object), init_states = tail(object$model$states,1), init_res = NULL, 
                                     init_sigma = NULL, ...)
{
    # check for negative values in simulation and exclude + re-simulate
    switch(object$spec$variance$type, 
           "constant" = .simulate_constant(object = object, nsim = nsim, h = h, newxreg = newxreg,
                                          sim_dates = sim_dates, bootstrap = bootstrap, innov = innov, 
                                          innov_type = innov_type, pars = pars, 
                                          init_states = init_states, seed = seed, ...),
           "dynamic" = .simulate_dynamic(object = object, nsim = nsim, h = h, newxreg = newxreg,
                                         sim_dates = sim_dates, bootstrap = bootstrap, innov = innov, 
                                         innov_type = innov_type, pars = pars, 
                                         init_states = init_states, init_res = init_res, 
                                         init_sigma = init_sigma, seed = seed, ...)
    )
}


.simulate_constant <- function(object, nsim = 1, h = 1, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, innov = NULL,
                               innov_type = "q", pars = coef(object), init_states = tail(object$model$states,1), seed = NULL, ...)
{
    group <- parameters <- NULL
    if (!is.null(seed)) set.seed(seed)
    if (is.null(h)) h <- object$spec$dims[2]
    
    if (is.null(init_states)) {
        row_sample <- sample(1:nrow(object$model$states), 1)
        xseed <- object$model$states[row_sample,, drop = FALSE]
    } else {
        if (length(as.numeric(init_states)) != object$spec$dims[1]) stop(paste0("\ninit.states must be of dimensions 1 x ",object$spec$dims[1]))
        xseed <- init_states
    }
    date_class <- attr(object$spec$target$sampling, "date_class")
    if (!is.null(pars)) {
        estimate <- NULL
        pmatrix <- copy(object$parmatrix)
        setkeyv(pmatrix, "parameters")
        pars <- na.omit(pars[pmatrix$parameters])
        if (length(pars) == 0) {
            stop("\npars names do not match any of the model parameters")
        }
        pmatrix[list(names(pars))]$value <- pars
        pmatrix <- pmatrix[list(object$parmatrix$parameters)]
        pars <- pmatrix$value
        names(pars) <- pmatrix$parameters
        Mnames <- na.omit(object$spec$S$pars)
        S <- object$spec$S
        S[which(!is.na(pars)),"values"] <- pars[Mnames]
    } else {
        pars <- object$parmatrix$value
        names(pars) <- object$parmatrix$parameters
        Mnames <- na.omit(object$spec$S$pars)
        S <- object$spec$S
        S[which(!is.na(pars)),"values"] <- pars[Mnames]
    }
    if (!object$spec$xreg$include_xreg) {
        newxreg <- matrix(0, ncol = 1, nrow = h)
        if (is.null(sim_dates)) {
            sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
        }
    } else {
        if (!is.null(newxreg)) {
            sim_dates <- index(newxreg)
        } else {
            if (is.null(sim_dates)) {
                sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
            newxreg <- xts(matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h), sim_dates)
            colnames(newxreg) <- colnames(object$spec$xreg$xreg)
        }
    }
    newxreg <- rbind(matrix(0, ncol = ncol(newxreg), nrow = 1), newxreg)
    
    res <- object$model$error[-1]
    sigma_res <- object$model$sigma
    skew <- object$parmatrix[parameters == "skew"]$value
    shape <- object$parmatrix[parameters == "shape"]$value
    
    if (bootstrap) {
        res <- na.omit(res)
        E <- matrix(sample(res, h * nsim, replace = TRUE), ncol = h, nrow = nsim)
    } else {
        if (is.null(innov)) {
            r <- rdist(object$spec$distribution$type, nsim * h, 0, sigma_res, skew = skew, shape = shape)
            E <- matrix(r, nrow = nsim, ncol = h)
        } else {
            if (innov_type == "q") {
                if (length(innov) != (h * nsim)) {
                    stop("\nlength innov must be nsim x h")
                }
                if (any(innov) == 0) {
                    innov[which(innov == 0)] <- 1e-12
                }
                if ( any(innov == 1)) {
                    innov[which(innov == 1)] <- 1 - 1e-12
                }
                innov <- matrix(innov, nrow = nsim, ncol = h)
                E <- matrix(qdist(object$spec$distribution$type, innov, 0, sigma = sigma_res, skew = skew, shape = shape), nrow = nsim, ncol = h)
            } else {
                E <- matrix(innov * sigma_res, nrow = nsim, ncol = h)
            }
        }
    }
    E <- cbind(matrix(0, ncol = 1, nrow = nsim), E)
    

    ##################################################################
    n_states <- object$spec$dims[1]
    mdim <- as.integer(c(n_states, nsim, h))
    F0 <- matrix(S[list("F0")]$values, n_states, n_states)
    F1 <- matrix(S[list("F1")]$values, n_states, n_states)
    F2 <- matrix(S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    kappa <- as.numeric(S[list("kappa")]$values)
    f <- .isspredict_constant(F0 = F0, F1 = F1, F2 = F2, w = w, g = g, E = E, xseed = xseed, X = newxreg, kappa = kappa, mdim = mdim)

    lambda <- object$parmatrix[parameters == "lambda"]$value
    ysim <- matrix(object$spec$transform$inverse(f$simulated[,-1], lambda = lambda), nrow = nsim, ncol = h)
    colnames(ysim) <- as.character(sim_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    xspec <- object$spec
    xspec$parmatrix <- object$parmatrix
    zList <- list(distribution = ysim, Error = E, dates = as.character(sim_dates),
                  spec = xspec,
                  seed = seed, pars = pars, sigma = sigma_res, states = f$states)
    class(zList) <- c("tsissm.simulate")
    return(zList)
}


.simulate_dynamic <- function(object, nsim = 1, seed = NULL, h = NULL, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, innov = NULL, 
                              innov_type = "q", pars = coef(object), init_states = tail(object$model$states,1), 
                              init_res = NULL, init_sigma = NULL, ...)
{
    group <- parameters <- NULL
    if (!is.null(seed)) set.seed(seed)
    if (is.null(h)) h <- object$spec$dims[2]
    if (is.null(init_states)) {
        row_sample <- sample(1:nrow(object$model$states), 1)
        xseed <- object$model$states[row_sample,, drop = FALSE]
    } else {
        if (length(as.numeric(init_states)) != object$spec$dims[1]) stop(paste0("\ninit.states must be of dimensions 1 x ",object$spec$dims[1]))
        xseed <- init_states
    }
    date_class <- attr(object$spec$target$sampling, "date_class")
    if (!is.null(pars)) {
        estimate <- NULL
        pmatrix <- copy(object$parmatrix)
        setkeyv(pmatrix, "parameters")
        pars <- na.omit(pars[pmatrix$parameters])
        if (length(pars) == 0) {
            stop("\npars names do not match any of the model parameters")
        }
        pmatrix[list(names(pars))]$value <- pars
        pmatrix <- pmatrix[list(object$parmatrix$parameters)]
        pars <- pmatrix$value
        names(pars) <- pmatrix$parameters
        Mnames <- na.omit(object$spec$S$pars)
        S <- object$spec$S
        S[which(!is.na(pars)),"values"] <- pars[Mnames]
    } else {
        pars <- object$parmatrix$value
        names(pars) <- object$parmatrix$parameters
        Mnames <- na.omit(object$spec$S$pars)
        S <- object$spec$S
        S[which(!is.na(pars)),"values"] <- pars[Mnames]
    }
    if (!object$spec$xreg$include_xreg) {
        newxreg <- matrix(0, ncol = 1, nrow = h)
        if (is.null(sim_dates)) {
            sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
        }
    } else {
        if (!is.null(newxreg)) {
            sim_dates <- index(newxreg)
        } else {
            if (is.null(sim_dates)) {
                sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
            newxreg <- xts(matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h), sim_dates)
            colnames(newxreg) <- colnames(object$spec$xreg$xreg)
        }
    }
    newxreg <- rbind(matrix(0, ncol = ncol(newxreg), nrow = 1), newxreg)
    
    res <- object$model$error[-1]
    garch_sigma <- object$model$sigma[-1]
    omega <- object$model$target_omega
    z_res <- res/garch_sigma
    skew <- object$parmatrix[parameters == "skew"]$value
    shape <- object$parmatrix[parameters == "shape"]$value
    
    if (bootstrap) {
        z_res <- na.omit(z_res)
        Z <- matrix(sample(z_res, h * nsim, replace = TRUE), ncol = h, nrow = nsim)
    } else {
        if (is.null(innov)) {
            r <- rdist(object$spec$distribution$type, nsim * h, 0, 1, skew = skew, shape = shape)
            Z <- matrix(r, ncol = h, nrow = nsim)
        } else {
            if (innov_type == "q") {
                if (length(innov) != (h * nsim)) {
                    stop("\nlength innov must be nsim x h")
                }
                if (any(innov) == 0) {
                    innov[which(innov == 0)] <- 1e-12
                }
                if ( any(innov == 1)) {
                    innov[which(innov == 1)] <- 1 - 1e-12
                }
                innov <- matrix(innov, nrow = nsim, ncol = h)
                Z <- matrix(qdist(object$spec$distribution$type, innov, 0, sigma = 1, skew = skew, shape = shape), nrow = nsim, ncol = h)
            } else {
                Z <- matrix(innov, nrow = nsim, ncol = h)
            }
        }
    }
    Z <- cbind(matrix(0, ncol = 1, nrow = nsim), Z)
    
    maxpq <- max(object$spec$variance$order)
    E <- rep(0, h + maxpq)
    if (!is.null(init_res)) {
        if (length(init_res) != maxpq) stop(paste0("\ninit_res must be of length ", maxpq))
        E[1:maxpq] <- init_res
    } else {
        init_res <- tail(object$model$error, maxpq)
        E[1:maxpq] <- init_res
    }
    
    if (is.null(init_sigma)) {
        init_garch <- tail(object$model$sigma^2, max(object$spec$variance$order))
    } else {
        if (length(init_sigma) != maxpq) stop(paste0("\ninit_sigma must be of length ", maxpq))
        init_garch <- init_sigma^2
    }
    
    variance_intercept <- object$model$target_omega
    init_arch <- init_res^2

    ##################################################################
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
                             Z = Z, xseed = xseed, X = newxreg, kappa = kappa, 
                             init_arch = init_arch, init_garch = init_garch,
                             alpha = alpha, beta = beta, 
                             variance_intercept = variance_intercept,
                             mdim = mdim, order = as.integer(object$spec$variance$order))
    lambda <- object$parmatrix[parameters == "lambda"]$value
    ysim <- matrix(object$spec$transform$inverse(f$simulated[,-1], lambda = lambda), nrow = nsim, ncol = h)
    colnames(ysim) <- as.character(sim_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    xspec <- object$spec
    xspec$parmatrix <- object$parmatrix
    
    E <- f$Error
    E <- E[,-1]
    colnames(E) <- as.character(sim_dates)
    attr(E, "date_class") <- date_class
    class(E) <- "tsmodel.distribution"
    
    sigma <- f$sigma[,-1]
    colnames(sigma) <- as.character(sim_dates)
    attr(sigma, "date_class") <- date_class
    class(sigma) <- "tsmodel.distribution"
    
    zList <- list(distribution = ysim, Error = E, sigma = sigma, dates = as.character(sim_dates),
                  spec = xspec, seed = seed, pars = pars, states = f$states)
    class(zList) <- c("tsissm.simulate")
    return(zList)
}


#' @method simulate tsissm.selection
#' @rdname simulate
#' @export
#'
#'
simulate.tsissm.selection <- function(object, nsim = 1, seed = NULL, h = 1, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, 
                                      pars = coef(object), init_states = tail(object$model$states,1), 
                                      init_res = NULL, init_sigma = NULL, ...)
{
    C <- object$correlation
    cop <- normalCopula(P2p(C), dim = NCOL(C), dispstr = "un")
    U <- rCopula(h * nsim, cop)
    out <- future_lapply(1:length(object$models), function(i){
        r <- matrix(U[,i], nrow = nsim, ncol = h)
        simulate(object$models[[i]], h = h, seed = seed, newxreg = newxreg, nsim = nsim, sim_dates = sim_dates, bootstrap = bootstrap, 
                 innov = r, pars = pars, innov_type = "q", init_states = init_states, init_res = init_res, ...)
    }, future.seed = TRUE, future.packages = "tsissm")
    out <- eval(out)
    L <- list(models = out, table = object$table, correlation = object$correlation)
    class(out) <- c("tsissm.selection_predict")
    return(out)
}