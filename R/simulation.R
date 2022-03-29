simulate.tsissm.estimate <- function(object, nsim = 1, seed = NULL, h = NULL, newxreg = NULL, sim_dates = NULL, bootstrap = FALSE, innov = NULL,
                                     sigma_scale = 1, pars = coef(object), init_states = object$spec$xseed, ...)
{
    parameters <- NULL
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (is.null(h)) h <- object$spec$dims[2]
    if (is.null(init_states)) {
        xseed <- object$spec$xseed
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
        pmatrix[list(names(pars))]$optimal <- pars
        pmatrix <- pmatrix[list(object$parmatrix$parameters)]
        pars <- pmatrix$optimal
        names(pars) <- pmatrix$parameters
        Mnames <- na.omit(object$spec$S$pars)
        S <- object$spec$S
        S[which(!is.na(pars)),"values"] <- pars[Mnames]
    } else {
        pars <- object$parmatrix$optimal
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
        }
        else {
            if (is.null(sim_dates)) {
                sim_dates <- future_dates(head(object$spec$target$index, 1), frequency = object$spec$target$sampling, n = h)
            }
            warning("\nxreg use in estimation but newxreg is NULL...setting to zero")
            newxreg <- xts(matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = h), sim_dates)
            colnames(newxreg) <- colnames(object$spec$xreg$xreg)
        }
    }
    act <- object$spec$transform$transform(object$spec$target$y_orig, lambda = object$parmatrix[parameters == "lambda"]$optimal)
    fit <- object$spec$transform$transform(object$model$fitted, lambda = object$parmatrix[parameters == "lambda"]$optimal)
    res <- act - fit
    sigma_res <- sd(res, na.rm = TRUE)
    
    if (bootstrap) {
        res <- na.omit(res) * sigma_scale
        E <- matrix(sample(res, h * nsim, replace = TRUE), ncol = h, nrow = nsim)
    } else {
        if (is.null(innov)) {
            E <- matrix(rnorm(h * nsim, 0, sigma_res * sigma_scale), ncol = h, nrow = nsim)
        } else {
            if (length(innov) != (h * nsim)) {
                stop("\nlength innov must be nsim x h")
            }
            if (any(innov) == 0) {
                innov[which(innov == 0)] <- 1e-12
            }
            if ( any(innov == 1)) {
                innov[which(innov == 1)] <- 1 - 1e-12
            }
            innov <- matrix(innov, h, nsim)
            E <- qnorm(innov) * (sigma_res * sigma_scale)
        }
    }
    ##################################################################
    mdim <- c(object$spec$dims[1], nsim, h, ncol(object$spec$xreg$xreg))
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
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    ysim <- matrix(object$spec$transform$inverse(f$simulated[,-1], lambda = lambda), nrow = nsim, ncol = h)
    colnames(ysim) <- as.character(sim_dates)
    class(ysim) <- "tsmodel.distribution"
    attr(ysim, "date_class") <- date_class
    L <- .tsdecompose.simulation(object, f$states, sim_dates, date_class)
    components <- names(L)
    zList <- list(Simulated = ysim, Error = E, dates = as.character(sim_dates),
                  seed = RNGstate, pars = object$parmatrix[estimate == 1]$optimal,
                  sigma = sigma_res, sigma_scale = sigma_scale, components = components)
    zList <- c(zList, L)
    class(zList) <- c("tsissm.simulate")
    return(zList)
}


.tsdecompose.simulation <- function(object, states, sim_dates, date_class, ...)
{
    estimate <- NULL
    idx <- object$spec$idmatrix
    rnames <- rownames(idx)
    indx <- object$spec$target$index
    pars <- object$parmatrix[estimate == 1]$optimal
    names(pars) <- object$spec$parmatrix[estimate == 1]$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    sim_dates <- as.character(sim_dates)
    ##################################################################
    mdim <- object$spec$dims
    nsim <- dim(states)[[3]]
    L <- list()
    k <- 1
    w <- matrix(S[list("w")]$values, nrow = mdim[1], ncol = 1)
    w_t <- t(w)
    s_names <- c()
    if (idx["Level","n"] > 0) {
        Level <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(states[-1, idx["Level","start"]:idx["Level","end"], i],nrow = 1)
        }))
        colnames(Level) <- sim_dates
        class(Level) <- "tsmodel.distribution"
        attr(Level, "date_class") <- date_class
        L[[1]] <- Level
        k <- k + 1
        s_names <- c(s_names,"Level")
    }
    if (idx["Slope","n"] > 0) {
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(w_t[,nstart:nend] * states[-1, idx["Slope","start"]:idx["Slope","end"], i], nrow = 1)
        }))
        colnames(Slope) <- sim_dates
        class(Slope) <- "tsmodel.distribution"
        attr(Slope, "date_class") <- date_class
        L[[k]] <- Slope
        k <- k + 1
        s_names <- c(s_names,"Slope")
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
                colnames(tmp) <- sim_dates
                class(tmp) <- "tsmodel.distribution"
                attr(tmp, "date_class") <- date_class
                L[[k]] <- tmp
                k <- k + 1
                s_names <- c(s_names,paste0("Seasonal",object$spec$seasonal$seasonal_frequency[j]))
            }
        } else {
            j <- 1
            for (i in 1:ns) {
                tmp <- do.call(rbind, lapply(1:nsim, function(i){
                    matrix(t(w_t[,nstart[j]:(nstart[j] + frequency[j] - 1)] %*% t(states[-1, nstart[j]:(nstart[j] + frequency[j] - 1), i])), nrow = 1)
                }))
                colnames(tmp) <- sim_dates
                class(tmp) <- "tsmodel.distribution"
                attr(tmp, "date_class") <- date_class
                L[[k]] <- tmp
                k <- k + 1
                s_names <- c(s_names,paste0("Seasonal",object$spec$seasonal$seasonal_frequency[j]))
            }
        }
    }
    if (idx["AR","n"] > 0) {
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[-1, idx["AR","start"]:idx["AR","end"], i])), nrow = 1)
        }))
        colnames(tmp) <- sim_dates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        L[[k]] <- tmp
        k <- k + 1
        s_names <- c(s_names,paste0("AR",object$spec$arma$order[1]))

    }
    if (idx["MA","n"] > 0) {
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[-1, idx["MA","start"]:idx["MA","end"], i])), nrow = 1)
        }))
        colnames(tmp) <- sim_dates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        L[[k]] <- tmp
        k <- k + 1
        s_names <- c(s_names,paste0("MA",object$spec$arma$order[2]))
    }
    names(L) <- s_names
    return(L)
}
