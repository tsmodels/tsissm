#hresiduals
hresiduals.tsissm.estimate <- function(object, h = 12, cores = 1, seed = NULL, trace = FALSE, raw = TRUE, index_start = 1, ...)
{
    n <- length(object$spec$target$y)
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (trace) {
        iterations <- n
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    if (index_start <= 0 | index_start >= n) stop("\nindex_start must be strictly positive and less than length of dataset.")
    i <- 0
    object$model$states <- rbind(object$model$xseed, object$model$states)
    h_residuals <- foreach(i = index_start:(n - 1), .packages = c("tsmethods","tsissm","xts","data.table"), .options.snow = opts) %dopar% {
        tmp <- tsissm:::hresiduals.issm(object, h = h, nsim = 2000, start = i, seed = seed, raw = raw)
        return(tmp)
    }
    if (trace == 1) {
        close(pb)
    }
    stopCluster(cl)
    h_residuals <- rbindlist(h_residuals)
    h_residuals <- dcast(h_residuals, date~horizon, value.var = "error")
    return(h_residuals)
}

hresiduals.issm <- function(object,  h = 12, nsim = 1000, start = 1, seed = NULL, raw = TRUE, ...)
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
    # truncate h to the available time periods
    N <- length(object$spec$target$y)
    if ((start + h) > N) {
        h <- N - start
    }
    if (raw) {
        actual <- object$spec$target$y[(start):(start + h - 1)]
    } else {
        actual <- object$spec$target$y_orig[(start):(start + h - 1)]
    }

    sigma.res <- sd(residuals(object, h = 1, raw = TRUE))
    E <- matrix(rnorm(h * nsim, 0, sigma.res), ncol = h, nrow = nsim)

    xseed <- object$model$states[start, ,drop = FALSE]
    pars <- object$parmatrix$optimal
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    if (!object$spec$xreg$include_xreg) {
        xreg = matrix(0, ncol = 1, nrow = h)
    } else {
        xreg <- object$spec$xreg$x[start:(start + h)]
    }
    ##################################################################
    mdim = c(object$spec$dims[1], nsim, h, ncol(object$spec$xreg$xreg))
    f <- isspredict(f0_ = S[list("F0")]$values,
                    f1_ = S[list("F1")]$values,
                    f2_ = S[list("F2")]$values,
                    w_ = as.numeric(S[list("w")]$values),
                    g_ = as.numeric(S[list("g")]$values),
                    error_ = as.vector(E),
                    xseed_ = xseed,
                    xreg_ = as.numeric(xreg),
                    kappa_ = S[list("kappa")]$values,
                    mdim = mdim)
    lambda <- object$parmatrix[parameters == "lambda"]$optimal
    if (!raw) {
        Y <- matrix(object$spec$transform$inverse(f$simulated[,-1], lambda = lambda), ncol = h, nrow = nsim)
    } else {
        Y <- f$simulated[,-1, drop = FALSE]
    }
    mu <- colMeans(Y)
    e <- actual - mu
    return(data.table(date = object$spec$target$index[start], error = e, horizon = 1:h))
}
