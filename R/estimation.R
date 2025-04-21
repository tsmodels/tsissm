log_parameters <- function(pars, file = "~/parameters.txt") {
    write.table(matrix(pars, nrow = 1), file = file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}

#' Model Estimation
#'
#' @description Estimates a model given a specification object using
#' maximum likelihood.
#' @param object an object of class \dQuote{tsissm.spec} or \dQuote{tsissm.autospec}.
#' @param control solver control parameters passed to the nloptr function.
#' @param scores whether to calculate the analytic scores (Jacobian) of the
#' likelihood. This is not available for the \dQuote{tsissm.autospec} object.
#' @param trace whether to show a progress bar for the automatic selection object.
#' @param debug_mode for development testing, will include the TMB object.
#' @param ... not used.
#' @details The maximum likelihood estimation for this model is described in the 
#' vignette.
#' @note When calculating the scores, a future promise is created so it is fastest if
#' a future plan is pre-created with at least 2 workers so that the function can 
#' run in the background without having to wait for the estimation object to be
#' returned.\cr
#' If the control argument is NULL, then we hybrid strategy is adopted whereby
#' the SLSQP algorithm is initially used and if it fails (status less than 0) 
#' then the Augmented Lagrange with MMA local solver is used which is slower but
#' may be more reliable.\cr
#' For the automatic selection estimation, this will benefit from the use of multiple
#' processes which can be set up with a \code{\link[future]{plan}}. For progress
#' tracing, use \code{\link[progressr]{handlers}}. The function will check for
#' number of parallel workers initialized (using \sQuote{nbrOfWorkers}), and if 
#' it finds only 1 then will revert to non-parallel execution of the code.
#' @returns An object of class \dQuote{tsissm.estimate} or \dQuote{tsissm.selection}.
#' In the case of automatic model selection an object of class \dQuote{tsissm.estimate} 
#' will be returned based on AIC (minimum) if \dQuote{top_n} is 1, else an object of
#' class \dQuote{tsissm.selection} with a list of length \dQuote{top_n}. This object
#' can then be used for filtering, prediction and simulation and then ensembled 
#' (based on user specified weights).
#' @aliases estimate
#' @method estimate tsissm.spec
#' @rdname estimate
#' @export
#'
#'
estimate.tsissm.spec <- function(object, control = issm_control(algorithm = "SLSQP", trace = 0), scores = TRUE, debug_mode = FALSE, ...)
{
    estimate <- NULL
    tic <- Sys.time()
    solver_env = new.env(hash = TRUE)
    assign("issm_llh", 1, envir = solver_env)
    object$solver_env <- solver_env
    pars <- object$parmatrix[estimate == 1]$initial
    assign("issm_llh", 1, envir = solver_env)
    assign("xseed", 100, envir = solver_env)
    
    if (object$variance$type == "constant") {
        opt <- .estimate_ad_constant(object, control = control, return_scores = scores, ...)
        object$xseed <- opt$xseed
        f <- iss_filter_constant(opt$pars, obj = object)
    } else if (object$variance$type == "dynamic") {
        opt <- .estimate_ad_dynamic(object, control = control, return_scores = scores, ...)
        object$xseed <- opt$xseed
        f <- iss_filter_dynamic(opt$pars, obj = object, opt = opt)
    }
    parameters <- NULL
    # transform will return original data when lambda = 1
    object$target$y <- object$transform$transform(object$target$y_orig, f$parmatrix[parameters == "lambda"]$value)
    
    f$spec <- object
    if (debug_mode) {
        f$tmb <- opt$tmb
        f$opt <- opt$solver_out
    }
    f$status <- opt$solver_out$status
    f$score_promise <- opt$scores
    f$elapsed <- difftime(Sys.time(), tic, units = "mins")
    f$hessian <- opt$hessian
    f$autodiff <- TRUE
    object$solver_env <- NULL
    class(f) <- "tsissm.estimate"
    return(f)
}

#' @method estimate tsissm.autospec
#' @rdname estimate
#' @export
#'
estimate.tsissm.autospec <- function(object, control = NULL, trace = FALSE, ...)
{
    parameters <- NULL
    args_grid <- object$grid
    gnames <- colnames(args_grid)
    if (is.null(control)) {
        control <- issm_control(algorithm = "SLSQP", trace = 0)
    } else {
        control$print_level <- 0
    }
    distribution <- object$distribution
    garch_order <- object$garch_order
    init_garch <- object$init_garch
    sample_n <- object$sample_n
    sampling <- object$sampling
    xreg <- object$xreg
    lambda <- object$lambda
    lower <- object$lower
    upper <- object$upper
    seasonal_frequency <- object$seasonal_frequency
    y <- object$y
    top_n <- object$top_n
    n <- NROW(args_grid)
    n_cores <- nbrOfWorkers()
    if (trace) {
        tmp_tic <- Sys.time()
        tmp_spec <- issm_modelspec(y, slope = args_grid[1,"slope"], slope_damped = args_grid[1,"slope_damped"], seasonal = args_grid[1,"seasonal"], 
                                   seasonal_frequency = seasonal_frequency, seasonal_harmonics = as.numeric(args_grid[1,grepl("^Seasonal",gnames)]), 
                                   ar = args_grid[1,"ar"], ma = args_grid[1,"ma"], xreg = xreg, lambda = lambda, lower = lower, 
                                   upper = upper, sampling = sampling, variance = args_grid[1,"variance"], distribution = distribution)
        tmp_mod <- try(estimate(tmp_spec, control = control, scores = FALSE), silent = TRUE)
        tmp_toc <- Sys.time() - tmp_tic    
        est_time_one <- tmp_toc
        estimated_time <- round(as.numeric(est_time_one/n_cores) * NROW(args_grid), 2) %/% 60
        print(paste0("no. of models to evaluate: ", NROW(args_grid)))
        print(paste0("estimated evaluation time (mins): ", estimated_time))
        prog_trace <- progressor(n)
    }
    tic <- Sys.time()
    if (n_cores <= 1) {
        b <- lapply(1:n, function(i) {
            if (trace) prog_trace()
            iter <- NULL
            spec <- issm_modelspec(y, slope = args_grid[i,"slope"], slope_damped = args_grid[i,"slope_damped"], seasonal = args_grid[i,"seasonal"], 
                                   seasonal_frequency = seasonal_frequency, seasonal_harmonics = as.numeric(args_grid[i,grepl("^Seasonal",gnames)]), 
                                   ar = args_grid[i,"ar"], ma = args_grid[i,"ma"], xreg = xreg, lambda = lambda, lower = lower, 
                                   upper = upper, sampling = sampling, variance = args_grid[i,"variance"], distribution = distribution)
            mod <- try(estimate(spec, control = control, scores = FALSE), silent = TRUE)
            if (inherits(mod, 'try-error')) {
                tab <- data.table(iter = i, lambda = as.numeric(NA), AIC = as.numeric(NA), MAPE = as.numeric(NA))
            } else {
                tab <- data.table(iter = i, lambda = mod$parmatrix[parameters == "lambda"]$value, AIC = AIC(mod), MAPE = mape(y, fitted(mod)))
            }
            out <- list(model = mod, table = tab)
            return(out)
        })
    } else {
        b <- future_lapply(1:n, function(i) {
            if (trace) prog_trace()
            iter <- NULL
            spec <- issm_modelspec(y, slope = args_grid[i,"slope"], slope_damped = args_grid[i,"slope_damped"], seasonal = args_grid[i,"seasonal"], 
                                   seasonal_frequency = seasonal_frequency, seasonal_harmonics = as.numeric(args_grid[i,grepl("^Seasonal",gnames)]), 
                                   ar = args_grid[i,"ar"], ma = args_grid[i,"ma"], xreg = xreg, lambda = lambda, lower = lower, 
                                   upper = upper, sampling = sampling, variance = args_grid[i,"variance"], distribution = distribution)
            mod <- try(estimate(spec, control = control, scores = FALSE), silent = TRUE)
            if (inherits(mod, 'try-error')) {
                tab <- data.table(iter = i, lambda = as.numeric(NA), AIC = as.numeric(NA), MAPE = as.numeric(NA))
            } else {
                tab <- data.table(iter = i, lambda = mod$parmatrix[parameters == "lambda"]$value, AIC = AIC(mod), MAPE = mape(y, fitted(mod)))
            }
            out <- list(model = mod, table = tab)
            return(out)
        }, future.seed = TRUE, future.packages = c("tsmethods","xts","tsissm","data.table"))
        b <- eval(b)
    }
    toc <- Sys.time()
    dtime <- difftime(toc, tic, units = "mins")
    tab <- lapply(b, function(x) x$table)
    tab <- rbindlist(tab)
    iter <- NULL
    args_grid <- as.data.table(args_grid)
    args_grid[,iter := 1:.N]
    tab <- merge(args_grid, tab, by = "iter")
    tab <- tab[!is.na(AIC)]
    tab <- tab[order(AIC)]
    if (top_n == 1) {
        model <- b[[tab[1]$iter]]$model
        return(model)
    } else {
        # one more check since we eliminated NA (unsuccessful solution)
        n <- NROW(tab)
        if (top_n > n) top_n <- n
        if (n == 0) {
            warning("\nno models successfully estimated, returning NULL")
            return(NULL)
        }
        L <- vector(mode = "list", length = top_n)
        for (i in 1:top_n) L[[i]] <- b[[tab[i]$iter]]$model
        e <- do.call(cbind, lapply(1:top_n, function(i) residuals(L[[i]], transformed = TRUE)))
        colnames(e) <- NULL
        C <- cor(e, method = "kendall")
        # convert back
        C <- sin(pi * C / 2)
        # Ensure diagonals are exactly 1
        diag(C) <- 1
        colnames(C) <- rownames(C) <- paste0("model_",1:top_n)
        out <- list(table = tab[1:top_n], models = L, correlation = C, elapsed = as.numeric(toc))
        class(out) <- "tsissm.selection"
        return(out)
    }
}


iss_filter_constant <- function(pars, obj)
{
    estimate <- parameters <- NULL
    obj$parmatrix[estimate == 1, "initial"] <- pars
    pars <- obj$parmatrix$initial
    names(pars) <- obj$parmatrix$parameters
    Mnames <- na.omit(obj$S$pars)
    S <- obj$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    parnames <- names(pars)
    mdim <- obj$dims
    n_states <- mdim[1]
    timesteps <- mdim[2]   
    n_xreg <- NCOL(obj$xreg$xreg)
    f0 <- matrix(S[list("F0")]$values, n_states, n_states)
    f1 <- matrix( S[list("F1")]$values, n_states, n_states)
    f2 <- matrix( S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    y <- as.numeric(obj$target$y)
    xreg <- rbind(matrix(0, nrow = 1, ncol = n_xreg), matrix(S[list("xreg")]$values, timesteps, n_xreg))
    kappa <- as.numeric(S[list("kappa")]$values)
    xseed <- as.numeric(obj$xseed)
    good <- c(0, as.numeric(obj$good))
    
    f <- .issfilter_constant(F0 = f0, F1 = f1, F2 = f2, w = w, g = g, y = c(1.0, y), X = xreg,  kappa = kappa, 
                    xseed = xseed, good = good, mdim = mdim, lambda = as.numeric(pars["lambda"]))
    L <- list()
    L$fitted <- obj$transform$inverse(f$fitted[-1], lambda = pars["lambda"])
    
    L$residuals <- obj$target$y_orig - L$fitted
    L$states <- f$states[-1,,drop = FALSE]
    L$xseed <- f$xseed
    L$error <- f$error[-1]
    good_index <- which(obj$good == 1)
    L$sigma <- sqrt(sum(L$error[good_index]^2)/length(L$error))
    n <- length(obj$target$y_orig)
    ngood <- sum(obj$good)
    if (obj$distribution$type == "norm") {
        nll <- -log(dnorm(L$error[good_index]/L$sigma, 0, 1)/L$sigma) - (pars["lambda"] - 1.0) * log(as.numeric(obj$target$y_orig[good_index]))
    } else if (obj$distribution$type == "std") {
        shape <- obj$parmatrix[parameters == "shape"]$initial
        nll <- -log(ddist("std", L$error[good_index]/L$sigma, 0, 1, shape = shape)/L$sigma) - (pars["lambda"] - 1.0) * log(as.numeric(obj$target$y_orig[good_index]))
    } else if (obj$distribution$type == "jsu") {
        skew <- obj$parmatrix[parameters == "skew"]$initial
        shape <- obj$parmatrix[parameters == "shape"]$initial
        nll <- -log(ddist("jsu", L$error[good_index]/L$sigma, 0, 1, skew = skew, shape = shape)/L$sigma) - (pars["lambda"] - 1.0) * log(as.numeric(obj$target$y_orig[good_index]))
    }
    loglik <- ngood/n * sum(nll)
    L$nll <- nll
    L$loglik <- loglik
    L$condition <- f$condition
    L$w <- f$w
    L$g <- f$g
    L$F <- f$F
    L$D <- f$D
    p <- obj$parmatrix
    colnames(p)[2] <- "value"
    return(list(model = L, parmatrix = p))
}


iss_filter_dynamic <- function(pars, obj, opt)
{
    estimate <- parameters <- group <- NULL
    ipars <- pars
    obj$parmatrix[estimate == 1, "initial"] <- pars
    pars <- obj$parmatrix$initial
    names(pars) <- obj$parmatrix$parameters
    Mnames <- na.omit(obj$S$pars)
    S <- obj$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    parnames <- names(pars)
    mdim <- obj$dims
    n_states <- mdim[1]
    timesteps <- mdim[2]   
    n_xreg <- NCOL(obj$xreg$xreg)
    f0 <- matrix(S[list("F0")]$values, n_states, n_states)
    f1 <- matrix( S[list("F1")]$values, n_states, n_states)
    f2 <- matrix( S[list("F2")]$values, n_states, n_states)
    w <- as.numeric(S[list("w")]$values)
    g <- as.numeric(S[list("g")]$values)
    y <- as.numeric(obj$target$y)
    xreg <- rbind(matrix(0, nrow = 1, ncol = n_xreg), matrix(S[list("xreg")]$values, timesteps, n_xreg))
    kappa <- as.numeric(S[list("kappa")]$values)
    xseed <- as.numeric(obj$xseed)
    good <- c(0, as.numeric(obj$good))
    eta <- obj$parmatrix[group == "eta"]$initial
    delta <- obj$parmatrix[group == "delta"]$initial
    lambda <- obj$parmatrix[group == "lambda"]$initial
    initv <- opt$tmb$report(ipars)$initial_arch
    initial_variance <- initv
    initial_arch <- initv
    omega <- opt$tmb$report(ipars)$target_omega
    vmodel <- opt$tmb$env$data$vmodel
    f <- .issfilter_dynamic(F0 = f0, F1 = f1, F2 = f2, w = w, g = g, y = c(1.0,y), X = xreg,  kappa = kappa, 
                            xseed = xseed, good = good, eta = eta, delta = delta, mdim = as.integer(mdim),
                            initial_arch = initial_arch, initial_variance = initial_variance, 
                            vmodel = as.integer(vmodel), omega = omega, lambda = lambda)
    
    L <- list()
    L$fitted <- obj$transform$inverse(f$fitted[-1], lambda = pars["lambda"])
    
    L$residuals <- obj$target$y_orig - L$fitted
    L$states <- f$states[-1,,drop = FALSE]
    L$xseed <- f$xseed
    L$error <- f$error[-1]
    good_index <- which(obj$good == 1)
    L$sigma <- f$sigma[-1]
    n <- length(obj$target$y_orig)
    ngood <- sum(obj$good)
    std_z <- L$error[good_index]/L$sigma[good_index]
    if (obj$distribution$type == "norm") {
        nll <- -1.0 * log(dnorm(std_z, 0, 1)/L$sigma[good_index]) - (pars["lambda"] - 1.0) * log(as.numeric(obj$target$y_orig[good_index]))
    } else if (obj$distribution$type == "std") {
        shape <- obj$parmatrix[parameters == "shape"]$initial
        nll <- -log(ddist("std", std_z, 0, 1, shape = shape)/L$sigma[good_index]) - (pars["lambda"] - 1.0) * log(as.numeric(obj$target$y_orig[good_index]))
    } else if (obj$distribution$type == "jsu") {
        skew <- obj$parmatrix[parameters == "skew"]$initial
        shape <- obj$parmatrix[parameters == "shape"]$initial
        nll <- -log(ddist("jsu", std_z, 0, 1, skew = skew, shape = shape)/L$sigma[good_index]) - (pars["lambda"] - 1.0) * log(as.numeric(obj$target$y_orig[good_index]))
    }
    loglik <- ngood/n * sum(nll)
    L$nll <- nll
    L$loglik <- loglik
    L$condition <- f$condition
    L$w <- f$w
    L$g <- f$g
    L$F <- f$F
    L$D <- f$D
    L$initial_variance <- initv
    L$initial_arch <- initv
    L$target_omega <- omega
    L$persistence <- opt$tmb$report(ipars)$persistence
    
    p <- obj$parmatrix
    colnames(p)[2] <- "value"
    return(list(model = L, parmatrix = p))
}


tmb_inputs_issm_dynamic <- function(spec)
{
    estimate <- include <- value <- group <- .N <- initial <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[, include := 0]
    parmatrix[estimate == 1, include := 1]
    
    parmatrix[group == "kappa", include := 1]
    parmatrix[group == "kappa" & estimate == 0, include := as.numeric(NA)]
    
    parmatrix[group == "distribution", include := 1]
    parmatrix[group == "distribution" & estimate == 0, include := as.numeric(NA)]
    
    parmatrix[group == "lambda", include := 1]
    parmatrix[group == "lambda" & estimate == 0, include := as.numeric(NA)]
    
    parmatrix[group == "eta", include := 1]
    parmatrix[group == "eta" & estimate == 0, include := as.numeric(NA)]
    
    parmatrix[group == "delta", include := 1]
    parmatrix[group == "delta" & estimate == 0, include := as.numeric(NA)]
    
    
    
    map <- lapply(split(parmatrix[include == 1 | is.na(include), list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    
    parameters <- lapply(split(parmatrix[include == 1 | is.na(include), list(initial, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$initial))
    
    S <- spec$S
    S <- S[matrix %in% c("F0","F1","F2","g","w")]
    S <- S[matrix != "xreg"]
    V <- S$values
    V[which(is.na(V))] <- 0
    findex <- which(!is.na(S$pars)) - 1
    p <- spec$parmatrix
    p <- p[group == "issmpars"]
    pars <- p[estimate == 1]
    allpars <- p$initial
    parnames_all <- p$parameters
    parnames_estimate <- pars$parameters
    ppindex <- match(pars$parameters, p$parameters) - 1
    
    fpindex <- match(na.omit(S$pars), p$parameters) - 1
    f0_index <- range(which(S$matrix == "F0"))
    f0_index <- c(f0_index[1] - 1, f0_index[2] - f0_index[1] + 1)
    f1_index <- range(which(S$matrix == "F1"))
    f1_index <- c(f1_index[1] - 1, f1_index[2] - f1_index[1] + 1)
    f2_index <- range(which(S$matrix == "F2"))
    f2_index <- c(f2_index[1] - 1, f2_index[2] - f2_index[1] + 1)
    modeli <- spec$dims
    windex <-  range(which(S$matrix == "w"))
    windex <- c(windex[1] - 1, windex[2] - windex[1] + 1)
    gindex <-  range(which(S$matrix == "g"))
    gindex <- c(gindex[1] - 1, gindex[2] - gindex[1] + 1)
    fshape <- c(f0_index, f1_index, f2_index, windex, gindex)
    y <- spec$target$y
    X <- spec$xreg$xreg
    
    par_list <- list(issmpars = parmatrix[group == "issmpars" & estimate == 1]$initial, lambda = parmatrix[group == "lambda"]$initial, 
                     kappa = parmatrix[group == "kappa"]$initial, 
                     eta = parmatrix[group == "eta"]$initial,
                     delta = parmatrix[group == "delta"]$initial,
                     distribution = parmatrix[group == "distribution"]$initial)
    
    tmb_names <- rep("pars", length(parnames_estimate))
    lower <- pars$lower
    upper <- pars$upper
    
    if (spec$transform$include_lambda) {
        lower <- c(lower, spec$parmatrix[parameters == "lambda"]$lower)
        upper <- c(upper, spec$parmatrix[parameters == "lambda"]$upper)
        parnames_estimate <- c(parnames_estimate, "lambda")
        tmb_names <- c(tmb_names, "lambda")
    }
    if (spec$xreg$include_xreg) {
        if (any(spec$parmatrix[grepl("kappa",parameters)]$estimate == 1)) {
            lower <- c(lower, spec$parmatrix[group == "kappa" & estimate == 1]$lower)
            upper <- c(upper, spec$parmatrix[group == "kappa" & estimate == 1]$upper)
            parnames_estimate <- c(parnames_estimate, spec$parmatrix[group == "kappa" & estimate == 1]$parameters)
            tmb_names <- c(tmb_names, rep("kappa",nrow(spec$parmatrix[group == "kappa" & estimate == 1])))
        }
    }

    if (spec$parmatrix[group == "eta"]$estimate == 1) {
        lower <- c(lower, spec$parmatrix[group == "eta" & estimate == 1]$lower)
        upper <- c(upper, spec$parmatrix[group == "eta" & estimate == 1]$upper)
        parnames_estimate <- c(parnames_estimate, spec$parmatrix[group == "eta" & estimate == 1]$parameters)
        tmb_names <- c(tmb_names, rep("eta",nrow(spec$parmatrix[group == "eta" & estimate == 1])))
    }
    
    if (spec$parmatrix[group == "delta"]$estimate == 1) {
        lower <- c(lower, spec$parmatrix[group == "delta" & estimate == 1]$lower)
        upper <- c(upper, spec$parmatrix[group == "delta" & estimate == 1]$upper)
        parnames_estimate <- c(parnames_estimate, spec$parmatrix[group == "delta" & estimate == 1]$parameters)
        tmb_names <- c(tmb_names, rep("delta",nrow(spec$parmatrix[group == "delta" & estimate == 1])))
    }
    
    if (spec$distribution$type != "norm") {
        lower <- c(lower, spec$parmatrix[group == "distribution" & estimate == 1]$lower)
        upper <- c(upper, spec$parmatrix[group == "distribution" & estimate == 1]$upper)
        parnames_estimate <- c(parnames_estimate, spec$parmatrix[group == "distribution" & estimate == 1]$parameters)
        tmb_names <- c(tmb_names, rep("distribution",nrow(spec$parmatrix[group == "distribution" & estimate == 1])))
    }
    # check for NA values
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
        y <- na.fill(y, fill = 0)
    }
    valid_index <- which(good > 0) - 1
    
    # create function for ARMA and non ARMA models
    llh_fun <- function(pars, fun, issmenv) {
        names(pars) <- issmenv$tmb_names
        lik <- fun$fn(pars)
        if (is.na(lik) | !is.finite(lik)) {
            lik <- issmenv$lik + 0.1 * abs(issmenv$lik)
            issmenv$lik <- lik
        } else {
            issmenv$lik <- lik
        }
        return(lik)
    }
    grad_fun <- function(pars, fun, issmenv)
    {
        names(pars) <- issmenv$tmb_names
        fun$gr(pars)
    }
    
    hess_fun <- function(pars, fun, issmenv)
    {
        names(pars) <- issmenv$tmb_names
        fun$he(pars, atomic = TRUE)
    }
    
    init_v <- spec$variance$init_variance
    vmodel <- as.integer(c(max(spec$variance$order), spec$variance$order[1], spec$variance$order[2], spec$variance$sample_n))
    
    
    if (spec$parmatrix[parameters == "lambda"]$estimate == 1 ) {
        data <- list(V = V, X = X, good = good, y = y, allpars = allpars, 
                     findex = findex, fpindex = fpindex, ppindex = ppindex, fshape = fshape, 
                     modeli = modeli, initmethod = init_v, vmodel = vmodel, valid_index = valid_index, 
                     model = "garchlambda")
    } else {
        data <- list(V = V, X = X, good = good, y = y, allpars = allpars, 
                     findex = findex, fpindex = fpindex, ppindex = ppindex, fshape = fshape, 
                     modeli = modeli, initmethod = init_v, vmodel = vmodel, valid_index = valid_index, 
                     model = "garch")
        xseed <- .initialize_states_cpp(data, par_list)
        data$xseed <- as.numeric(xseed)
    }
    
    L <- list(data = data, par_list = par_list, map = map, lower = lower, upper = upper, 
              llh_fun = llh_fun, grad_fun = grad_fun, hess_fun = hess_fun, 
              parnames_estimate = parnames_estimate, tmb_names = tmb_names, parnames_all = parnames_all)
    return(L)
}

tmb_inputs_issm_constant <- function(spec)
{
    estimate <- include <- value <- group <- .N <- initial <- NULL
    parmatrix <- spec$parmatrix
    parmatrix[, include := 0]
    parmatrix[estimate == 1, include := 1]
    
    parmatrix[group == "kappa", include := 1]
    parmatrix[group == "kappa" & estimate == 0, include := as.numeric(NA)]
    
    parmatrix[group == "distribution", include := 1]
    parmatrix[group == "distribution" & estimate == 0, include := as.numeric(NA)]
    

    parmatrix[group == "lambda", include := 1]
    parmatrix[group == "lambda" & estimate == 0, include := as.numeric(NA)]
    
    
    map <- lapply(split(parmatrix[include == 1 | is.na(include), list(P = 1:.N * include), by = "group"], by = "group", keep.by = FALSE, drop = T), function(x) as.factor(x$P))
    
    parameters <- lapply(split(parmatrix[include == 1 | is.na(include), list(initial, group)], by = "group", keep.by = FALSE), function(x) as.numeric(x$initial))
    
    S <- spec$S
    S <- S[matrix %in% c("F0","F1","F2","g","w")]
    S <- S[matrix != "xreg"]
    V <- S$values
    V[which(is.na(V))] <- 0
    findex <- which(!is.na(S$pars)) - 1
    p <- spec$parmatrix
    p <- p[group == "issmpars"]
    pars <- p[estimate == 1]
    allpars <- p$initial
    parnames_all <- p$parameters
    parnames_estimate <- pars$parameters
    ppindex <- match(pars$parameters, p$parameters) - 1

    fpindex <- match(na.omit(S$pars), p$parameters) - 1
    f0_index <- range(which(S$matrix == "F0"))
    f0_index <- c(f0_index[1] - 1, f0_index[2] - f0_index[1] + 1)
    f1_index <- range(which(S$matrix == "F1"))
    f1_index <- c(f1_index[1] - 1, f1_index[2] - f1_index[1] + 1)
    f2_index <- range(which(S$matrix == "F2"))
    f2_index <- c(f2_index[1] - 1, f2_index[2] - f2_index[1] + 1)
    modeli <- spec$dims
    windex <-  range(which(S$matrix == "w"))
    windex <- c(windex[1] - 1, windex[2] - windex[1] + 1)
    gindex <-  range(which(S$matrix == "g"))
    gindex <- c(gindex[1] - 1, gindex[2] - gindex[1] + 1)
    fshape <- c(f0_index, f1_index, f2_index, windex, gindex)
    y <- spec$target$y
    X <- spec$xreg$xreg
    
    par_list <- list(issmpars = parmatrix[group == "issmpars" & estimate == 1]$initial, lambda = parmatrix[group == "lambda"]$initial, kappa = parmatrix[group == "kappa"]$initial, distribution = parmatrix[group == "distribution"]$initial)
    tmb_names <- rep("pars", length(parnames_estimate))
    lower <- pars$lower
    upper <- pars$upper
    
    if (spec$transform$include_lambda) {
        lower <- c(lower, spec$parmatrix[parameters == "lambda"]$lower)
        upper <- c(upper, spec$parmatrix[parameters == "lambda"]$upper)
        parnames_estimate <- c(parnames_estimate, "lambda")
        tmb_names <- c(tmb_names, "lambda")
    }
    if (spec$xreg$include_xreg) {
        if (any(spec$parmatrix[grepl("kappa",parameters)]$estimate == 1)) {
            lower <- c(lower, spec$parmatrix[group == "kappa" & estimate == 1]$lower)
            upper <- c(upper, spec$parmatrix[group == "kappa" & estimate == 1]$upper)
            parnames_estimate <- c(parnames_estimate, spec$parmatrix[group == "kappa" & estimate == 1]$parameters)
            tmb_names <- c(tmb_names, rep("kappa",nrow(spec$parmatrix[group == "kappa" & estimate == 1])))
        }
    }
    if (spec$distribution$type != "norm") {
        lower <- c(lower, spec$parmatrix[group == "distribution" & estimate == 1]$lower)
        upper <- c(upper, spec$parmatrix[group == "distribution" & estimate == 1]$upper)
        parnames_estimate <- c(parnames_estimate, spec$parmatrix[group == "distribution" & estimate == 1]$parameters)
        tmb_names <- c(tmb_names, rep("distribution",nrow(spec$parmatrix[group == "distribution" & estimate == 1])))
    }
    # check for NA values
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
        y <- na.fill(y, fill = 0)
    }
    valid_index <- which(good > 0) - 1
    # create function for ARMA and non ARMA models
    llh_fun <- function(pars, fun, issmenv) {
        log_parameters(pars, file = "~/parameters.txt")
        names(pars) <- issmenv$tmb_names
        lik <- fun$fn(pars)
        if (is.na(lik) | !is.finite(lik)) {
            lik <- issmenv$lik + 0.1 * abs(issmenv$lik)
            issmenv$lik <- lik
        } else {
            issmenv$lik <- lik
        }
        return(lik)
    }
    grad_fun <- function(pars, fun, issmenv)
    {
        names(pars) <- issmenv$tmb_names
        fun$gr(pars)
    }
    
    hess_fun <- function(pars, fun, issmenv)
    {
        names(pars) <- issmenv$tmb_names
        fun$he(pars, atomic = TRUE)
    }
    
    if (spec$parmatrix[parameters == "lambda"]$estimate == 1) {
        data <- list(V = V, X = X, good = good, y = y, allpars = allpars, 
                     findex = findex, fpindex = fpindex, ppindex = ppindex, fshape = fshape, 
                     valid_index = valid_index, modeli = modeli, model = "constantlambda")
    } else {
        data <- list(V = V, X = X, good = good, y = y, allpars = allpars, 
                     findex = findex, fpindex = fpindex, ppindex = ppindex, fshape = fshape, 
                     valid_index = valid_index, modeli = modeli, model = "constant")
        xseed <- .initialize_states_cpp(data, par_list)
        data$xseed <- as.numeric(xseed)
    }
    L <- list(data = data, par_list = par_list, map = map, lower = lower, upper = upper, 
              llh_fun = llh_fun, grad_fun = grad_fun, hess_fun = hess_fun, 
              parnames_estimate = parnames_estimate, tmb_names = tmb_names, parnames_all = parnames_all)
    return(L)
}

.estimate_ad_constant <- function(object, control, return_scores = TRUE, ...)
{
    parameters <- group <- NULL
    spec_list <- tmb_inputs_issm_constant(object)
    other_opts <- list(...)
    if (!is.null(other_opts$silent)) {
        silent <- other_opts$silent
    } else {
        silent <- TRUE
    }
    # box cox transpose init_states
    # append to spec_list$data$xseed as matrix
    fun <- try(MakeADFun(data = spec_list$data, parameters = spec_list$par_list, DLL = "tsissm_TMBExports", 
                         map = spec_list$map, trace = FALSE, silent = silent, checkParameterOrder = T), silent = FALSE)
    fun$env$tracemgc <- FALSE
    if (inherits(fun, 'try-error')) {
        stop("\nestimate found an error. Please use non ad version of estimator and contact developer with reproducible example.")
    }
    
    issmenv <- new.env()
    issmenv$lik <- 1
    issmenv$grad <- NULL
    issmenv$parameter_names <- spec_list$parnames_all
    issmenv$parnames_estimate <- spec_list$parnames_estimate
    issmenv$tmb_names <- spec_list$tmb_names
    issmenv$parmatrix <- object$parmatrix
    issmenv$map <- spec_list$map
    issmenv$parameters <- spec_list$par_list
    
    if (object$arma$order[1] > 1 & object$arma$order[2] <= 1) {
        issmenv$arindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("theta",parameters)]$parameters)
        issmenv$arcons <- ar_jconstraint(fun$par, fun, issmenv)
    } else if (object$arma$order[1] <= 1 & object$arma$order[2] > 1) {
        issmenv$maindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("psi",parameters)]$parameters)
        issmenv$macons <- ma_jconstraint(fun$par, fun, issmenv)
    } else if (object$arma$order[1] > 0 & object$arma$order[2] > 0) {
        issmenv$arindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("theta",parameters)]$parameters)
        issmenv$maindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("psi",parameters)]$parameters)
        issmenv$macons <- ma_jconstraint(fun$par, fun, issmenv)
        issmenv$arcons <- ar_jconstraint(fun$par, fun, issmenv)
    }
    case_id <- model_case(object)
    if (case_id[2] > 0) {
        issmenv$constraint <- issm_constraint(fun$par, fun, issmenv)
    }
    cfun <- make_constraint(object, fun, issmenv)
    if (is.null(control)) {
        control <- issm_control(algorithm = "SLSQP", trace = 0)
        sol <- nloptr(x0 = fun$par, eval_f = spec_list$llh_fun, eval_grad_f = spec_list$grad_fun, eval_g_ineq = cfun,
                      lb = spec_list$lower, ub = spec_list$upper, opts = control, fun = fun, issmenv = issmenv)
        if (sol$status < 0) {
            control <- issm_control(algorithm = "AUGLAG/CCSAQ", trace = 0)
            par_iter <- sol$solution
            sol <- nloptr(x0 = par_iter, eval_f = spec_list$llh_fun, eval_grad_f = spec_list$grad_fun, eval_g_ineq = cfun,
                          lb = spec_list$lower, ub = spec_list$upper, opts = control, fun = fun, issmenv = issmenv)
        }
    } else {
        sol <- nloptr(x0 = fun$par, eval_f = spec_list$llh_fun, eval_grad_f = spec_list$grad_fun, eval_g_ineq = cfun,
                      lb = spec_list$lower, ub = spec_list$upper, opts = control, fun = fun, issmenv = issmenv) 
    }
    spec_list$data$xseed <- as.numeric(fun$env$data$xseed)
    pars <- sol$solution
    if (return_scores) {
        scores <- score_function(pars, spec_list)
    } else {
        scores <- NULL
    }
    names(pars) <- issmenv$tmb_names
    llh <- spec_list$llh_fun(pars, fun, issmenv)
    gradient <- spec_list$grad_fun(pars, fun, issmenv)
    hessian <- spec_list$hess_fun(pars, fun, issmenv)
    names(pars) <- issmenv$parnames_estimate
    parmatrix <- object$parmatrix
    parmatrix[parameters %in% issmenv$estimation_names]$initial <- pars
    colnames(gradient) <- issmenv$estimation_names
    colnames(hessian) <- rownames(hessian) <- issmenv$estimation_names
    xseed <- fun$report()$states[1,,drop = FALSE]
    D <- abs(eigen(fun$report(pars)$D, only.values = TRUE, symmetric = FALSE)$values)
    out <- list(pars = pars, llh = llh, gradient = gradient, hessian = hessian, scores = scores, xseed = xseed, solver_out = sol, tmb = fun, D = D)
    return(out)
}


.estimate_ad_dynamic <- function(object, control = list(trace = 0), return_scores = TRUE, ...)
{
    parameters <- group <- NULL
    spec_list <- tmb_inputs_issm_dynamic(object)
    other_opts <- list(...)
    if (!is.null(other_opts$silent)) {
        silent <- other_opts$silent
    } else {
        silent <- TRUE
    }
    fun <- try(MakeADFun(data = spec_list$data, hessian = TRUE, parameters = spec_list$par_list, DLL = "tsissm_TMBExports", 
                         map = spec_list$map, trace = FALSE, silent = silent, checkParameterOrder = F), silent = FALSE)
    fun$env$tracemgc <- FALSE
    
    if (inherits(fun, 'try-error')) {
        stop("\nestimate found an error. Please use non ad version of estimator and contact developer with reproducible example.")
    }

    issmenv <- new.env()
    issmenv$lik <- 1
    issmenv$grad <- NULL
    issmenv$parameter_names <- spec_list$parnames_all
    issmenv$parnames_estimate <- spec_list$parnames_estimate
    issmenv$tmb_names <- spec_list$tmb_names
    issmenv$parmatrix <- object$parmatrix
    garch_index <- NULL
    if (object$variance$order[1] > 0) {
        garch_index <- c(garch_index, which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[group == "eta"]$parameters))
    }
    if (object$variance$order[2] > 0) {
        garch_index <- c(garch_index, which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[group == "delta"]$parameters))
    }
    issmenv$garchindex <- garch_index
    issmenv$map <- spec_list$map
    issmenv$parameters <- spec_list$par_list
    
    if (object$arma$order[1] > 1 & object$arma$order[2] <= 1) {
        issmenv$arindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("theta",parameters)]$parameters)
        issmenv$arcons <- ar_jconstraint(fun$par, fun, issmenv)
    } else if (object$arma$order[1] <= 1 & object$arma$order[2] > 1) {
        issmenv$maindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("psi",parameters)]$parameters)
        issmenv$macons <- ma_jconstraint(fun$par, fun, issmenv)
    } else if (object$arma$order[1] > 0 & object$arma$order[2] > 0) {
        issmenv$arindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("theta",parameters)]$parameters)
        issmenv$maindex <- which(object$parmatrix[estimate == 1]$parameters %in% object$parmatrix[grepl("psi",parameters)]$parameters)
        issmenv$macons <- ma_jconstraint(fun$par, fun, issmenv)
        issmenv$arcons <- ar_jconstraint(fun$par, fun, issmenv)
    }

    case_id <- model_case(object)
    if (case_id[2] > 0) {
        issmenv$constraint <- issm_constraint(fun$par, fun, issmenv)
    }
    cfun <- make_constraint(object, fun, issmenv)
    if (is.null(control)) {
        control <- issm_control(algorithm = "SLSQP", trace = 0)
        sol <- nloptr(x0 = fun$par, eval_f = spec_list$llh_fun, eval_grad_f = spec_list$grad_fun, eval_g_ineq = cfun,
                      lb = spec_list$lower, ub = spec_list$upper, opts = control, fun = fun, issmenv = issmenv)
        if (sol$status < 0) {
            control <- issm_control(algorithm = "AUGLAG/CCSAQ", trace = 0)
            par_iter <- sol$solution
            sol <- nloptr(x0 = par_iter, eval_f = spec_list$llh_fun, eval_grad_f = spec_list$grad_fun, eval_g_ineq = cfun,
                          lb = spec_list$lower, ub = spec_list$upper, opts = control, fun = fun, issmenv = issmenv)
        }
    } else {
        sol <- nloptr(x0 = fun$par, eval_f = spec_list$llh_fun, eval_grad_f = spec_list$grad_fun, eval_g_ineq = cfun,
                      lb = spec_list$lower, ub = spec_list$upper, opts = control, fun = fun, issmenv = issmenv) 
    }
    spec_list$data$xseed <- as.numeric(fun$env$data$xseed)
    pars <- sol$solution
    
    if (return_scores) {
        scores <- score_function(pars, spec_list)
    } else {
        scores <- NULL
    }
    names(pars) <- issmenv$tmb_names
    llh <- spec_list$llh_fun(pars, fun, issmenv)
    gradient <- spec_list$grad_fun(pars, fun, issmenv)
    hessian <- spec_list$hess_fun(pars, fun, issmenv)
    names(pars) <- issmenv$parnames_estimate

    parmatrix <- object$parmatrix
    parmatrix[parameters %in% issmenv$estimation_names]$initial <- pars
    colnames(gradient) <- issmenv$estimation_names
    colnames(hessian) <- rownames(hessian) <- issmenv$estimation_names
    xseed <- fun$report()$states[1,,drop = FALSE]
    D <- abs(eigen(fun$report(pars)$D, only.values = TRUE, symmetric = FALSE)$values)
    out <- list(pars = pars, llh = llh, gradient = gradient, hessian = hessian, scores = scores, xseed = xseed, solver_out = sol, tmb = fun, D = D)
    return(out)
}

score_function <- function(pars, spec_list)
{
    promise <- future({
        fun_scores <- try(MakeADFun(data = spec_list$data, parameters = spec_list$par_list, 
                                    DLL = "tsissm_TMBExports", map = spec_list$map, trace = FALSE, 
                                    silent = TRUE, checkParameterOrder = FALSE, ADreport = TRUE), 
                          silent = FALSE)
        fun_scores$gr(pars)
    }, packages = c("TMB","tsissm"))
    return(list(promise = promise))
}


.initialize_states_cpp <- function(data, parameters, parmatrix = NULL, pars = NULL)
{
    initial <- NULL
    if (!is.null(pars)) {
        group <- NULL
        tmp <- copy(parmatrix)
        tmp[estimate == 1, initial := pars]
        parameters$issmpars <- tmp[group == "issmpars"]$initial
        parameters$kappa <- tmp[group == "kappa"]$initial
        parameters$lambda <- tmp[group == "lambda"]$initial
        parameters$distribution <- tmp[group == "distribution"]$initial
    }
    init_states <- .initialize_states(y = data$y, good = as.numeric(data$good), valid_index = as.integer(data$valid_index), 
                                      issmpars = parameters$issmpars, X = data$X, kappa = parameters$kappa,
                                      V = data$V, allpars = data$allpars, findex = as.integer(data$findex), 
                                      fpindex = as.integer(data$fpindex), ppindex = as.integer(data$ppindex), 
                                      fshape = as.integer(data$fshape), modeli = as.integer(data$modeli), lambda = parameters$lambda)
    return(init_states)
}


seed_transform <- function(x, lambda)
{
    if (lambda < 1e-12) {
        x <- log(x)
    } else {
        x <- (x^lambda - 1) / lambda
    }
    return(x)

}

seed_inverse <- function(x, lambda)
{
    if (lambda < 1e-12) {
        x <- exp(x)
    } else {
        x <- (x * lambda + 1)^(1/lambda)
    }
    return(x)
}