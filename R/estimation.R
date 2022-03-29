estimate.tsissm.spec <- function(object, solver = "nlminb", control = list(trace = 0), autodiff = FALSE, ...)
{
    # initialize parameters
    if (object$seasonal$include_seasonal & object$seasonal$seasonal_type == "regular") {
        if (autodiff) {
            autodiff <- FALSE
            warning("\nautodiff only currently supported for trigonometric seasonality (switching to non autodiff)")
        }
    }
    estimate <- NULL
    tic <- Sys.time()
    solver_env = new.env(hash = TRUE)
    assign("issm_llh", 1, envir = solver_env)
    object$solver_env <- solver_env
    pars <- object$parmatrix[estimate == 1]$initial
    assign("issm_llh", 1, envir = solver_env)
    if (autodiff) {
        opt <- estimate_ad(object, solver = solver, control = control, ...)
        object$xseed <- opt$xseed
        f <- iss_filter(opt$pars, obj = object)
        parameters <- NULL
        object$target$y <- object$transform$transform(object$target$y_orig, f$parmatrix[parameters == "lambda"]$optimal)
        f$spec <- object
        f$opt <- opt$solver_out
        f$opt$elapsed <- difftime(Sys.time(), tic, units = "mins")
        f$hessian <- opt$hessian
        f$autodiff <- TRUE
        class(f) <- "tsissm.estimate"
    } else {
        if (solver == "nlminb") {
            if (is.null(control$iter.max)) control$iter.max <- 5000
            if (is.null(control$eval.max)) control$eval.max <- 5000
            control2 <- list()
            control2$maxit <- 100
            control2$trace <- 0
            opt <- optim(par = pars, control = control2, method = "Nelder-Mead", fn = loglik_fun, obj = object)
            opt <- nlminb(start = opt$par, control = control, objective = loglik_fun, lower = object$parmatrix[estimate == 1]$lower,
                          upper = object$parmatrix[estimate == 1]$upper, obj = object)
            pars <- opt$par
        } else if (solver == "optim") {
            if (is.null(control$maxit)) control$maxit <- 5000
            control$parscale <- object$parmatrix[estimate == 1]$scale
            pars <- transform_pars(pars, lower =  object$parmatrix[estimate == 1]$lower, upper = object$parmatrix[estimate == 1]$upper, inverse = TRUE)
            opt <- optim(par = pars, control = control, method = "BFGS", fn = loglik_fun_unc, obj = object)
            pars <- opt$par
            pars <- transform_pars(pars, object$parmatrix[estimate == 1]$lower, object$parmatrix[estimate == 1]$upper)
        } else if (solver == "solnp") {
            control2 <- list()
            control2$maxit <- 100
            control2$trace <- 0
            opt <- optim(par = pars, control = control2, method = "Nelder-Mead", fn = loglik_fun, obj = object)
            opt <- solnp(pars = opt$par, fun = loglik_fun, LB = object$parmatrix[estimate == 1]$lower,
                         UB = object$parmatrix[estimate == 1]$upper, obj = object, control = control)
            pars <- opt$pars
        }
        object$xseed <- get("xseed", solver_env)
        if (is.null(object$xseed)) {
            if (solver == "optim" | solver == "solnp") {
                opt <- nlminb(start = pars, control = control, objective = loglik_fun, lower = object$parmatrix[estimate == 1]$lower,
                              upper = object$parmatrix[estimate == 1]$upper, obj = object)
                pars <- opt$par
            } else {
                control2 <- list()
                control2$maxit <- 5000
                control2$parscale <- object$parmatrix[estimate == 1]$scale
                opt <- optim(par = pars, control = control2, method = "Nelder-Mead", fn = loglik_fun, obj = object)
                pars <- opt$par
            }
            if (is.null(object$xseed)) {
                stop("\nproblem estimating model...exiting")
            }
        }
        f <- iss_filter(pars, obj = object)
        parameters <- NULL
        object$target$y <- object$transform$transform(object$target$y_orig, f$parmatrix[parameters == "lambda"]$optimal)
        f$spec <- object
        f$opt <- opt
        f$opt$elapsed <- difftime(Sys.time(), tic, units = "mins")
        f$autodiff <- FALSE
        class(f) <- "tsissm.estimate"
    }
    return(f)
}

iss_filter <- function(pars, obj)
{
    estimate <- NULL
    obj$parmatrix[estimate == 1, "initial"] <- pars
    pars <- obj$parmatrix$initial
    names(pars) <- obj$parmatrix$parameters
    Mnames <- na.omit(obj$S$pars)
    S <- obj$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    parnames <- names(pars)
    mdim = obj$dims
    f <- issfilter(f0_ = S[list("F0")]$values,
                      f1_ = S[list("F1")]$values,
                      f2_ = S[list("F2")]$values,
                      w_ = as.numeric(S[list("w")]$values),
                      g_ = as.numeric(S[list("g")]$values),
                      y_ = as.numeric(obj$target$y),
                      lambda_ = as.numeric(pars["lambda"]),
                      xreg_ = S[list("xreg")]$values,
                      kappa_ = S[list("kappa")]$values,
                      mdim = mdim, xseed_ = obj$xseed, good_ = as.numeric(obj$good))
    L <- list()
    L$fitted <- obj$transform$inverse(f$fitted[-1], lambda = pars["lambda"])
    L$residuals <- obj$target$y_orig - L$fitted
    L$states <- f$states[-1,,drop = FALSE]
    L$xseed <- f$xseed
    L$loglik <- f$loglik
    L$condition <- f$condition
    L$w <- f$w
    L$g <- f$g
    L$F <- f$F
    L$D <- f$D
    p <- obj$parmatrix
    colnames(p)[2] <- "optimal"
    return(list(model = L, parmatrix = p))
}

loglik_fun <- function(pars, obj)
{
    estimate <- NULL
    pp <- NULL
    pp <<- pars
    if (any(is.na(pars))) return(Inf)
    #print(pars)
    solver_env <- obj$solver_env
    obj$parmatrix[estimate == 1, "initial"] <- pars
    pars <- obj$parmatrix$initial
    names(pars) <- obj$parmatrix$parameters
    Mnames <- na.omit(obj$S$pars)
    S <- obj$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    parnames <- names(pars)
    if (sum(obj$arma$order) > 0) {
        if (arma_conditions(pars, parnames)) {
            llh <- get("issm_llh", solver_env) + 0.25*(abs(get("issm_llh", solver_env)))
            return(llh)
        }
    }
    mdim = obj$dims
    f <- issestimation(f0_ = S[list("F0")]$values,
                      f1_ = S[list("F1")]$values,
                      f2_ = S[list("F2")]$values,
                      w_ = as.numeric(S[list("w")]$values),
                      g_ = as.numeric(S[list("g")]$values),
                      y_ = as.numeric(obj$target$y),
                      lambda_ = as.numeric(pars["lambda"]),
                      xreg_ = S[list("xreg")]$values,
                      kappa_ = S[list("kappa")]$values,
                      mdim = mdim, good_ = as.numeric(obj$good))

    if (!is.finite(f$loglik) | f$condition == 1 ) {
        llh <- get("issm_llh", solver_env) + 0.25*(abs(get("issm_llh", solver_env)))
    } else{
        llh <- f$loglik
        assign("issm_llh", llh, envir = solver_env)
    }
    assign("xseed", f$xseed, envir = solver_env)
    return(llh)
}

loglik_fun_unc <- function(pars, obj)
{
    estimate <- NULL
    pars <- transform_pars(pars, obj$parmatrix[estimate == 1]$lower, obj$parmatrix[estimate == 1]$upper)
    pp <- NULL
    pp <<- pars
    if (any(is.na(pars))) return(Inf)
    #print(pars)
    solver_env <- obj$solver_env
    obj$parmatrix[estimate == 1, "initial"] <- pars
    pars <- obj$parmatrix$initial
    names(pars) <- obj$parmatrix$parameters
    Mnames <- na.omit(obj$S$pars)
    S <- obj$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    parnames <- names(pars)
    if (sum(obj$arma$order) > 0) {
        if (arma_conditions(pars, parnames)) {
            llh <- get("issm_llh", solver_env) + 0.25*(abs(get("issm_llh", solver_env)))
            return(llh)
        }
    }
    mdim = obj$dims
    f <- issestimation(f0_ = S[list("F0")]$values,
                       f1_ = S[list("F1")]$values,
                       f2_ = S[list("F2")]$values,
                       w_ = as.numeric(S[list("w")]$values),
                       g_ = as.numeric(S[list("g")]$values),
                       y_ = as.numeric(obj$target$y),
                       lambda_ = as.numeric(pars["lambda"]),
                       xreg_ = S[list("xreg")]$values,
                       kappa_ = S[list("kappa")]$values,
                       mdim = mdim, good_ = as.numeric(obj$good))
    
    if (!is.finite(f$loglik) | f$condition == 1 ) {
        llh <- get("issm_llh", solver_env) + 0.25*(abs(get("issm_llh", solver_env)))
    } else{
        llh <- f$loglik
        assign("issm_llh", llh, envir = solver_env)
    }
    assign("xseed", f$xseed, envir = solver_env)
    return(llh)
}

arma_conditions <- function(pars, parnames) {
    if (any(grepl("theta",parnames))) {
        ar <- pars[grepl("theta",parnames)]
        if (all(ar == 0)) {

        } else {
            if (min(Mod(polyroot(c(1, -ar)))) < 1.01) {
                return(TRUE)
            }
        }
    }
    if (any(grepl("psi",parnames))) {
        ma <- pars[grepl("psi",parnames)]
        if (all(ma == 0)) {

        } else {
            if (min(Mod(polyroot(c(1, ma)))) < 1.01) {
                return(TRUE)
            }
        }
    }
    return(FALSE)
}

transform_pars <- function(pars, lower, upper, inverse = FALSE)
{
    if (!inverse) {
        ans <- lower + (upper - lower)/(1 + exp(-1 * pars))
    } else{
        ans <- -1 * log(-(upper - pars)/(-pars + lower))
    }
    return(ans)
}

