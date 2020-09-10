estimate.tsissm.spec <- function(object, solver = "nlminb", control = list(trace = 0), ...)
{
    # initialize parameters
    estimate <- NULL
    tic <- Sys.time()
    solver_env = new.env(hash = TRUE)
    assign("issm_llh", 1, envir = solver_env)
    object$solver_env <- solver_env
    pars <- object$parmatrix[estimate == 1]$initial
    assign("issm_llh", 1, envir = solver_env)
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
        opt <- optim(par = pars, control = control, method = "Nelder-Mead", fn = loglik_fun, obj = object)
        pars <- opt$par
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
    f <- iss_filter(pars, obj = object)
    parameters <- NULL
    object$target$y <- object$transform$transform(object$target$y_orig, f$parmatrix[parameters == "lambda"]$optimal)
    f$spec <- object
    f$opt <- opt
    f$opt$elapsed <- difftime(Sys.time(), tic, units = "mins")
    class(f) <- "tsissm.estimate"
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
                      y_ = as.numeric(obj$target$y_orig),
                      lambda_ = as.numeric(pars["lambda"]),
                      xreg_ = S[list("xreg")]$values,
                      kappa_ = S[list("kappa")]$values,
                      mdim = mdim, xseed_ = obj$xseed)
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
            llh <- get("issm_llh", solver_env) + 0.1*(abs(get("issm_llh", solver_env)))
            return(llh)
        }
    }
    mdim = obj$dims
    f <- issestimation(f0_ = S[list("F0")]$values,
                      f1_ = S[list("F1")]$values,
                      f2_ = S[list("F2")]$values,
                      w_ = as.numeric(S[list("w")]$values),
                      g_ = as.numeric(S[list("g")]$values),
                      y_ = as.numeric(obj$target$y_orig),
                      lambda_ = as.numeric(pars["lambda"]),
                      xreg_ = S[list("xreg")]$values,
                      kappa_ = S[list("kappa")]$values,
                      mdim = mdim)

    if (!is.finite(f$loglik) | f$condition == 1 ) {
        llh <- get("issm_llh", solver_env) + 0.1*(abs(get("issm_llh", solver_env)))
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
