#' Model Simulation Based Profiling
#'
#' @description Profiling of model dynamics using simulation/estimation/prediction.
#' @details The function profiles an estimated model by simulating and then estimating
#' multiple paths from the assumed DGP while leaving h values out for prediction
#' evaluation. Each simulated path is equal to the size of the original dataset
#' plus h additional values, and initialized with the initial state vector from
#' the model. The resulting output contains the distribution of the MAPE,
#' percent bias (BIAS) and mean squared log relative error (MSLRE) per horizon h.
#' Since these matrices are of class \dQuote{tsmodel.distribution} they can be
#' readily plotted with the special purpose \dQuote{plot} function for this class
#' from the \dQuote{tsmethods} package. Additionally, a data.table matrix is
#' also returned with the distribution of the coefficients from each path estimation.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param h the forecast horizon on which to evaluate performance metrics.
#' @param nsim the number of paths to generate.
#' @param seed an object specifying if and how the random number generator
#' should be initialized. See the simulate documentation for more details.
#' @param solver only \dQuote{nlortr} currently available.
#' @param trace whether to show the progress bar and additionally output verbose messages.
#' @param control optional control parameters.
#' The user is expected to have set up appropriate handlers for this using 
#' the \code{\link[progressr]{handlers}} function from the \dQuote{progressr} package.
#' @param ... not currently used.
#' @note The function can use parallel functionality as long as the user has set
#' up a \code{\link[future]{plan}} using the future package.
#' The simulated states are checked for positivity and any paths which have negative
#' values are excluded. If more than half of the paths have negative values then
#' an error is raised and the function will stop.
#' @return An object of class \dQuote{tsissm.profile}.
#' @aliases tsprofile
#' @method tsprofile tsissm.estimate
#' @rdname tsprofile
#' @export
#'
tsprofile.tsissm.estimate <- function(object, h = 1, nsim = 100, seed = NULL, solver = "nloptr", control = NULL, trace = FALSE, ...)
{
    solver <- "nloptr"
    if (is.null(control)) {
        control <- issm_control(solver = "nloptr", algorithm = "SLSQP", trace = 0)
    }
    sim <- simulate(object, seed = seed, nsim = nsim, h = length(object$spec$target$y_orig) + h)
    sim <- .check_positivity(sim)
    if (is.null(sim)) return(sim)
    profile <- profile_fun(sim$distribution, object, h = h, solver = solver, control = control, trace = trace)
    class(profile) <- "tsissm.profile"
    return(profile)
}

profile_fun <- function(sim, object, h, solver = "nloptr", control, trace)
{
    Horizon <- Simulation <- NULL
    if (trace) {
        prog_trace <- progressor(nrow(sim))
    }
    tic <- Sys.time()
    date_class <- attr(object$spec$target$sampling, "date_class")
    date_fun <- match.fun(paste0("as.",date_class))
    prof <- future_lapply(1:nrow(sim), function(i) {
        if (trace) prog_trace()
        parameters <- NULL
        y <- xts(sim[i,], date_fun(colnames(sim)))
        yin <- y[1:(nrow(y) - h)]
        spec <- tsspec(object, yin)
        # add try catch
        mod <- try(estimate(spec, solver = solver, control = control, scores = FALSE), silent = TRUE)
        if (inherits(mod, 'try-error')) {
            return(list(L1 = NULL, L2 = NULL, L3 = NULL))
        }
        if (mod$status < 0) {
            return(list(L1 = NULL, L2 = NULL, L3 = NULL))
        }
        p <- predict(mod, h = h)
        L1 <- data.table("Variable" = names(coef(mod)), "Value" = coef(mod), "Simulation" = i)
        L2 <- data.table("Predicted" = as.numeric(p$mean), "Actual" = as.numeric(tail(y, h)), "Simulation" = i, "Horizon" = 1:h)
        L3 <- p
        # L3 <- p$distribution
        return(list(L1 = L1, L2 = L2, L3 = L3))
    }, future.packages = c("tsmethods","tsissm","xts","data.table","tsaux"), future.seed = TRUE)
    prof <- eval(prof)
    toc <- Sys.time()
    dtime <- difftime(toc, tic, units = "secs")
    if (trace) {
        cat(paste0("\nCompleted Profiling in ", round(as.numeric(dtime),2)," secs."))
        cat(paste0("\nCompiling Performance Metrics..."))
    }
    C <- rbindlist(lapply(1:length(prof), function(i) prof[[i]]$L1))
    M <- rbindlist(lapply(1:length(prof), function(i) prof[[i]]$L2))
    
    # Metrics
    Actual <- Predicted <- Simulation <- NULL
    MAPE <- M[,list(MAPE = mape(Actual, Predicted)), by = c("Simulation","Horizon")]
    BIAS <- M[,list(BIAS = bias(Actual, Predicted)), by = c("Simulation","Horizon")]
    MASE <- lapply(1:length(prof), function(i){
        ans <- M[Simulation == i, list(MASE = mase(Actual, Predicted,  prof[[i]]$L3$original_series, frequency = object$spec$seasonal$seasonal_frequency[1])), by = "Horizon"]
        ans[,Simulation := as.integer(i)]
        ans <- ans[,list(Simulation, Horizon, MASE)]
        return(ans)
    })
    MASE <- rbindlist(MASE)
    CRPS <- lapply(1:length(prof), function(i){
        d <- lapply(1:h, function(j){
            ans <- crps(M[Simulation == i & Horizon %in% (1:j)]$Actual, prof[[i]]$L3$distribution[,1:j,drop = FALSE])
            data.table(Simulation = as.integer(i), Horizon = as.integer(j), CRPS = ans)
        })
        rbindlist(d)
    })
    CRPS <- rbindlist(CRPS)
    M <- M[,list(Simulation, Horizon, Actual, Predicted)]
    L <- list(coef = C, true_coef = coef(object), Predictions = M, MAPE = MAPE, MASE = MASE, CRPS = CRPS)
    class(L) <- "tsissm.profile"
    return(L)
}


.check_positivity <- function(sim) {
    if (any(sim$distribution < 0)) {
        exc <- sort(unique(which(sim$distribution < 0, arr.ind = TRUE)[,1]))
        n <- length(exc)
        if (n > (NROW(sim$distribution)/2)) {
            warning("\nmore than half of the simulated paths have negative values. Exiting and returning NULL.")
            return(NULL)
        }
        # exclude negative paths
        sim$distribution <- sim$distribution[-exc, , drop = FALSE]
        sim$Error <- sim$Error[-exc, , drop = FALSE]
        sim$states <- sim$states[-exc]
        return(sim)
    } else {
        return(sim)
    }
}


