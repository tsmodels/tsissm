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
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param sigma_scale a scaling factor for the innovations standard deviation.
#' @param solver choice of solver to use for the estimation of the paths.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsissmad package.
#' @param ... not currently used.
#' @note The function can use parallel functionality as long as the user has set
#' up a \code{\link[future]{plan}} using the future package.
#' @return An object of class \dQuote{tsissm.profile}.
#' @aliases tsprofile
#' @method tsprofile tsissm.estimate
#' @rdname tsprofile
#' @export
#'
tsprofile.tsissm.estimate <- function(object, h = 1, nsim = 100, seed = NULL, trace = FALSE, sigma_scale = 1, solver = "nlminb", autodiff = TRUE, ...)
{
    if (object$spec$seasonal$include_seasonal & object$spec$seasonal$seasonal_type == "regular") {
        if (autodiff) {
            autodiff <- FALSE
            warning("\nautodiff only currently supported for trigonometric seasonality (switching to non autodiff)")
        }
    }
    sim <- simulate(object, seed = seed, nsim = nsim, h = length(object$spec$target$y_orig) + h, sigma_scale = sigma_scale)
    profile <- profile_fun(sim$Simulated, object, h, trace = trace, solver = solver, autodiff = autodiff)
    profile$sigma <- sim$sigma * sim$sigma_scale
    class(profile) <- "tsissm.profile"
    return(profile)
}

profile_fun <- function(sim, object, h, trace, solver, autodiff)
{

    if (trace) {
        prog_trace <- progressor(nrow(sim))
    }
    date_class <- attr(object$spec$target$sampling, "date_class")
    date_fun <- match.fun(paste0("as.",date_class))
    prof %<-% future_lapply(1:nrow(sim), function(i) {
        if (trace) prog_trace()
        parameters <- NULL
        y <- xts(sim[i,], date_fun(colnames(sim)))
        yin <- y[1:(nrow(y) - h)]
        spec <- tsspec(object, yin, lambda = object$parmatrix[parameters == "lambda"]$optimal)
        # add try catch
        mod <- try(estimate(spec, solver = solver, autodiff = autodiff), silent = TRUE)
        if (inherits(mod, 'try-error')) {
            return(list(L1 = NULL, L2 = NULL))
        }
        p <- predict(mod, h = h)
        L1 <- data.table("Variable" = names(coef(mod)), "Value" = coef(mod), "Simulation" = i)
        L2 <- data.table("Predicted" = as.numeric(p$mean), "Actual" = as.numeric(tail(y, h)), "Simulation" = i, "Horizon" = 1:h)
        return(list(L1 = L1, L2 = L2))
    }, future.packages = c("tsmethods","tsissm","xts","data.table"), future.seed = TRUE)
    prof <- eval(prof)
    C <- rbindlist(lapply(1:length(prof), function(i) prof[[i]]$L1))
    M <- rbindlist(lapply(1:length(prof), function(i) prof[[i]]$L2))

    Actual <- NULL
    Predicted <- NULL
    Simulation <- NULL
    # create distribution for all performance metrics
    maped <- M[,list(MAPE = mape(Actual, Predicted)), by = c("Horizon","Simulation")]
    maped <- dcast(maped, Simulation~Horizon, value.var = "MAPE")
    maped[,Simulation := NULL]
    maped <- as.matrix(maped)
    class(maped) <- "tsmodel.distribution"
    attr(maped, "date_class") <- "numeric"

    biasd <- M[,list(Bias = bias(Actual, Predicted)), by = c("Horizon","Simulation")]
    biasd <- dcast(biasd, Simulation~Horizon, value.var = "Bias")
    biasd[,Simulation := NULL]
    biasd <- as.matrix(biasd)
    class(biasd) <- "tsmodel.distribution"
    attr(biasd, "date_class") <- "numeric"

    mslred <- M[,list(MSLRE = mslre(Actual, Predicted)), by = c("Horizon","Simulation")]
    mslred <- dcast(mslred, Simulation~Horizon, value.var = "MSLRE")
    mslred[,Simulation := NULL]
    mslred <- as.matrix(mslred)
    class(mslred) <- "tsmodel.distribution"
    attr(mslred, "date_class") <- "numeric"

    L <- list(MAPE = maped, BIAS = biasd, MSLRE = mslred, coef = C, true.coef = coef(object))
    return(L)
}
