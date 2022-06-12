#' Automatic Model Selection
#'
#' @description Automatic model selection based on the AIC criterion.
#' @param y an xts vector.
#' @param slope whether to include or not the slope. Having both TRUE and FALSE 
#' will estimate all possible combinations with and without a slope.
#' @param slope_damped whether to include or not the slope dampening. Having 
#' both TRUE and FALSE will estimate all possible combinations with and without 
#' slope dampening.
#' @param seasonal whether to include or not a seasonal component. Having both 
#' TRUE and FALSE  will estimate all possible combinations with and without a 
#' seasonal component.
#' @param seasonal_frequency a vector of seasonal frequencies.
#' @param seasonal_type The type of seasonality to include. Trigonometric is the 
#' preferred type and the only one accepting multiple seasonal frequencies.
#' @param seasonal_harmonics a list with slots for each seasonal frequency listing 
#' the sequence of harmonics to test for each (see details).
#' @param ar a vector of the ar terms to test for (see details).
#' @param ma a vector of the ma terms to test for (see details).
#' @param xreg an optional xts matrix of regressors (pre-lagged).
#' @param transformation The transformation to use (defaults to box-cox and 
#' is effectively NULL if lambda is NULL or 1).
#' @param lambda the box-cox lambda parameter.
#' @param lower the lower bound for the transformation.
#' @param upper the upper bound for the transformation.
#' @param sampling the sampling frequency of the data.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param return_table Whether to return the table with the enumerated options 
#' and the AIC and MAPE for each combination of those options used.
#' @param solver the solver to use for estimation.
#' @param autodiff whether to use automatic differentiation
#' (see \code{\link{estimate.tsissm.spec}}).
#' @param ... not used.
#' @return An object of class \dQuote{tsissm.estimate} returning the best model 
#' based on AIC (minimum).
#' @details The user is responsible for passing reasonable options for the options. 
#' For instance, the harmonics must be strictly positive and less than one half 
#' the seasonal frequency. The \code{\link{expand.grid}} function is used to 
#' enumerate all possible combinations of the options with some sanity checks 
#' and eliminations if testing for both seasonal and non-seasonal models, 
#' or for the case of slope and dampening. If the user prefers to ensemble 
#' models or use instead the MAPE criterion, the \sQuote{return_table} should 
#' be set to TRUE and the returned table used to estimate an alternative to the 
#' top AIC model.
#' @note The function can use parallel functionality as long as the user has set up a
#' \code{\link[future]{plan}} using the future package.
#' @aliases auto_issm
#' @rdname auto_issm
#' @export
#'
auto_issm <- function(y, slope = c(TRUE,FALSE), slope_damped = c(TRUE,FALSE), seasonal = c(TRUE,FALSE), 
                      seasonal_frequency = 1, seasonal_type = "trigonometric", seasonal_harmonics = list(),
                      ar = 0:2, ma = 0:2, xreg = NULL, transformation = "box-cox", lambda = 1, 
                      lower = 0, upper = 1, sampling = NULL, trace = FALSE, 
                      return_table = FALSE, solver = "nlminb", autodiff = TRUE, ...)
{
    if (any(seasonal)) {
        if (!all(seasonal_frequency > 1)) stop("seasonal_frequency must be > 1 if any seasonal is set to TRUE")
        nf <- length(seasonal_frequency)
        if (!is.list(seasonal_harmonics)) stop("\nseasonal_harmonics must be a list equal to the length of seasonal_frequency")
        if (length(seasonal_harmonics) != nf) stop("\nseasonal_harmonics must be a list equal to the length of seasonal_frequency")
        # check that harmonics < seasonal_frequency/2
        for (i in 1:length(seasonal_harmonics)) if (any(seasonal_harmonics[[i]] >= 0.5 * seasonal_frequency[i])) stop("\nseasonal_harmonics must be less then 1/2 of seasonal_frequency")
        sk <- paste0("K",1:length(seasonal_frequency))
        names(seasonal_harmonics) <- sk
    } else {
        seasonal_frequency <- NULL
        seasonal_harmonics <- NULL
    }
    args_grid_nonseasonal <- NULL
    args_grid_seasonal <- NULL
    if (!all(seasonal)) {
        args_list <- list()
        args_list$slope <- unique(as.logical(slope))
        args_list$slope_damped <- unique(as.logical(slope_damped))
        args_list$seasonal <- FALSE
        args_list$ar <- ar
        args_list$ma <- ma
        args_grid_nonseasonal <- as.data.frame(expand.grid(args_list))
    }
    if (any(seasonal)) {
        args_list <- list()
        args_list$slope <- unique(as.logical(slope))
        args_list$slope_damped <- unique(as.logical(slope_damped))
        args_list$seasonal <- TRUE
        args_list$ar <- ar
        args_list$ma <- ma
        for (i in 1:length(seasonal_harmonics)) args_list <- c(args_list, seasonal_harmonics[i])
        args_grid_seasonal <- as.data.frame(expand.grid(args_list))
        if (!is.null(args_grid_nonseasonal)) {
            tmp <- matrix(0, ncol = length(sk), nrow = nrow(args_grid_nonseasonal))
            colnames(tmp) <- sk
            tmp <- as.data.frame(tmp)
            args_grid_nonseasonal <- cbind(args_grid_nonseasonal, tmp)
        }
    }
    args_grid <- rbind(args_grid_nonseasonal, args_grid_seasonal)
    if (any(args_grid$slope)) {
        if (any(args_grid$slope_damped)) {
            exc <- which(args_grid$slope == FALSE & args_grid$slope_damped)
            args_grid <- args_grid[-exc,]
        }
    }
    n <- NROW(args_grid)
    gnames <- colnames(args_grid)
    if (trace) {
        prog_trace <- progressor(n)
    }
    b %<-% future_lapply(1:n, function(i) {
        if (trace) prog_trace()
        iter <- NULL
        spec <- issm_modelspec(y, slope = args_grid[i,"slope"], slope_damped = args_grid[i,"slope_damped"], seasonal = args_grid[i,"seasonal"], 
                               seasonal_frequency = seasonal_frequency, seasonal_harmonics = as.numeric(args_grid[i,grepl("K",gnames)]), 
                               seasonal_type = seasonal_type, ar = args_grid[i,"ar"], ma = args_grid[i,"ma"], xreg = xreg, 
                               transformation = transformation, lambda = lambda, lower = lower, upper = upper, sampling = sampling)
        mod <- try(estimate(spec, solver = solver, autodiff = autodiff), silent = TRUE)
        if (inherits(mod, 'try-error')) {
            tab <- data.table(iter = i, AIC = as.numeric(NA), MAPE = as.numeric(NA))
        } else {
            tab <- data.table(iter = i, AIC = AIC(mod), MAPE = mape(y, fitted(mod)))
        }
        return(tab)
    }, future.packages = c("tsmethods","tsaux","xts","tsissm","data.table"))
    b <- eval(b)
    b <- rbindlist(b)
    iter <- NULL
    args_grid <- as.data.table(args_grid)
    args_grid[,iter := 1:.N]
    b <- merge(args_grid, b, by = "iter")
    b <- b[!is.na(AIC)]
    b <- b[order(AIC)]
    # estimate once more the best model
    use <- as.data.frame(b[1])
    use <- use[,gnames]
    spec <- issm_modelspec(y, slope = use[1,"slope"], slope_damped = use[1,"slope_damped"], seasonal = use[1,"seasonal"], 
                           seasonal_frequency = seasonal_frequency, seasonal_harmonics = as.numeric(use[1,grepl("^K[0-9]",gnames)]), 
                           seasonal_type = seasonal_type, ar = use[1,"ar"], ma = use[1,"ma"], xreg = xreg, 
                           transformation = transformation, lambda = lambda, lower = lower, upper = upper, sampling = sampling)
    mod <- estimate(spec, solver = solver, autodiff = autodiff)
    # return table? (merge with results)
    if (return_table) {
        mod$selection <- b
    }
    class(mod) <- c(class(mod), "tsissm.select")
    return(mod)
}
