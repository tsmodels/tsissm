#' Model Specification
#'
#' @description Specifies an ISSM model prior to estimation with option for automatic
#' model selection.
#' @details The specification object holds the information and data which is
#' then passed to the maximum likelihood estimation routines. Depending on whether
#' automatic selection is chosen, it will dispatch to the appropriate estimation
#' routine.
#' @param y an xts vector.
#' @param auto whether to use automatic model selection.
#' @param slope (Logical) slope component. If \dQuote{auto} is TRUE, then this can
#' be a vector of size 2 with TRUE and FALSE.
#' @param slope_damped (Logical) slope dampening component. If \dQuote{auto} is TRUE, then this can
#' be a vector of size 2 with TRUE and FALSE.
#' @param seasonal (Logical) seasonal component(s).
#' @param seasonal_frequency vector of numeric seasonal frequencies (can be fractional).
#' @param seasonal_harmonics the number of harmonics per seasonal frequency for 
#' the trigonometric seasonality. If \dQuote{auto} is TRUE, this can be a list 
#' with slots for each seasonal frequency listing the sequence of harmonics to test for each.
#' @param ar AR order.
#' @param ma MA order.
#' @param xreg an xts matrix of external regressors.
#' @param variance either \dQuote{constant} or \dQuote{dynamic}. In the latter
#' case a GARCH model will be used. For the automatic selection case, both can
#' be provided as a vector.
#' @param garch_order the order of the GARCH model (no automatic selection of
#' GARCH order allowed).
#' @param lambda the Box Cox lambda. If not NULL (no transformation), then 
#' either a numeric value or NA denoting automatic estimation.
#' @param lower lower bound for the transformation (defaults to 0).
#' @param upper upper bound for the transformation (defaults to 1.5).
#' @param distribution a choice of the Normal (\dQuote{norm}), Student (\dQuote{std})
#' or Johnson's SU (\dQuote{jsu}) distributions. There is no choice for selecting
#' multiple choices for automatic selection.
#' @param sampling (optional) sampling frequency of the dataset. If NULL, will 
#' try to identify from the timestamps of y. This is useful for plotting and 
#' extending the timestamps in the prediction horizon.
#' @param init_garch GARCH variance initialization method with options to use
#' the \dQuote{unconditional} variance or a \dQuote{sample} of length \dQuote{sample_n}.
#' @param sample_n the sample length to use if choosing to initialize the GARCH variance
#' using the \dQuote{sample} method.
#' @param top_n how many models to return from the top when using automatic model selection. 
#' If this is equal to 1 then the best selected model based on lowest AIC will be returned 
#' else a list of the top_n estimated models (based on AIC).
#' @param ... not used.
#' @details The specification performs some sanity checks on the arguments provided 
#' and sets up the required state space matrices and parameters which are used in 
#' the estimation stage.
#' @returns An object of class \dQuote{tsissm.spec} or \dQuote{tsissm.autospec}.
#' @references De Livera, Alysha M and Hyndman, Rob J and Snyder, Ralph D, 2011, 
#' Forecasting time series with complex seasonal patterns using exponential smoothing, 
#' \emph{Journal of the American Statistical Association}, \bold{106(496)}, 1513--1527.
#' @aliases issm_modelspec
#' @rdname issm_modelspec
#' @export
#'
#'
#'
#'
issm_modelspec <- function(y, auto = FALSE, slope = TRUE, slope_damped = FALSE, seasonal = FALSE, 
                           seasonal_frequency = 1, seasonal_harmonics = NULL, ar = 0, ma = 0, xreg = NULL,  
                           variance = "constant", garch_order = c(1,1), lambda = NULL, lower = 0, upper = 1, 
                           distribution = "norm", sampling = NULL, init_garch = "unconditional", 
                           sample_n = 10, top_n = 1, ...)
{
    if (!is.xts(y)) {
        stop("\ny is not an xts object...refusing to continue. Please fix and resubmit.")
    }
    if (any(y <= 0, na.rm = TRUE)) 
        stop("\nThe issm model only supports strictly positive data.")
    
    if (auto) {
        spec <- .auto_modelspec(y, slope = slope, slope_damped = slope_damped, seasonal = seasonal,
                               seasonal_frequency = seasonal_frequency, seasonal_harmonics = seasonal_harmonics,
                               ar = ar, ma = ma, xreg = xreg, variance = variance, garch_order = garch_order, 
                               lambda = lambda, lower = lower, upper = upper, distribution = distribution,
                               sampling = sampling, init_garch = init_garch, sample_n = sample_n, top_n = top_n)
    } else {
        spec <- .single_modelspec(y, slope = slope, slope_damped = slope_damped, seasonal = seasonal,
                                seasonal_frequency = seasonal_frequency, seasonal_harmonics = seasonal_harmonics,
                                ar = ar, ma = ma, xreg = xreg, variance = variance, garch_order = garch_order, 
                                lambda = lambda, lower = lower, upper = upper, distribution = distribution,
                                sampling = sampling, init_garch = init_garch, sample_n = sample_n) 
    }
    return(spec)
}



.single_modelspec <- function(y, slope = TRUE, slope_damped = FALSE, seasonal = FALSE, 
                           seasonal_frequency = 1, seasonal_harmonics = NULL, ar = 0, ma = 0, xreg = NULL,  
                           variance = "constant", garch_order = c(1,1), lambda = NULL, lower = 0, upper = 1, 
                           distribution = "norm", sampling = NULL, init_garch = "unconditional", 
                           sample_n = 10, ...)
{
    initial <- parameters <- NULL
    n <- NROW(y)
    if (is.null(sampling)) {
        sampling <- sampling_frequency(index(y))
    }
    good <- rep(1, NROW(y))
    if (any(is.na(y))) {
        good[which(is.na(y))] <- 0
    }
    if (!is.null(xreg)) {
        if (nrow(xreg) != n) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (n/2)) warning("\nnumber of regressors greater than half the length of y!")
        if (is.null(colnames(xreg))) 
            colnames(xreg) <- paste0("reg.", 1:ncol(xreg))
        include_xreg <- TRUE
        xreg <- coredata(xreg)
    } else {
        xreg <- matrix(0, ncol = 1, nrow = n)
        colnames(xreg) <- paste0("reg.", 1)
        include_xreg <- FALSE
    }
    
    spec <- list()
    spec$target$y_orig <- as.numeric(y)
    spec$target$index <- index(y)
    spec$target$sampling <- sampling
    spec$slope$include_slope <- as.logical(slope)
    if (slope) {
        if (slope_damped) {
            include_damped <- TRUE
        } else {
            include_damped <- FALSE
        }
        spec$slope$include_damped <- include_damped
    } else {
        spec$slope$include_damped <- FALSE
    }
    spec$seasonal$include_seasonal <- as.logical(seasonal)
    if (seasonal) {
        if (any(is.null(seasonal_frequency))) { 
            stop("\nseasonal_frequency cannot be NULL if seasonal is set to TRUE")
        }
        if (any(seasonal_frequency < 2)) {
            stop("\nseasonal_frequency must be 2 or greater")
        }
        seasonal_type <- "trigonometric"
        if (any(seasonal_frequency > (n/2))) {
            warning("\nseasonal_frequency contains a value(s) which are greater than half the length of y!")
        }
        if (length(seasonal_frequency) > 2) {
            seasonal_frequency <- sort(seasonal_frequency)
            if (any(diff(seasonal_frequency) == 0)) {
                stop("\nfound duplicate values in seasonal_frequency.")
            }
        }
    } else {
        seasonal_harmonics <- NULL
        seasonal_frequency <- 1
    }
    spec$seasonal$seasonal_harmonics <- seasonal_harmonics[1:length(seasonal_frequency)]
    spec$seasonal$seasonal_frequency <- seasonal_frequency
    spec$seasonal$seasonal_type <- "trigonometric"
    spec$xreg$include_xreg <- include_xreg
    spec$xreg$xreg <- xreg
    transformation <- "box-cox"
    if (is.null(lambda)) lambda <- 1
    if (is.na(lambda)) {
        include_lambda <- TRUE
        transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
        tmp <- transform$transform(y = y, frequency = seasonal_frequency[1])
        transform$lambda <- attr(tmp, "lambda")
        transform$include_lambda <- TRUE
        lambda <- transform$lambda
        transform$lower <- lower
        transform$upper <- upper
    } else {
        include_lambda <- FALSE
        transform <- tstransform(transformation, lambda = lambda, lower = lower, upper = upper)
        transform$include_lambda <- FALSE
        transform$lambda <- lambda
        transform$lower <- lower
        transform$upper <- upper
    }
    transform$name <- "box-cox"
    spec$target$y <- as.numeric(y)
    spec$transform <- transform
    spec$arma$order <- c(ar, ma)
    spec$variance$type <- match.arg(variance[1], c("constant","dynamic"))
    spec$variance$order <- c(pmax(0, garch_order[1]), pmax(0, garch_order[2]))
    if (spec$variance$type == "constant") {
        spec$variance$order <- c(0,0)
    }
    spec$variance$init_variance <- match.arg(init_garch[1], c("unconditional","sample"))
    spec$variance$sample_n <- as.integer(pmax(1, sample_n))
    # vmodel: [maxpq, arch, garch, sample_n]
    vmodel <- c(max(garch_order), garch_order[1], garch_order[2], sample_n)
    spec$variance$vmodel <- vmodel
    
    spec$distribution$type <- match.arg(distribution[1], c("norm","std","jsu"))
    if (distribution == "norm") {
        spec$distribution$shape <- FALSE
        spec$distribution$skew <- FALSE
        dist_flag <- 1
    } else if (distribution == "std") {
        distribution_table <- rbind(
            data.table(parameters = "skew", initial = 0, lower = 0, upper = 0, estimate = 0, scale = 1),
            data.table(parameters = "shape", initial = 5, lower = 4.1, upper = 50, estimate = 1, scale = 1))
        spec$distribution$shape <- TRUE
        spec$distribution$skew <- FALSE
        dist_flag <- 2
    } else if (distribution == "jsu") {
        distribution_table <- rbind(
            data.table(parameters = "shape", initial = 2, lower = 0.1, upper = 10, estimate = 1, scale = 1),
            data.table(parameters = "skew", initial = 0, lower = -20, upper = 20, estimate = 1, scale = 1))
        spec$distribution$shape <- TRUE
        spec$distribution$skew <- TRUE
        dist_flag <- 7
    } else {
        spec$distribution$shape <- FALSE
        spec$distribution$skew <- FALSE
        distribution_table <- rbind(
            data.table(parameters = "skew", initial = 0, lower = 0, upper = 0, estimate = 0, scale = 1),
            data.table(parameters = "shape", initial = 0, lower = 0, upper = 0, estimate = 0, scale = 1))
        dist_flag <- 1
    }
    spec$distribution$number <- dist_flag
    
    M <- ss_matrices(y = y, slope = spec$slope$include_slope, 
                     damped_slope = spec$slope$include_damped, frequency = spec$seasonal$seasonal_frequency, 
                     type = spec$seasonal$seasonal_type, K = spec$seasonal$seasonal_harmonics, 
                     ar = spec$arma$order[1], ma = spec$arma$order[2], include_xreg = include_xreg, 
                     xreg = xreg)
    S <- M$setup
    S <- rbind(S, data.table(matrix = "xreg", values = as.vector(coredata(xreg)), pars = NA))
    if (include_xreg) {
        S <- rbind(S, data.table(matrix = "kappa", values = rep(as.numeric(NA), ncol(xreg)), pars = paste0("kappa.",1:ncol(xreg))))
    } else {
        S <- rbind(S, data.table(matrix = "kappa", values = rep(as.numeric(0), ncol(xreg)), pars = paste0("kappa")))
    }
    # dims [states timesteps seasonal_type seasonal_start seasonal_end arma_terms xreg distribution_number]
    dims <- c(M$dims, "distribution" = dist_flag)
    setkeyv(S, "matrix")
    spec$S <- S
    spec$dims <- dims
    pars <- M$pars
    if (spec$variance$type == "constant" && spec$distribution$type == "norm" && max(spec$arma$order) == 0) {
        # we don't need res
        res <- NULL
    } else {
        res <- init_res(y, spec)
    }
    parmatrix <- issm_init_pars(pars, spec)
    garch_parmatrix <- garch_init_pars(res, spec)
    distribution_parmatrix <- distribution_init_pars(res, spec)
    parmatrix <- rbind(parmatrix, garch_parmatrix, distribution_parmatrix)
    if (max(spec$arma$order) > 0) {
        init_arma <- try(arima(res, order = c(spec$arma$order[1], 0, spec$arma$order[2]), method = "ML", transform.pars = TRUE), silent = TRUE)
        if (!inherits(init_arma,'try-error')) {
            arma_coef <- coef(init_arma)
            arma_coef <- arma_coef[-length(arma_coef)]
            parmatrix[estimate == 1 & grepl("^theta", parameters) | grepl("^psi", parameters), initial := unname(as.numeric(arma_coef))]
            # check bounds since arma parameter bounds are based on -1 and 1 but this is not a requirement for stationarity
            tmp <- parmatrix[estimate == 1 & (grepl("^theta", parameters) | grepl("^psi", parameters)) & (initial >= upper)]
            if (NROW(tmp) > 0) {
                parmatrix[estimate == 1 & (grepl("^theta", parameters) | grepl("^psi", parameters)) & (initial >= upper), initial := sign(initial) * 0.5]
            }
            # Check and update if 'initial' is less than or equal to 'lower'
            tmp <- parmatrix[estimate == 1 & (grepl("^theta", parameters) | grepl("^psi", parameters)) & (initial <= lower)]
            if (NROW(tmp) > 0) {
                parmatrix[estimate == 1 & (grepl("^theta", parameters) | grepl("^psi", parameters)) & (initial <= lower), initial := sign(initial) * 0.5]
            }
        }
    }
    spec$parmatrix <- parmatrix
    spec$idmatrix <- M$idmatrix
    spec$good <- good
    spec$valid_index <- which(good != 0)
    class(spec) <- c("tsissm.spec", "tsmodel.spec")
    return(spec)
}


.auto_modelspec <- function(y, slope = TRUE, slope_damped = FALSE, seasonal = FALSE, 
                            seasonal_frequency = 1, seasonal_harmonics = NULL, ar = 0, ma = 0, xreg = NULL,  
                            variance = "constant", garch_order = c(1,1), lambda = NULL, lower = 0, upper = 1, 
                            distribution = "norm", sampling = NULL, init_garch = "unconditional", 
                            sample_n = 10, top_n = 1)
{
    # logic for seasonal frequency is complex
    # TRUE: all(seasonal_frequency > 1)
    # FALSE: seasonal_frequency and seasonal_harmonics ignored
    # (FALSE, TRUE): check for TRUE case, do not include anything for FALSE case
    
    n <- NROW(y)
    if (is.null(sampling)) {
        sampling <- sampling_frequency(index(y))
    }
    if (!is.null(xreg)) {
        if (nrow(xreg) != n) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (n/2)) warning("\nnumber of regressors greater than half the length of y!")
    }
    
    valid_logical <- c(TRUE, FALSE)
    if (!all(slope %in% valid_logical)) stop("slope must be TRUE, FALSE, or both.")
    if (!all(slope_damped %in% valid_logical)) stop("slope_damped must be TRUE, FALSE, or both.")
    if (!all(seasonal %in% valid_logical)) stop("seasonal must be TRUE, FALSE, or both.")
    if (seasonal_frequency[1] == 1 | all(seasonal == FALSE)) {
        seasonal <- FALSE
        seasonal_harmonics = list()
        seasonal_frequency <- 1
    }
    # Validate seasonal_frequency
    if (!is.numeric(seasonal_frequency) || any(seasonal_frequency <= 0)) 
        stop("seasonal_frequency must be a positive.")
    if (length(seasonal_frequency) > 1 && (any(seasonal_frequency <= 1) || anyDuplicated(seasonal_frequency))) 
        stop("If seasonal_frequency has length > 1, all values must be > 1 and unique.")
    
    # Validate seasonal_harmonics
    if (any(seasonal) & seasonal_frequency[1] > 1 & !is.list(seasonal_harmonics)) stop("seasonal_harmonics must be a list.")
    if (any(seasonal_frequency > 1) && length(seasonal_harmonics) != length(seasonal_frequency)) 
        stop("seasonal_harmonics must have length equal to seasonal_frequency.")
    if (length(seasonal_frequency) == 1 && seasonal_frequency == 1 && length(seasonal_harmonics) > 0) {
        warning("seasonal_harmonics is ignored when seasonal_frequency = 1.")
        seasonal_harmonics <- list()
    }
    if (any(seasonal_frequency > 1)) {
        for (i in seq_along(seasonal_frequency)) {
            if (length(seasonal_harmonics[[i]]) == 0) {
                stop("seasonal_harmonics cannot be empty for frequency > 1.")
            } else {
                if (max(seasonal_harmonics[[i]]) > 0.5 * seasonal_frequency[i]) {
                    stop("Seasonal_harmonics cannot be > 0.5 * seasonal_frequency.")
                }
            }
        }
    }
    # Validate ar and ma
    if (!is.numeric(ar) || any(ar < 0)) stop("ar must be a numeric vector with non-negative values.")
    if (!is.numeric(ma) || any(ma < 0)) stop("ma must be a numeric vector with non-negative values.")
    
    # Validate lambda
    if (!is.null(lambda) && !is.na(lambda) && (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0 || lambda > 2)) 
        stop("lambda must be NULL, NA, or a numeric value between 0 and 2.")
    
    # Validate variance
    valid_variance <- c("constant", "dynamic")
    if (!all(variance %in% valid_variance)) 
        stop('variance must be either "constant", "dynamic", or both.')
    if (any(variance == "dynamic")) {
        if (length(garch_order) != 2) stop("\ngarch_order must be a vector of length 2 (no selection for GARCH)")
        garch_order <- pmax(garch_order, 0, na.rm = TRUE)
        init_garch <- match.arg(init_garch[1], c("unconditional","sample"))
        sample_n <- max(sample_n[1], 1)
    } else {
        garch_order <- c(1,1)
        init_garch <- "unconditional"
        sample_n <- 10
    }
    
    # Validate lower and upper
    if (!is.numeric(lower) || lower < 0) stop("lower must be positive.")
    if (!is.numeric(upper) || upper <= lower) stop("upper must be greater than lower.")
    
    # Validate distribution
    valid_distribution <- c("norm", "std", "jsu")
    if (!distribution %in% valid_distribution) 
        stop('distribution must be one of "norm", "std", or "jsu".')
    
    # Validate top_n
    if (!is.numeric(top_n) || length(top_n) != 1 || top_n <= 0 || top_n > 20) 
        stop("top_n must be a positive number <= 20.")
    
    parameters <- NULL
    args_grid_nonseasonal <- NULL
    args_grid_seasonal <- NULL
    if (any(seasonal == FALSE)) {
        args_list <- list()
        args_list$slope <- unique(as.logical(slope))
        args_list$slope_damped <- unique(as.logical(slope_damped))
        args_list$seasonal <- FALSE
        args_list$ar <- ar
        args_list$ma <- ma
        args_list_non_seasonal <- args_list
        args_grid_nonseasonal <- as.data.frame(expand.grid(args_list))
    }
    if (any(seasonal == TRUE)) {
        args_list <- list()
        args_list$slope <- unique(as.logical(slope))
        args_list$slope_damped <- unique(as.logical(slope_damped))
        args_list$seasonal <- TRUE
        args_list$ar <- ar
        args_list$ma <- ma
        s_list <- vector(mode = "list", length = length(seasonal_frequency))
        for (i in 1:length(seasonal_frequency)) {
            s_list[[i]] <- seasonal_harmonics[[i]]
        }
        names(s_list) <- paste0("Seasonal", seasonal_frequency)
        args_list <- c(args_list, s_list)
        args_grid_seasonal <- as.data.frame(expand.grid(args_list))
        if (!is.null(args_grid_nonseasonal)) {
            tmp <- matrix(0, ncol = length(seasonal_frequency), nrow = nrow(args_grid_nonseasonal))
            colnames(tmp) <- paste0("Seasonal", seasonal_frequency)
            args_grid_nonseasonal <- cbind(args_grid_nonseasonal, as.data.frame(tmp))
        }
    }
    args_grid <- rbind(args_grid_nonseasonal, args_grid_seasonal)
    if (any(args_grid$slope)) {
        if (any(args_grid$slope_damped)) {
            exc <- which(args_grid$slope == FALSE & args_grid$slope_damped)
            if (length(exc) > 0) args_grid <- args_grid[-exc,]
        }
    }
    variance <- check_variance(variance)
    
    if (length(variance) == 2) {
        vtype <- data.frame(variance = rep("constant", NROW(args_grid)))
        tmp1 <- cbind(args_grid, vtype)
        vtype <- data.frame(variance = rep("dynamic", NROW(args_grid)))
        tmp2 <- cbind(args_grid, vtype)
        args_grid <- rbind(tmp1, tmp2)
    } else {
        vtype <- data.frame(variance = rep(variance[1], NROW(args_grid)))
        args_grid <- cbind(args_grid, vtype)
    }
    n <- NROW(args_grid)
    
    if (top_n > n) {
        warning("\ntop_n is greater than no. of models to evaluate. Will return all models instead.")
        top_n <- n
    }
    args_grid <- cbind(args_grid, data.frame(distribution = rep(distribution, nrow(args_grid))))
    gnames <- colnames(args_grid)
    spec <- list(y = y, slope = slope, slope_damped = slope_damped, seasonal = seasonal, seasonal_frequency = seasonal_frequency, 
                 seasonal_harmonics = seasonal_harmonics, ar = ar, ma = ma, lambda = lambda, lower = lower,  
                 upper = upper, variance = variance, lower = lower, upper = upper, distribution = distribution, 
                 sampling = sampling, xreg = xreg, top_n = top_n, garch_order = garch_order, init_garch = init_garch, sample_n = sample_n,
                 grid = args_grid)
    class(spec) <- c("tsissm.autospec","tsmodel.spec") 
    return(spec)
}


issm_init_pars <- function(pars, spec)
{
    tmp <- data.table(parameters = "alpha", initial = 0.11, lower = 1e-12, upper = 0.99, estimate = 1, scale = 1, group = "issmpars", equation = "[M]", symbol = "\\alpha")
    if (spec$slope$include_slope) {
        tmp <- rbind(tmp, data.table(parameters = "beta", initial = 0.01, lower = 1e-12, upper = 0.9, estimate = 1, scale = 1, group = "issmpars", equation = "[M]", symbol = "\\beta"))
        if (spec$slope$include_damped) {
            tmp <- rbind(tmp, data.table(parameters = "phi", initial = 0.98, lower = 0.5, upper = 1, estimate = 1, scale = 1, group = "issmpars", equation = "[M]", symbol = "\\phi"))
        }
    }
    if (spec$seasonal$include_seasonal) {
        snames <- pars[grepl("gamma", pars)]
        sequation <- sapply(snames, seasonal_symbol_maker)
        stmp <- data.table(parameters = snames, initial = rep(1e-8, length(snames)), lower = -5e-2, upper = 1, estimate = 1, scale = 1, group = "issmpars", equation = "[M]", symbol = sequation)
        tmp <- rbind(tmp, stmp)
    }
    
    
    if (spec$arma$order[1] > 0) {
        tmp <- rbind(tmp, data.table(parameters = paste0("theta",1:spec$arma$order[1]), initial = 1e-03, lower = -0.99, upper = 0.99, estimate = 1, scale = 1, group = "issmpars", equation = "[M]", symbol = paste0("\\theta_{",1:spec$arma$order[1],"}")))
    }
    
    if (spec$arma$order[2] > 0) {
        tmp <- rbind(tmp, data.table(parameters = paste0("psi",1:spec$arma$order[2]), initial = 1e-03, lower = -0.99, upper = 0.99, estimate = 1, scale = 1, group = "issmpars", equation = "[M]", symbol = paste0("\\psi_{",1:spec$arma$order[2],"}")))
    }
    
    if (spec$transform$include_lambda) {
        # pre-estrimate lambda (initial value)
        lam <- lambda_init(spec$target$y, frequency = spec$seasonal$seasonal_frequency[1])
        tmp <- rbind(tmp, data.table(parameters = "lambda", initial = lam, lower = 1e-12, upper = 1.5, estimate = 1, scale = 1, group = "lambda", equation = "[M]", symbol = "\\lambda"))
    } else {
        tmp <- rbind(tmp, data.table(parameters = "lambda", initial = spec$transform$lambda, lower = 1e-12, upper = 1.5, estimate = 0, scale = 1, group = "lambda", equation = "[M]", symbol = "\\lambda"))
    }
    lower_x <- -1e8
    upper_x <-  1e8
    if (spec$xreg$include_xreg) {
        tmp <- rbind(tmp, data.table(parameters = paste0("kappa.",1:ncol(spec$xreg$xreg)), initial = 1, lower = lower_x, upper = upper_x, estimate = 1, scale = 1, group = "kappa", equation = "[R]", symbol = paste0("\\kappa_{",1:ncol(spec$xreg$xreg),"}")))
    } else {
        tmp <- rbind(tmp, data.table(parameters = "kappa", initial = 0, lower = lower_x, upper = upper_x, estimate = 0, scale = 1, group = "kappa", equation = "[R]", symbol = "\\kappa"))
    }
    return(tmp)
}

init_res <- function(y, spec) {
    lambda <- spec$transform$lambda
    seasonal_frequency <- spec$seasonal$seasonal_frequency[1]
    if (spec$variance$type == "constant" | length(y) > 10000L) {
        seasonal <- ifelse(seasonal_frequency > 1, TRUE, FALSE)
        if (is.null(lambda) | lambda == 1.0) {
            yt <- y
            res <- residuals(tslinear(yt, trend = TRUE, seasonal = seasonal, frequency = seasonal_frequency))
        } else if (is.na(lambda)) {
            B <- box_cox(lambda = NA)
            yt <- B$transform(y, frequency = seasonal_frequency)
            res <- residuals(tslinear(yt, trend = TRUE, seasonal = seasonal, frequency = seasonal_frequency))
        } else {
            B <- box_cox(lambda = lambda)
            yt <- B$transform(y, frequency = seasonal_frequency)
            res <- residuals(tslinear(yt, trend = TRUE, seasonal = seasonal, frequency = seasonal_frequency))
        }
    } else {
        if (is.null(lambda) | lambda == 1.0) {
            yt <- y
            res <- try(residuals(arima(yt, order = c(1,1,2), seasonal = list(order = c(1, 1, 1), period = seasonal_frequency))), silent = TRUE)
        } else if (is.na(lambda)) {
            B <- box_cox(lambda = NA)
            yt <- B$transform(y, frequency = seasonal_frequency)
            res <- try(residuals(arima(yt, order = c(1,1,2), seasonal = list(order = c(1, 1, 1), period = seasonal_frequency))), silent = TRUE)
        } else {
            B <- box_cox(lambda = lambda)
            yt <- B$transform(y, frequency = seasonal_frequency)
            if (seasonal_frequency > 1) slist <- c(0,0,0) else slist <- c(1,1,1)
            res <- try(residuals(arima(yt, order = c(1,1,2), seasonal = list(order = slist, period = seasonal_frequency))), silent = TRUE)
        }
        if (inherits(res,'try-error')) {
            res <- residuals(tslinear(yt, trend = TRUE, seasonal = ifelse(seasonal_frequency>1, TRUE, FALSE), frequency = seasonal_frequency))
        }
    }
    res <- as.numeric(res)
    return(res)
}

garch_init_pars <- function(res, spec)
{
    garch_order <- spec$variance$order
    if (sum(garch_order) > 0) {
        alpha <- 0.05/garch_order[1]
        beta <- 0.8/garch_order[2]
        omega <- mean(res^2, na.rm = TRUE) * (1 - sum(alpha) - sum(beta))
    } else {
        alpha <- beta <- omega <- 0
    }
    if (spec$variance$type == "dynamic") {
        garch_table <- NULL
        if (garch_order[1] > 0) {
            if (garch_order[1] == 1) {
                tmp <- data.table(parameters = "eta", initial = 0.01, lower = 0, upper = 0.99, estimate = 1, scale = 1, group = "eta", equation = "[V]", symbol = "\\eta")
            } else {
                tmp <- data.table(parameters = paste0("eta",1:garch_order[1]), initial = 0.01, lower = 0, upper = 0.99, estimate = 1, scale = 1, group = "eta", equation = "[V]", symbol = paste0("\\eta_",1:garch_order[1]))
            } 
        } else {
            tmp <- data.table(parameters = "eta", initial = 0, lower = 0, upper = 0.99, estimate = 0, scale = 1, group = "eta", equation = "[V]", symbol = "\\eta")
        }
        garch_table <- rbind(garch_table, tmp)
        if (garch_order[2] > 0) {
            if (garch_order[2] == 1) {
                tmp <- data.table(parameters = "delta", initial = 0.8, lower = 0, upper = 0.99, estimate = 1, scale = 1, group = "delta", equation = "[V]", symbol = "\\delta")
            } else {
                tmp <- data.table(parameters = paste0("delta",1:garch_order[2]), initial = c(0.8)/garch_order[2], lower = 0, upper = 0.99, estimate = 1, scale = 1, group = "delta", equation = "[V]", symbol = paste0("\\delta_",1:garch_order[2]))
            } 
        } else {
            tmp <- data.table(parameters = "delta", initial = 0, lower = 0, upper = 0.99, estimate = 0, scale = 1, group = "delta", equation = "[V]", symbol = "\\delta")
        }
        garch_table <- rbind(garch_table, tmp)
    } else {
        garch_table <- rbind(
            tmp <- data.table(parameters = "eta", initial = 0, lower = 0, upper = 0.99, estimate = 0, scale = 1, group = "eta", equation = "[V]", symbol = "\\eta"),
            tmp <- data.table(parameters = "delta", initial = 0, lower = 0, upper = 0.99, estimate = 0, scale = 1, group = "delta", equation = "[V]", symbol = "\\delta")
        )
    }
    return(garch_table)
}

distribution_init_pars <- function(res, spec) {
    distribution <- spec$distribution$type
    if (distribution != "norm") {
        spec <- distribution_modelspec(na.omit(res), distribution = distribution)
        mod <- try(estimate(spec, use_hessian = FALSE), silent = TRUE)
        if (inherits(mod,'try-error')) {
            cf <- c("skew" = 1, "shape" = ifelse(distribution == "std", 4.1, 2.1))
        } else {
            cf <- coef(mod)
        }
        if (distribution == "std") {
            skew <- 0
            shape <- cf["shape"]
            tmp <- rbind(
                data.table(parameters = "skew", initial = skew, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                data.table(parameters = "shape", initial = shape, lower = 2.01, upper = 100, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"))
        } else if (distribution == "jsu") {
            skew <- cf["skew"]
            shape <- cf["shape"]
            tmp <- rbind(
                data.table(parameters = "skew", initial = skew, lower = -20, upper = 20, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
                data.table(parameters = "shape", initial = shape, lower = 0.1, upper = 100, estimate = 1, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"))
        }
    } else {
        skew <- shape <- 0.0
        tmp <- rbind(
            data.table(parameters = "skew", initial = skew, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\zeta"),
            data.table(parameters = "shape", initial = shape, lower = 0, upper = 0, estimate = 0, scale = 1, group = "distribution", equation = "[D]", symbol = "\\nu"))
    }
    return(tmp)
}


lambda_init <- function(y, frequency = 1)
{
    f <- tstransform(lambda = NA)
    x <- f$transform(y, frequency = frequency)
    lambda <- attr(x,"lambda")
    return(lambda)
}

#' Model Specification Extractor
#'
#' @description Extracts a model specification (class \dQuote{tsissm.spec}) from
#' an object of class \dQuote{tsissm.estimate}.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param y an optional new xts vector.
#' @param lambda an optional lambda parameter for the Box Cox transformation (if
#' previously used).
#' @param xreg an optional matrix of regressors.
#' @param parameters whether to set the parameters to their initial values, \dQuote{initial}
#' or the optimally estimated \dQuote{optimal}.
#' @param ... not currently used.
#' @note This function is used by other functions in the package such as the
#' backtest which requires rebuilding the specification for each re-estimation
#' step with updated data but keeping all else equal.
#' @return An object of class \dQuote{tsissm.spec}.
#' @aliases tsspec
#' @method tsspec tsissm.estimate
#' @rdname tsspec
#' @export
#'
#'
tsspec.tsissm.estimate <- function(object, y = NULL, lambda = NULL, xreg = NULL, parameters = "initial", ...)
{
    if (is.null(y)) {
        y <- xts(object$spec$target$y_orig, object$spec$target$index)
    }
    if (!is.null(xreg)) {
        xreg <- coredata(xreg)
        if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
        if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
    } else {
        if (object$spec$xreg$include_xreg) {
            xreg <- object$spec$xreg$xreg
            if (nrow(xreg) != NROW(y)) stop("\nxreg should have the same number of rows as y")
            if (ncol(xreg) > (NROW(y)/2)) warning("\nnumber of regressors greater than half the length of y!")
        } else {
            xreg <- NULL
        }
    }
    if (is.null(lambda)) {
        # check whether lambda was originally estimated
        if (object$spec$parmatrix[parameters == "lambda"]$estimate == 1) {
            lambda <- NA
        } else {
            lambda <- object$parmatrix[parameters == "lambda"]$value
        }
    }
    spec <- issm_modelspec(y, slope = object$spec$slope$include_slope,
                   slope_damped = object$spec$slope$include_damped,
                   seasonal = object$spec$seasonal$include_seasonal,
                   seasonal_frequency = object$spec$seasonal$seasonal_frequency,
                   lower = object$spec$transform$lower,
                   upper = object$spec$transform$upper,
                   seasonal_harmonics = object$spec$seasonal$seasonal_harmonics,
                   ar = object$spec$arma$order[1],
                   ma = object$spec$arma$order[2],
                   xreg = xreg, lambda = lambda, 
                   variance = object$spec$variance$type, 
                   garch_order = object$spec$variance$order, 
                   init_garch = object$spec$variance$init_garch, 
                   sample_n = object$spec$variance$sample_n,
                   distribution = object$spec$distribution$type,
                   sampling = object$spec$target$sampling)
    if (parameters != "initial") {
        spec$parmatrix$initial <- object$parmatrix$value
    }
    return(spec)
}

get_pars.tsissm.spec <- function(object)
{
    return(object$parmatrix)
}

check_variance <- function(variance) {
    allowed_values <- c("constant", "dynamic")
    if (!is.null(variance)) {
        if (!is.character(variance) || !all(variance %in% allowed_values)) {
            stop("'variance' must be NULL or a character vector containing 'constant' and/or 'dynamic'.")
        }
    } else {
        variance <- "constant"
    }
    return(variance)
}
