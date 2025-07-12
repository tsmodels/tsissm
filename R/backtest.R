#' Walk Forward Model Backtest
#'
#' @description Generates an expanding window walk forward backtest.
#' @param object an object of class \dQuote{tsissm.spec} \dQuote{tsissm.autospec}.
#' @param start numeric data index from which to start the backtest.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param estimate_every number of periods at which the model is re-estimated
#' and new predictions are generated (defaults to 1).
#' @param rolling this indicates whether forecasts are made only on the estimation
#' date (FALSE) or whether to filter the data 1 period at a time and forecast
#' from the filtered data (TRUE).
#' @param weights_scheme the weighting scheme to use when using ensembling (see note).
#' @param weights a vector of fixed user supplied weights of length \code{top_n} when choosing
#' \dQuote{U} in the \dQuote{weights_scheme}.
#' @param seed an value specifying if and how the random number generator should 
#' be initialized (\sQuote{seeded}). Either NULL or an integer that will be used in a 
#' call to set.seed before simulating the response vectors.
#' @param solver only \dQuote{nlortr} currently supported.
#' @param control optional control parameters.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... not used.
#' @return A list with the following data.tables:
#' \itemize{
#' \item prediction : the backtest table with forecasts and actuals
#' \item metrics: a summary performance table showing metrics by
#' forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
#' }
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package. Model ensembling
#' is used when the input object is of class \dQuote{tsissm.autospec} and
#' top_n is greeter than 1. The following weighting schemes are available:
#' \describe{
#'   \item{U:}{User supplied fixed weights.}
#'   \item{AIC:}{Models are weighted based on their Akaike Information Criterion (AIC), 
#'         favoring models with lower AIC values. The weights are computed as:}
#'   \deqn{w_i = \frac{\exp(-0.5 \times \Delta_i)}{\sum_{j} \exp(-0.5 \times \Delta_j)}}
#'   where \deqn{\Delta_i = AIC_i - \min(AIC)}
#'   \item{BIC:}{Similar to AIC-based weighting but uses the Bayesian Information Criterion (BIC) 
#'         instead, which penalizes model complexity more strongly.}
#' }
#' The final ensemble prediction is computed as:
#' \deqn{\hat{y}_{\text{ensemble}} = \sum_{i} w_i \hat{y}_i}
#' where \eqn{w_i} are the model weights and \eqn{\hat{y}_i} are the individual model predictions.
#' AIC-based weighting is preferred when models have different complexities but are estimated on the same dataset, 
#' while BIC-based weighting may be better suited for large sample sizes due to its stronger penalty on complexity.
#' @aliases tsbacktest
#' @seealso [tsensemble()]
#' @method tsbacktest tsissm.autospec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.tsissm.autospec <- function(object, start = floor(length(object$y)/2), end = length(object$y),
                                   h = 1, estimate_every = 1, rolling = FALSE, weights_scheme = c("AIC","BIC","U"), weights = NULL,
                                   seed = NULL, solver = "nloptr", control = NULL, trace = FALSE, ...)
{
    solver <- "nloptr"
    if (is.null(control)) {
        control <- issm_control(solver = "nloptr", algorithm = "SLSQP", trace = 0)
    }
    if (object$top_n > 1) {
        top_n <- object$top_n
        weights_scheme <- match.arg(weights_scheme[1], c("AIC","BIC","U"))
        if (!is.null(weights)) {
            if (length(as.numeric(weights)) != top_n) stop(paste0("weights must be a vector of length (top_n): ", top_n))
            if (sum(weights) != 1) warning("\nweights do not sum to 1.")
            weights_scheme <- "U"
        }
        out <- .backtest_ensemble(object, start = start, end = end, h = h, estimate_every = estimate_every, rolling = rolling, weights_scheme = weights_scheme, 
                                  weights = weights, seed = seed, solver = solver, control = control, trace = trace)
    } else {
        out <- .backtest_top(object, start = start, end = end, h = h, estimate_every = estimate_every, rolling = rolling, seed = seed, solver = solver, 
                             control = control, trace = trace)
    }
    return(out)
}

#' @method tsbacktest tsissm.spec
#' @rdname tsbacktest
#' @export
tsbacktest.tsissm.spec <- function(object, start = floor(length(object$target$y_orig)/2), end = length(object$target$y_orig),
                                   h = 1, estimate_every = 1, rolling = FALSE, seed = NULL, solver = "nloptr", control = NULL, trace = FALSE, ...)
{
    solver <- "nloptr"
    if (is.null(control)) {
        control <- issm_control(solver = "nloptr", algorithm = "SLSQP", trace = 0)
    }
    parameter <- b <- forecast_dates <- NULL
    data <- xts(object$target$y_orig, object$target$index)
    
    if (start > NROW(data)) stop("\nstart is greater than length of data")
    if (end > NROW(data)) stop("\nend is greater than length of data")
    
    
    if (object$xreg$include_xreg) {
        use_xreg <- TRUE
        xreg <- xts(object$xreg$xreg, object$target$index)
    } else {
        use_xreg <- FALSE
        xreg <- NULL
    }
    start_date <- index(data)[start]
    end <- min(NROW(data), end)
    end_date <- index(data)[end - 1]
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    if (estimate_every != 1) {
        estimate_every <- max(1, as.integer(estimate_every))
        ns <- length(seqdates)
        seqdates <- seqdates[seq(1, ns, by = estimate_every)]
    }
    idates <- index(data)
    # check for ending period to avoid overlap
    index_table <- rbindlist(lapply(1:length(seqdates), function(i) {
        if (i == length(seqdates)) {
            s <- index(data[paste0(seqdates[i],"/",max(idates[end]))])
            s <- s[-length(s)]
        } else {
            s <- index(data[paste0(seqdates[i],"/",seqdates[i + 1])])
            s <- s[-length(s)]
        }
        if (length(s) == 0) return(NULL)
        rbindlist(lapply(1:length(s), function(j){
            idx <- which(idates == s[j])
            maxh <- min(idx + h, NROW(data[1:end]))
            tmp <- na.omit(index(data[idx:maxh]))
            data.table(estimate_date = seqdates[i], filter_date = tmp[1], forecast_dates = tmp[-1], h = length(tmp[-1]))
        }))
    }))
    max_index <- max(idates[end])
    index_table <- index_table[forecast_dates <= max_index]
    index_table <- split(index_table, by = "estimate_date")
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    i <- 1
    b <- future_lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        y_train <- data[paste0("/", seqdates[i])]
        if (use_xreg) {
            xreg_train <- xreg[index(y_train)]
        } else {
            xreg_train <- NULL
        }
        if (object$transform$include_lambda) {
            lambda <- NA
        } else {
            lambda <- object$transform$lambda
        }
        spec <- issm_modelspec(y_train, slope = object$slope$include_slope, slope_damped = object$slope$include_damped,
                               seasonal = object$seasonal$include_seasonal,
                               seasonal_frequency = object$seasonal$seasonal_frequency,
                               seasonal_harmonics = object$seasonal$seasonal_harmonics,
                               ar = object$arma$order[1], ma = object$arma$order[2], 
                               xreg = xreg_train, lambda = lambda, lower = object$transform$lower, upper = object$transform$upper, 
                               sampling = object$target$sampling, sample_n = object$variance$sample_n, 
                               init_garch =  object$variance$init_variance, garch_order = object$variance$order, 
                               variance = object$variance$type, distribution = object$distribution$type)
        mod <- try(estimate(spec, solver = solver, control = control, scores = FALSE), silent = TRUE)
        model_coef <- coef(mod)
        log_lik <- as.numeric(logLik(mod))
        aic <- as.numeric(AIC(mod))
        L <- index_table[[i]]
        M <- split(L, by = "filter_date")
        rolls <- length(M)
        # first index is the estimation so we start with the prediction
        # xreg
        if (use_xreg) {
            xreg_test <- xreg[M[[1]]$forecast_dates]
        } else {
            xreg_test <- NULL
        }
        # actual
        P <- vector(mode = "list", length = rolls)
        tmp <- predict(mod, h = M[[1]]$h[1], newvreg = xreg_test, forc_dates = M[[1]]$forecast_date, nsim = 2000, seed = seed)
        if (rolls > 1 & rolling) {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            P[[1]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                 "filter_date" = M[[1]]$filter_date,
                                 "horizon" =  1:M[[1]]$h[1],
                                 "size" = rep(nrow(y_train), nh),
                                 "forecast_date" =  M[[1]]$forecast_dates,
                                 "forecast" = as.numeric(tmp$analytic_mean),
                                 "sigma" = apply(tmp$distribution, 2, sd),
                                 "p10" = apply(tmp$distribution, 2, quantile, 0.1),
                                 "p90" = apply(tmp$distribution, 2, quantile, 0.9),
                                 "min" = apply(tmp$distribution, 2, min),
                                 "max" = apply(tmp$distribution, 2, max),
                                 "actual" = as.numeric(y_test))
            for (j in 2:rolls) {
                roll_dates <- M[[j]]$filter_date[1]
                updated_y <- data[paste0(seqdates[i],"/",roll_dates)][-1,]
                if (use_xreg) {
                    newxreg <- xreg[paste0(seqdates[i],"/",roll_dates)][-1,]
                } else {
                    newxreg <- NULL
                }
                f <- tsfilter(mod, y = updated_y, newxreg = newxreg)
                if (use_xreg) {
                    xreg_test <- xreg[M[[j]]$forecast_dates]
                } else {
                    xreg_test <- NULL
                }
                tmp <- predict(f, h = M[[j]]$h[1], newxreg = xreg_test, forc_dates = M[[j]]$forecast_dates, nsim = 2000)
                y_test <- data[M[[j]]$forecast_dates]
                nh <- length(M[[j]]$h)
                P[[j]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                     "filter_date" = M[[j]]$filter_date,
                                     "horizon" =  1:M[[j]]$h[1],
                                     "size" = rep(nrow(y_train), nh),
                                     "forecast_date" =  M[[j]]$forecast_dates,
                                     "forecast" = as.numeric(tmp$analytic_mean),
                                     "sigma" = apply(tmp$distribution, 2, sd),
                                     "p10" = apply(tmp$distribution, 2, quantile, 0.1),
                                     "p90" = apply(tmp$distribution, 2, quantile, 0.9),
                                     "min" = apply(tmp$distribution, 2, min),
                                     "max" = apply(tmp$distribution, 2, max),
                                     "actual" = as.numeric(y_test))
            }
            out <- rbindlist(P)
        } else {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            out <- data.table("estimation_date" = rep(seqdates[i], nh),
                              "filter_date" = M[[1]]$filter_date,
                              "horizon" =  1:M[[1]]$h[1],
                              "size" = rep(nrow(y_train), nh),
                              "forecast_date" =  M[[1]]$forecast_dates,
                              "forecast" = as.numeric(tmp$analytic_mean),
                              "sigma" = apply(tmp$distribution, 2, sd),
                              "p10" = apply(tmp$distribution, 2, quantile, 0.1),
                              "p90" = apply(tmp$distribution, 2, quantile, 0.9),
                              "min" = apply(tmp$distribution, 2, min),
                              "max" = apply(tmp$distribution, 2, max),
                              "actual" = as.numeric(y_test))
        }
        return(out)
    }, future.packages = c("tsmethods","tsissm","xts","data.table"), future.stdout	= FALSE, future.seed = FALSE)
    b <- eval(b)
    b <- rbindlist(b)
    out <- list(table = b, h = h, estimate_every = estimate_every, rolling = rolling, type = "simple")
    class(out) <- "tsissm.backtest"
    return(out)

}

.backtest_top <- function(object, start = floor(length(object$y)/2), end = length(object$y), h = 1,
                          estimate_every = 1, rolling = FALSE, seed = NULL, solver = "nloptr", control = NULL, 
                          trace = FALSE, ...)
{
    parameter <- b <- forecast_dates <- NULL
    data <- object$y
    if (!is.null(object$xreg)) {
        use_xreg <- TRUE
        xreg <- xts(coredata(object$xreg), index(data))
    } else {
        use_xreg <- FALSE
        xreg <- NULL
    }
    start_date <- index(data)[start]
    end <- min(NROW(data), end)
    end_date <- index(data)[end - 1]
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    if (estimate_every != 1) {
        estimate_every <- max(1, as.integer(estimate_every))
        ns <- length(seqdates)
        seqdates <- seqdates[seq(1, ns, by = estimate_every)]
    }
    idates <- index(data)
    # check for ending period to avoid overlap
    index_table <- rbindlist(lapply(1:length(seqdates), function(i) {
        if (i == length(seqdates)) {
            s <- index(data[paste0(seqdates[i],"/",max(idates[end]))])
            s <- s[-length(s)]
        } else {
            s <- index(data[paste0(seqdates[i],"/",seqdates[i + 1])])
            s <- s[-length(s)]
        }
        if (length(s) == 0) return(NULL)
        rbindlist(lapply(1:length(s), function(j){
            idx <- which(idates == s[j])
            maxh <- min(idx + h, NROW(data[1:end]))
            tmp <- na.omit(index(data[idx:maxh]))
            data.table(estimate_date = seqdates[i], filter_date = tmp[1], forecast_dates = tmp[-1], h = length(tmp[-1]))
        }))
    }))
    max_index <- max(idates[end])
    index_table <- index_table[forecast_dates <= max_index]
    index_table <- split(index_table, by = "estimate_date")
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    i <- 1
    b <- lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        y_train <- data[paste0("/", seqdates[i])]
        if (use_xreg) {
            xreg_train <- xreg[index(y_train)]
        } else {
            xreg_train <- NULL
        }
        spec <- issm_modelspec(y_train, auto = TRUE, slope = object$slope, slope_damped = object$slope_damped,
                               seasonal = object$seasonal,
                               seasonal_frequency = object$seasonal_frequency,
                               seasonal_harmonics = object$seasonal_harmonics,
                               ar = object$ar, ma = object$ma, 
                               xreg = xreg_train, lambda = object$lambda, lower = object$lower, upper = object$upper, 
                               sampling = object$sampling, sample_n = object$sample_n, 
                               init_garch =  object$init_garch, garch_order = object$garch_order, 
                               variance = object$variance, distribution = object$distribution, top_n = 1)
        mod <- try(estimate(spec, solver = solver, control = control, trace = FALSE), silent = TRUE)
        model_coef <- coef(mod)
        log_lik <- as.numeric(logLik(mod))
        aic <- as.numeric(AIC(mod))
        L <- index_table[[i]]
        M <- split(L, by = "filter_date")
        rolls <- length(M)
        # first index is the estimation so we start with the prediction
        # xreg
        if (use_xreg) {
            xreg_test <- xreg[M[[1]]$forecast_dates]
        } else {
            xreg_test <- NULL
        }
        # actual
        P <- vector(mode = "list", length = rolls)
        tmp <- predict(mod, h = M[[1]]$h[1], newvreg = xreg_test, forc_dates = M[[1]]$forecast_date, nsim = 2000, seed = seed)
        if (rolls > 1 & rolling) {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            P[[1]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                 "filter_date" = M[[1]]$filter_date,
                                 "horizon" =  1:M[[1]]$h[1],
                                 "size" = rep(nrow(y_train), nh),
                                 "forecast_date" =  M[[1]]$forecast_dates,
                                 "forecast" = as.numeric(tmp$analytic_mean),
                                 "sigma" = apply(tmp$distribution, 2, sd),
                                 "p10" = apply(tmp$distribution, 2, quantile, 0.1),
                                 "p90" = apply(tmp$distribution, 2, quantile, 0.9),
                                 "min" = apply(tmp$distribution, 2, min),
                                 "max" = apply(tmp$distribution, 2, max),
                                 "actual" = as.numeric(y_test))
            for (j in 2:rolls) {
                roll_dates <- M[[j]]$filter_date[1]
                updated_y <- data[paste0(seqdates[i],"/",roll_dates)][-1,]
                if (use_xreg) {
                    newxreg <- xreg[paste0(seqdates[i],"/",roll_dates)][-1,]
                } else {
                    newxreg <- NULL
                }
                f <- tsfilter(mod, y = updated_y, newxreg = newxreg)
                if (use_xreg) {
                    xreg_test <- xreg[M[[j]]$forecast_dates]
                } else {
                    xreg_test <- NULL
                }
                tmp <- predict(f, h = M[[j]]$h[1], newxreg = xreg_test, forc_dates = M[[j]]$forecast_dates, nsim = 2000)
                y_test <- data[M[[j]]$forecast_dates]
                nh <- length(M[[j]]$h)
                P[[j]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                     "filter_date" = M[[j]]$filter_date,
                                     "horizon" =  1:M[[j]]$h[1],
                                     "size" = rep(nrow(y_train), nh),
                                     "forecast_date" =  M[[j]]$forecast_dates,
                                     "forecast" = as.numeric(tmp$analytic_mean),
                                     "sigma" = apply(tmp$distribution, 2, sd),
                                     "p10" = apply(tmp$distribution, 2, quantile, 0.1),
                                     "p90" = apply(tmp$distribution, 2, quantile, 0.9),
                                     "min" = apply(tmp$distribution, 2, min),
                                     "max" = apply(tmp$distribution, 2, max),
                                     "actual" = as.numeric(y_test))
            }
            out <- rbindlist(P)
        } else {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            out <- data.table("estimation_date" = rep(seqdates[i], nh),
                              "filter_date" = M[[1]]$filter_date,
                              "horizon" =  1:M[[1]]$h[1],
                              "size" = rep(nrow(y_train), nh),
                              "forecast_date" =  M[[1]]$forecast_dates,
                              "forecast" = as.numeric(tmp$analytic_mean),
                              "sigma" = apply(tmp$distribution, 2, sd),
                              "p10" = apply(tmp$distribution, 2, quantile, 0.1),
                              "p90" = apply(tmp$distribution, 2, quantile, 0.9),
                              "min" = apply(tmp$distribution, 2, min),
                              "max" = apply(tmp$distribution, 2, max),
                              "actual" = as.numeric(y_test))
        }
        return(out)
    })
    b <- rbindlist(b)
    out <- list(table = b, h = h, estimate_every = estimate_every, rolling = rolling, type = "simple")
    class(out) <- "tsissm.backtest"
    return(out)
}

.backtest_ensemble <- function(object, start = floor(length(object$y)/2), end = length(object$y),
                               h = 1, estimate_every = 1, rolling = FALSE, weights_scheme = c("U","AIC","BIC"), weights = NULL,
                               seed = NULL, solver = "nloptr", control = NULL, trace = FALSE, ...)
{
    if (weights_scheme == "U") {
        wfun <- function(x) {
            fixed_weights(x, weights)
        }
    } else if (weights_scheme == "AIC") {
        wfun <- function(x) {
            aic_weighting(x)
        }
    } else {
        wfun <- function(x) {
            bic_weighting(x)
        }
    }
    parameter <- b <- forecast_dates <- NULL
    data <- object$y
    if (!is.null(object$xreg)) {
        use_xreg <- TRUE
        xreg <- xts(coredata(object$xreg), index(data))
    } else {
        use_xreg <- FALSE
        xreg <- NULL
    }
    start_date <- index(data)[start]
    end <- min(NROW(data), end)
    end_date <- index(data)[end - 1]
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    if (estimate_every != 1) {
        estimate_every <- max(1, as.integer(estimate_every))
        ns <- length(seqdates)
        seqdates <- seqdates[seq(1, ns, by = estimate_every)]
    }
    idates <- index(data)
    # check for ending period to avoid overlap
    index_table <- rbindlist(lapply(1:length(seqdates), function(i) {
        if (i == length(seqdates)) {
            s <- index(data[paste0(seqdates[i],"/",max(idates[end]))])
            s <- s[-length(s)]
        } else {
            s <- index(data[paste0(seqdates[i],"/",seqdates[i + 1])])
            s <- s[-length(s)]
        }
        if (length(s) == 0) return(NULL)
        rbindlist(lapply(1:length(s), function(j){
            idx <- which(idates == s[j])
            maxh <- min(idx + h, NROW(data[1:end]))
            tmp <- na.omit(index(data[idx:maxh]))
            data.table(estimate_date = seqdates[i], filter_date = tmp[1], forecast_dates = tmp[-1], h = length(tmp[-1]))
        }))
    }))
    max_index <- max(idates[end])
    index_table <- index_table[forecast_dates <= max_index]
    index_table <- split(index_table, by = "estimate_date")
    top_n <- object$top_n
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    i <- 1
    b <- lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        y_train <- data[paste0("/", seqdates[i])]
        if (use_xreg) {
            xreg_train <- xreg[index(y_train)]
        } else {
            xreg_train <- NULL
        }
        spec <- issm_modelspec(y_train, auto = TRUE, slope = object$slope, slope_damped = object$slope_damped,
                               seasonal = object$seasonal,
                               seasonal_frequency = object$seasonal_frequency,
                               seasonal_harmonics = object$seasonal_harmonics,
                               ar = object$ar, ma = object$ma, 
                               xreg = xreg_train, lambda = object$lambda, lower = object$lower, upper = object$upper, 
                               sampling = object$sampling, sample_n = object$sample_n, 
                               init_garch =  object$init_garch, garch_order = object$garch_order, 
                               variance = object$variance, distribution = object$distribution, top_n = top_n)
        mod <- try(estimate(spec, solver = solver, control = control, trace = FALSE), silent = TRUE)
        w <- wfun(mod)
        L <- index_table[[i]]
        M <- split(L, by = "filter_date")
        rolls <- length(M)
        # first index is the estimation so we start with the prediction
        # xreg
        if (use_xreg) {
            xreg_test <- xreg[M[[1]]$forecast_dates]
        } else {
            xreg_test <- NULL
        }
        # actual
        P <- vector(mode = "list", length = rolls)
        tmp <- predict(mod, h = M[[1]]$h[1], newvreg = xreg_test, forc_dates = M[[1]]$forecast_date, nsim = 2000, seed = seed)
        w_tmp <- tsensemble(tmp, weights = w, start = 1)
        if (rolls > 1 & rolling) {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            P[[1]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                 "filter_date" = M[[1]]$filter_date,
                                 "horizon" =  1:M[[1]]$h[1],
                                 "size" = rep(nrow(y_train), nh),
                                 "forecast_date" =  M[[1]]$forecast_dates,
                                 "forecast" = as.numeric(w_tmp$mean),
                                 "sigma" = apply(w_tmp$distribution, 2, sd),
                                 "p10" = apply(w_tmp$distribution, 2, quantile, 0.1),
                                 "p90" = apply(w_tmp$distribution, 2, quantile, 0.9),
                                 "min" = apply(w_tmp$distribution, 2, min),
                                 "max" = apply(w_tmp$distribution, 2, max),
                                 "actual" = as.numeric(y_test))
            for (j in 2:rolls) {
                roll_dates <- M[[j]]$filter_date[1]
                updated_y <- data[paste0(seqdates[i],"/",roll_dates)][-1,]
                if (use_xreg) {
                    newxreg <- xreg[paste0(seqdates[i],"/",roll_dates)][-1,]
                } else {
                    newxreg <- NULL
                }
                f <- tsfilter(mod, y = updated_y, newxreg = newxreg)
                if (use_xreg) {
                    xreg_test <- xreg[M[[j]]$forecast_dates]
                } else {
                    xreg_test <- NULL
                }
                tmp <- predict(f, h = M[[j]]$h[1], newxreg = xreg_test, forc_dates = M[[j]]$forecast_dates, nsim = 2000)
                w_tmp <- tsensemble(tmp, weights = w, start = 1)
                y_test <- data[M[[j]]$forecast_dates]
                nh <- length(M[[j]]$h)
                P[[j]] <- data.table("estimation_date" = rep(seqdates[i], nh),
                                     "filter_date" = M[[j]]$filter_date,
                                     "horizon" =  1:M[[j]]$h[1],
                                     "size" = rep(nrow(y_train), nh),
                                     "forecast_date" =  M[[j]]$forecast_dates,
                                     "forecast" = as.numeric(w_tmp$mean),
                                     "sigma" = apply(w_tmp$distribution, 2, sd),
                                     "p10" = apply(w_tmp$distribution, 2, quantile, 0.1),
                                     "p90" = apply(w_tmp$distribution, 2, quantile, 0.9),
                                     "min" = apply(w_tmp$distribution, 2, min),
                                     "max" = apply(w_tmp$distribution, 2, max),
                                     "actual" = as.numeric(y_test))
            }
            out <- rbindlist(P)
        } else {
            y_test <- data[M[[1]]$forecast_dates]
            nh <- length(M[[1]]$h)
            out <- data.table("estimation_date" = rep(seqdates[i], nh),
                              "filter_date" = M[[1]]$filter_date,
                              "horizon" =  1:M[[1]]$h[1],
                              "size" = rep(nrow(y_train), nh),
                              "forecast_date" =  M[[1]]$forecast_dates,
                              "forecast" = as.numeric(w_tmp$mean),
                              "sigma" = apply(w_tmp$distribution, 2, sd),
                              "p10" = apply(w_tmp$distribution, 2, quantile, 0.1),
                              "p90" = apply(w_tmp$distribution, 2, quantile, 0.9),
                              "min" = apply(w_tmp$distribution, 2, min),
                              "max" = apply(w_tmp$distribution, 2, max),
                              "actual" = as.numeric(y_test))
        }
        return(out)
    })
    b <- rbindlist(b)
    out <- list(table = b, h = h, estimate_every = estimate_every, rolling = rolling, type = "simple")
    class(out) <- "tsissm.backtest"
    return(out)
}

aic_weighting <- function(x)
{
    n <- length(x$models)
    aic <- AIC(x)
    if (any(is.na(aic))) {
        inc <- which(!is.na(aic))
    } else {
        inc <- 1:n
    }
    aic_w <- rep(0, n)
    aic_w[inc] <- aic[inc] - min(aic[inc])
    weights <- rep(0, n)
    weights[inc] <- exp(-0.5 * aic_w[inc])/sum(exp(-0.5 * aic_w[inc]))
    return(weights)
}

bic_weighting <- function(x)
{
    n <- length(x$models)
    bic <- BIC(x)
    if (any(is.na(bic))) {
        inc <- which(!is.na(bic))
    } else {
        inc <- 1:n
    }
    bic_w <- rep(0, n)
    bic_w[inc] <- bic[inc] - min(bic[inc])
    weights <- rep(0, n)
    weights[inc] <- exp(-0.5 * bic_w[inc])/sum(exp(-0.5 * bic_w[inc]))
    return(weights)
}

fixed_weights <- function(x, user_weights)
{
    n <- length(x$models)
    if (is.null(user_weights)) {
        user_weights <- rep(1/n, n)
    }
    bic <- BIC(x)
    if (any(is.na(bic))) {
        inc <- which(!is.na(bic))
    } else {
        inc <- 1:n
    }
    weights <- rep(0, n)
    weights[inc] <- user_weights[inc]
    weights <- weights/sum(weights)
    return(weights)
}