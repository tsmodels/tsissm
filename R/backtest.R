#' Walk Forward Model Backtest
#'
#' @description Generates an expanding window walk forward backtest.
#' @param object an object of class \dQuote{tsissm.spec}.
#' @param start numeric data index from which to start the backtest.
#' @param end numeric data index on which to end the backtest. The backtest will
#' end 1 period before that date in order to have at least 1 out of sample value
#' to compare against.
#' @param h forecast horizon. As the expanding window approaches the \dQuote{end},
#' the horizon will automatically shrink to the number of available out of sample
#' periods.
#' @param estimate_every number of periods at which the model is re-estimated
#' and new predictions are generated (defaults to 1).
#' @param FUN optional function which is applied across all horizons for each
#' draw (i.e. operating on each row of the distribution which represents a
#' single path). For example, using the max function will return the distribution
#' of the maximum across all horizons, whereas sum (for flow variables) would
#' represent the aggregate value distribution. The P50 of this distribution is
#' returned and aligned with the terminal horizon for each re-estimation period,
#' and if alpha is not NULL, then the quantiles of this distributions with
#' respect to the coverage (alpha) chosen.
#' @param alpha optional numeric vector of coverage rates for which to calculate
#' the quantiles.
#' @param solver solver to use.
#' @param autodiff whether to use automatic differentiation for estimation.
#' This makes use of the tsissmad package.
#' @param autoclean whether to perform automatic cleaning on the training data
#' prior to prediction as per the \sQuote{auto_clean} function in the tsaux package.
#' @param trace whether to show the progress bar. The user is expected to have
#' set up appropriate handlers for this using the \dQuote{progressr} package.
#' @param ... additional arguments passed to the \dQuote{auto_clean} function.
#' @return A list with the following data.tables:
#' \itemize{
#' \item prediction : the backtest table with forecasts and actuals
#' \item metrics: a summary performance table showing metrics by
#' forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
#' }
#' @note The function can use parallel functionality as long as the user has
#' set up a \code{\link[future]{plan}} using the future package.
#' @aliases tsbacktest
#' @method tsbacktest tsissm.spec
#' @rdname tsbacktest
#' @export
#'
#'
tsbacktest.tsissm.spec <- function(object, start = floor(length(object$target$y_orig)/2), end = length(object$target$y_orig),
                                   h = 1, estimate_every = 1, FUN = NULL, alpha = NULL, solver = "nlminb", 
                                   autodiff = TRUE, autoclean = FALSE, trace = FALSE, ...)
{
    if (object$seasonal$include_seasonal & object$seasonal$seasonal_type == "regular") {
        if (autodiff) {
            autodiff <- FALSE
            warning("\nautodiff only currently supported for trigonometric seasonality (switching to non autodiff)")
        }
    }
    data <- xts(object$target$y_orig, object$target$index)
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
    if ((start/2) < max(object$seasonal$seasonal_frequency)) {
        warning(paste0("\nmaximum seasonal frequency > 1/2 initial data length (based on start)"))
    }
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    if (estimate_every != 1) {
        estimate_every <- max(1, as.integer(estimate_every))
        ns <- length(seqdates)
        seqdates <- seqdates[seq(1, ns, by = estimate_every)]
    }
    elapsed_time <- function(idx, end_date, start_date) {
        min(h, which(end_date == idx) - which(start_date == idx))
    }
    if (!is.null(alpha)) {
        if (any(alpha <= 0)) {
            stop("\nalpha must be strictly positive")
        }
        if (any(alpha >= 1)) {
            stop("\nalpha must be less than 1")
        }
        quantiles <- as.vector(sapply(1:length(alpha), function(k) c(alpha[k]/2, 1 - alpha[k]/2)))
    } else {
        quantiles <- NULL
    }
    # setup backtest indices
    horizon <- sapply(1:length(seqdates), function(i){
        min(h, elapsed_time(index(data), index(data)[end], seqdates[i]))
    })
    i <- 1
   
    if (trace) {
        prog_trace <- progressor(length(seqdates))
    }
    
    extra_args <- list(...)
    if (length(extra_args) > 0 & any(names(extra_args) == "use_hessian")) {
        if (extra_args$use_hessian) {
            use_hessian <- TRUE
            extra_args$use_hessian <- NULL
        } else {
            use_hessian <- FALSE
        }
    } else {
        use_hessian <- FALSE
    }
    if (is.null(object$target$frequency)) {
        if (is.null(object$seasonal$seasonal_frequency[1])) {
            frequency <- 1
        } else {
            frequency <- object$seasonal$seasonal_frequency[1]
        }
    } else {
        frequency <- object$target$frequency
    }
    b %<-% future_lapply(1:length(seqdates), function(i) {
        if (trace) prog_trace()
        y_train <- data[paste0("/", seqdates[i])]
        ix <- which(index(data) == seqdates[i])
        y_test <- data[(ix + 1):(ix + horizon[i])]
        if (use_xreg) {
            xreg_train <- xreg[index(y_train)]
            xreg_test <- xreg[index(y_test)]
        } else {
            xreg_train <- NULL
            xreg_test <- NULL
        }
        if (object$transform$include_lambda) {
            lambda <- NA
        } else {
            lambda <- object$transform$lambda
        }
        if (autoclean) {
            if (object$transform$name == "box-cox") {
                if (is.na(lambda)) {
                    xlambda <- box_cox(lambda = NA)
                    xlambda <- attr(xlambda$transform(y_train, frequency = frequency),"lambda")
                } else {
                    xlambda <- lambda
                }
            } else {
                xlambda <- 1
            }
            args_x <- c(list(y = y_train), list(frequency = frequency), list(lambda = xlambda), extra_args)
            y_train <- do.call(auto_clean, args = args_x, quote = TRUE)
        }
        spec <- issm_modelspec(y_train, slope = object$slope$include_slope, slope_damped = object$slope$include_damped,
                               seasonal = object$seasonal$include_seasonal,
                               seasonal_frequency = object$seasonal$seasonal_frequency,
                               seasonal_type = object$seasonal$seasonal_type,
                               seasonal_harmonics = object$seasonal$seasonal_harmonics,
                               ar = object$arma$order[1], ma = object$arma$order[2], 
                               xreg = xreg_train, transformation = object$transform$name,
                               lambda = lambda, lower = object$transform$lower, upper = object$transform$upper, 
                               sampling = object$target$sampling)
        mod <- try(estimate(spec, solver = solver, autodiff = autodiff, use_hessian = use_hessian), silent = TRUE)
        if (inherits(mod, 'try-error') | is.null(mod)) {
            if (!is.null(quantiles)) {
                qp <- matrix(NA, ncol = length(quantiles), nrow = horizon[i])
                colnames(qp) <- paste0("P", round(quantiles*100,1))
            }
            out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]),
                              "horizon" = 1:horizon[i],
                              "size" = rep(nrow(y_train), horizon[i]),
                              "forecast_dates" = as.character(index(y_test)),
                              "forecast" = rep(as.numeric(NA), horizon[i]), "actual" = as.numeric(y_test))
            if (!is.null(quantiles)) out <- cbind(out, qp)
            if (!is.null(FUN)) {
                funp <- data.table(estimation_date = seqdates[i], horizon = horizon[i], fun_P50 = as.numeric(NA), fun_actual = FUN(as.numeric(y_test)))
                if (!is.null(quantiles)) {
                    qp <- matrix(rep(as.numeric(NA)), length(quantiles), nrow = 1)
                    colnames(qp) <- paste0("fun_P", round(quantiles*100,1))
                    funp <- cbind(funp, qp)
                }
                out <- merge(out, funp, by = c("estimation_date","horizon"), all.x = T)
            }
            return(out)
        } else {
            p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(y_test))
            if (!is.null(quantiles)) {
                qp <- apply(p$distribution, 2, quantile, quantiles)
                if (length(quantiles) == 1) {
                    qp <- matrix(qp, ncol = 1)
                } else{
                    qp <- t(qp)
                }
                colnames(qp) <- paste0("P", round(quantiles*100,1))
            }
            out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]),
                              "horizon" = 1:horizon[i],
                              "size" = rep(nrow(y_train), horizon[i]),
                              "forecast_dates" = as.character(index(y_test)),
                              "forecast" = as.numeric(p$mean), "actual" = as.numeric(y_test))
            if (!is.null(quantiles)) out <- cbind(out, qp)
            if (!is.null(FUN)) {
                pd <- apply(p$distribution, 1, FUN)
                funp <- data.table(estimation_date = seqdates[i], horizon = horizon[i], fun_P50 = median(pd), fun_actual = FUN(as.numeric(y_test)))
                if (!is.null(quantiles)) {
                    qp <- matrix(quantile(pd, quantiles), nrow = 1)
                    colnames(qp) <- paste0("fun_P", round(quantiles*100,1))
                    funp <- cbind(funp, qp)
                }
                out <- merge(out, funp, by = c("estimation_date","horizon"), all = TRUE)
            }
            return(out)
        }
    }, future.packages = c("tsmethods","tsaux","xts","tsissm","data.table"), future.seed = TRUE)
    b <- eval(b)
    b <- rbindlist(b)
    if (is.null(data_name)) data_name <- "y"
    actual <- NULL
    forecast <- NULL
    z <- copy(b)
    z <- na.omit(z)
    metrics <- z[,list(variable = data_name, MAPE = mape(actual, forecast), MSLRE = mslre(actual, forecast), BIAS = bias(actual, forecast), n = .N), by = "horizon"]
    if (!is.null(alpha)) {
        q_names <- matrix(paste0("P", round(quantiles*100,1)), ncol = 2, byrow = TRUE)
        q <- do.call(cbind, lapply(1:length(alpha), function(i){
            z[,list(mis = mis(actual, get(q_names[i,1]), get(q_names[i,2]), alpha[i])), by = "horizon"]
        }))
        q <- q[,which(grepl("mis",colnames(q))), with = FALSE]
        colnames(q) <- paste0("MIS[",alpha,"]")
        metrics <- cbind(metrics, q)
    }
    return(list(prediction = b, metrics = metrics))
}
