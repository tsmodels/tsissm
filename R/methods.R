# fitted  ---------------------------------------------------

#' Model Fitted Values
#'
#' @description Extract the fitted values from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @returns an xts object of the fitted values.
#' @aliases fitted
#' @method fitted tsissm.estimate
#' @rdname fitted
#' @export
#'
#'
fitted.tsissm.estimate <- function(object, ...)
{
    return(xts(object$model$fitted, object$spec$target$index))
}


# residuals ---------------------------------------------------

#' Model Residuals
#'
#' @description Extract the residual values from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param transformed residuals based values in transformed space (Box Cox).
#' @param standardize whether to scale the residuals by the estimated standard deviation.
#' @param ... not currently used.
#' @details For h>1, this is like performing an in-sample backtest starting at
#' time 1 with fixed coefficients. The purpose of having the matrix of h-step
#' ahead residuals is in order to calculate the 1:h covariance matrix as well
#' as the cross 1:h covariance matrix when ensembling series at multiple horizons.
#' @return An xts vector of the model residuals for h = 1, else a data.table
#' with rows representing the first prediction date and columns the h-ahead
#' forecast residuals.
#' @note The function can use parallel functionality (for h>1) as long as the
#' user has set up a \code{\link[future]{plan}} using the future package.
#' @aliases residuals
#' @method residuals tsissm.estimate
#' @rdname residuals
#' @export
#'
#'
residuals.tsissm.estimate <- function(object, standardize = FALSE, transformed = TRUE, ...)
{
    if (transformed) {
        out <- xts(object$model$error, object$spec$target$index)
        if (standardize) out <- out/sigma(object)
    } else {
        out <- xts(object$model$residuals, object$spec$target$index)
        if (standardize) out <- scale(out, center = FALSE, scale = TRUE)
    }
    return(out)
}

#' Multi-Step Ahead In-Sample Residuals
#'
#' @description Extract the multi-step ahead in-sample residual values from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param h the forecast horizon
#' @param transformed residuals based values in transformed space (Box Cox).
#' @param index_start the time point from which to initiate the in-sample rolling 
#' forecasts. This is zero based to enable the first forecast to be t=1.
#' @param simplify whether to return a matrix type data.table of error against
#' date and horizon, else the long for data.table with the forecasts, actuals and
#' errors.
#' @param ... not currently used.
#' @details For each time point t (t>=index_start), in the data sample, an h-steps ahead 
#' forecast  (predicting the observation at time t + h) is made, using the full sample 
#' estimated parameters and observed data up to t. These h-step-ahead fitted residuals,
#' in either transformed or untransformed space, can sometimes be used for diagnosing
#' the multi-step ahead in-sample performance of the model. This is not a substitute 
#' for a proper rolling out of sample forecast, but a quick method which may still
#' be useful.
#' @returns A data.table in either long or wide format.
#' @aliases hresiduals
#' @method hresiduals tsissm.estimate
#' @rdname hresiduals
#' @export
#'
#'
hresiduals.tsissm.estimate <- function(object, h = 12, transformed = TRUE, index_start = 0, simplify = TRUE, ...)
{
    n <- length(object$spec$target$y)
    if (index_start < 0 | index_start >= n) stop("\nindex_start must be positive and less than length of dataset.")
    object$model$states <- rbind(object$model$xseed, object$model$states)
    # since we use zero indexing we increment index_start and append a zero to each vector/matrix passed
    index_start <- index_start + 1
    n <- n + 1
    fun <- hresiduals_issm
    h_residuals <- lapply(index_start:(n - 1), function(i) {
        tmp <- fun(object, h = h, nsim = 2000, start = i, transformed = transformed)
        return(tmp)
    })
    h_residuals <- rbindlist(h_residuals)
    if (simplify) {
        h_residuals <- dcast(h_residuals, date~horizon, value.var = "error")
    }
    # error <- c(0, object$model$error)
    # test all.equal(h_residuals$`1`, error[-c(1:index_start)])
    return(h_residuals)
}

#' @aliases hresiduals
#' @rdname hresiduals
#' @export
#'
hresiduals <- function(object, ...)
{
    UseMethod("hresiduals")
}



# summary ---------------------------------------------------

.make_standard_errors <- function(object, vcov_type = "QMLE")
{
    pmatrix <- object$parmatrix
    pars <- pmatrix[estimate == 1]$value
    V <- vcov(object, type = vcov_type)
    se <- suppressWarnings(sqrt(diag(V)))
    tvalues <- pars/se
    pvalues <- 2*(1 - pnorm(abs(tvalues)))
    return(data.frame("Std. Error" = se,"t value" = tvalues, "Pr(>|t|)" = pvalues, check.names = FALSE))
}

.make_model_description <- function(object) {
    model <- "A"
    if (object$spec$slope$include_slope) {
        model <- c(model,"A")
            if (object$spec$slope$include_damped) {
                model <- c(model,"d")                
            }
    } else {
        model = c(model, "N")
    }
    if (object$spec$seasonal$include_seasonal) {
        model <- c(model, "A")
        if (object$spec$seasonal$seasonal_type == "trigonometric") {
            s <- sapply(1:length(object$spec$seasonal$seasonal_frequency), function(i) {
                if (i < length(object$spec$seasonal$seasonal_frequency)) {
                    paste0(object$spec$seasonal$seasonal_frequency[i],"{",object$spec$seasonal$seasonal_harmonics[i],"}/")
                } else {
                    paste0(object$spec$seasonal$seasonal_frequency[i],"{",object$spec$seasonal$seasonal_harmonics[i],"}")
                }
            })
        } else {
            s <- sapply(1:length(object$spec$seasonal$seasonal_frequency), function(i) {
                paste0(object$spec$seasonal$seasonal_frequency[i],"{}")
            })
        }
        model <- c(model,"(",s,")")
    } else {
        model <- c(model,"N")
    }
    model <- c(model,"+ARMA(",object$spec$arma$order[1],",",object$spec$arma$order[2],")")
    if (object$spec$xreg$include_xreg) {
        model <- c(model,"+X(",NCOL(object$spec$xreg$xreg),")")
    }
    if (object$spec$variance$type == "dynamic") {
        model <- c(model, "+GARCH(",object$spec$variance$order[1],",",object$spec$variance$order[1],")")
    }
    model <- paste0(model,collapse = "")
    return(model)
}

.distribution_name <- function(x) {
    switch(x, 
           "norm" = "Gaussian",
           "std" = "Student's T",
           "jsu" = "Johnson SU")
}

#' Model Estimation Summary
#'
#' @description Summary method for class \dQuote{tsissm.estimate}
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param vcov_type the type of standard errors based on the vcov estimate (see \code{\link{vcov}}).
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param ... not currently used.
#' @return A printout of the parameter summary, model type and some model metrics.
#' @aliases summary
#' @method summary tsissm.estimate
#' @rdname summary
#' @export
#'
#'
#'
summary.tsissm.estimate <- function(object, vcov_type = "H", digits = 4, ...)
{
    estimate <- NULL
    V <- vcov(object, type = vcov_type)
    est <- object$parmatrix[estimate == 1]$value
    par_names <- object$parmatrix[estimate == 1]$parameters
    se <- sqrt(diag(V))
    tval <- est/se
    coefficients <- cbind(Estimate = est, `Std. Error` = se,`t value` = tval, `Pr(>|t|)` = 2*(1 - pnorm(abs(tval))))
    n_obs <- nobs(object)
    n_missing <- length(which(object$spec$good == 0))
    n_parameters <- length(coef(object))
    llh <- -object$model$loglik
    elapsed <- object$opt$elapsed
    #conditions <- object$conditions[c("kkt1","kkt2","evratio")]
    distribution <- object$spec$distribution$type
    equation <- tsequation(object)
    coefficients <- as.data.table(coefficients, keep.rownames = TRUE)
    init_variance <- ifelse(object$spec$variance$type == "constant", object$model$sigma, object$model$initial_variance)
    setnames(coefficients, "rn","term")
    syms <- object$parmatrix[estimate == 1]$symbol
    y_actual <- valid_data(object$spec$target$y_orig, object$spec$good)
    y_fitted <- valid_data(fitted(object), object$spec$good)
    y_change <- na.omit(diff(y_actual))
    f_change <- na.omit(diff(y_fitted))
    DAC <- accuracy(y_change, f_change)
    MAPE <- mape(y_actual, y_fitted) * 100
    model <- .make_model_description(object)
    init_states <- object$model$xseed
    
    out <- list(coefficients = coefficients, distribution = distribution,
                loglikelihood = llh, n_obs = n_obs, n_missing = n_missing, n_parameters = n_parameters,
                AIC = AIC(object),
                BIC = BIC(object),
                MAPE = MAPE,
                DAC = DAC, 
                init_variance = init_variance,
                init_states = init_states,
                init_method = object$spec$model$init,
                elapsed = elapsed, conditions = NULL, equation = equation,
                model = model, symbol = syms,
                variance_type = object$spec$variance$type,
                equation = object$parmatrix[estimate == 1]$equation)
    class(out) <- "summary.tsissm.estimate"
    return(out)
}

#' Model Estimation Summary Print method
#'
#' @description Print method for class \dQuote{summary.tsissm.estimate}
#' @param x an object of class \dQuote{summary.tsissm.estimate}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param ... not currently used.
#' @return Invisibly returns the original summary object.
#' @aliases print.summary.tsissm.estimate
#' @method print summary.tsissm.estimate
#' @rdname print
#' @export
#'
#'
print.summary.tsissm.estimate <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...)
{
    .print_screen(x, digits = digits, signif.stars = signif.stars, ...)
}


#' Transform a summary object into flextable
#' @description
#' Transforms a \dQuote{summary.tsissm.estimate} object into a flextable
#' with options on symbolic representation and model equation.
#' @param x an object of class \dQuote{summary.tsissm.estimate}.
#' @param digits integer, used for number formatting. Optionally, to avoid
#' scientific notation, set \sQuote{options(scipen=999)}.
#' @param signif.stars logical. If TRUE, ‘significance stars’ are printed for each coefficient.
#' @param include.symbols logical. If TRUE, replaces parameter names with their symbols (if they exist).
#' @param include.equation logical. If TRUE, adds a section with the symbolic model equation.
#' @param include.statistics logical. If TRUE, adds a section with summary statistics on the model.
#' @param table.caption an optional string for the table caption.
#' @param ... additional arguments passed to flextable method.
#' @importFrom flextable as_flextable
#' @return A flextable object.
#' @aliases as_flextable.summary.tsissm.estimate
#' @method as_flextable summary.tsissm.estimate
#' @rdname as_flextable.summary
#' @export
as_flextable.summary.tsissm.estimate <- function(x, digits = max(3L, getOption("digits") - 3L),
                                                  signif.stars = getOption("show.signif.stars"),
                                                  include.symbols = TRUE, include.equation = TRUE,
                                                  include.statistics = TRUE,
                                                  table.caption = paste0("ISSM Model: ", x$model), ...)
{
    out <- .print_flextable(x, digits = digits, signif.stars = signif.stars,
                            include.symbols = include.symbols, include.equation = include.equation,
                            include.statistics = include.statistics,
                            table.caption = table.caption, ...)
    return(out)
}


.resid.variance <- function(object)
{
    parameters <- NULL
    td <- tsdecompose(object)
    snames <- colnames(td)
    n <- length(snames)
    v <- rep(0, ncol(td))
    r <- xts(matrix(0, ncol = ncol(td), nrow = nrow(td)), object$spec$target$index)
    A <- xts(object$spec$transform$transform(object$spec$target$y_orig, lambda = object$parmatrix[parameters == "lambda"]$value), object$spec$target$index)
    r[,1] = A - td[,1]
    v[1] <- var(r[,1])
    if (n > 1) {
        for (i in 2:n) {
            v[i] <- var(r[,i - 1] - td[,i])
        }
    }
    names(v) <- c("V[Actual - Level]", paste0("V[...-",snames[-1],"])"))
    attr(v, "states") <- snames
    return(v)
}

# coef ---------------------------------------------------

#' Extract Model Coefficients
#'
#' @description Extract the estimated coefficients of a model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @return a numeric named vector.
#' @aliases coef
#' @method coef tsissm.estimate
#' @rdname coef
#' @export
#'
#'
coef.tsissm.estimate <- function(object, ...)
{
    estimate <- NULL
    return(structure(object$parmatrix[estimate == 1]$value, names = object$parmatrix[estimate == 1]$parameters))
}

# tsdecompose ---------------------------------------------------

#' Model Decomposition
#'
#' @description Decomposes the estimated model or prediction into its component
#' parts (states).
#' @param object an object of class \dQuote{tsissm.estimate} or \dQuote{tsissm.predict}
#' @param simplify simplification of the returned states aggregates the level 
#' and slope (if present) into a Trend component, all Seasonal components, all 
#' ARMA components and the error terms into an Irregular component. This may be 
#' useful when bagging  is carried out (assuming equal lambda in the box-cox 
#' transform). This also simplifies the ability to created custom overrides of 
#' the Trend and rebuilt the predictive distribution.
#' @param start whether to return the predicted states from t=1 to h or the
#' states from t=0 to (h-1). The latter is sometimes useful as the sum of the
#' states equals the predicted value (since the predictions are based on the lagged
#' state).
#' @param ... not currently used.
#' @return For the estimated object, returns an xts matrix of the state components  
#' (including Irregular term). For the predicted object, a list with the simulated 
#' state components of class \dQuote{tsmodel.predict} which includes the predictive 
#' distribution and original (estimated) series components.
#' @details The 1-step ahead prediction is given by the following equation:
#' \deqn{y_{t} = x_{t-1} w + \varepsilon_{t}} 
#' Because the decomposition pre lags the states so that the seed state is aligned 
#' with the error term, then summing the state distribution of each component with 
#' the returned error distribution will ensure that the exact same predicted 
#' value matched to the correct date is returned.
#' @aliases tsdecompose
#' @method tsdecompose tsissm.estimate
#' @rdname tsdecompose
#' @export
#'
#'
tsdecompose.tsissm.estimate <- function(object, simplify = FALSE, start = 1, ...)
{
    estimate <- NULL
    idx <- object$spec$idmatrix
    rnames <- rownames(idx)
    states <- rbind(object$model$xseed, object$model$states)
    indx <- object$spec$target$index
    pars <- object$parmatrix[estimate == 1]$value
    names(pars) <- object$parmatrix[estimate == 1]$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    mdim = object$spec$dims
    w <- matrix(S[list("w")]$values, nrow = mdim[1], ncol = 1)
    w_t <- t(w)
    if (simplify) {
        zero_matrix <- xts(rep(0, length(indx)), indx)
        Trend <- S <- ARMA <- Irregular <- X <- zero_matrix
        colnames(Trend) <- "Trend"
        colnames(S) <- "Seasonal"
        colnames(ARMA) <- "ARMA"
        colnames(Irregular) <- "Irregular"
        colnames(X) <- "X"
    }
    if (start == 1) {
        state_index <- 2:dim(states)[1]
    } else {
        state_index <- 1:(dim(states)[1] - 1)
    }
    
    if (idx["Level","n"] > 0) {
        Level <- xts(states[state_index, idx["Level","start"]:idx["Level","end"]], indx)
        if (simplify) Trend <- Trend + Level
    } else {
        Level <- NULL
    }
    if (idx["Slope","n"] > 0) {
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- xts(w_t[,nstart:nend] * states[state_index, idx["Slope","start"]:idx["Slope","end"]], indx)
        colnames(Slope) <- "Slope"
        if (simplify) Trend <- Trend + Slope
    } else {
        Slope <- NULL
    }
    if (any(idx[grepl("Seasonal",rnames),"n"] > 0)) {
        ns <- length(idx[grepl("Seasonal",rnames),"n"])
        nstart <- idx[grepl("Seasonal",rnames),"start"]
        nend <- idx[grepl("Seasonal",rnames),"end"]
        frequency <- object$spec$seasonal$seasonal_frequency
        Seasonal <- matrix(0, ncol = ns, nrow = length(indx))
        colnames(Seasonal) <- names(nstart)
        if (object$spec$seasonal$seasonal_type == "trigonometric") {
            K <- object$spec$seasonal$seasonal_harmonics
            for (i in 1:ns) {
                Seasonal[,i] <- t(w_t[,nstart[i]:(nstart[i] + K[i] - 1)] %*% t(states[state_index, nstart[i]:(nstart[i] + K[i] - 1)]))
            }
        } else {
            for (i in 1:ns) {
                Seasonal[,i] <- t(w_t[,nstart[i]:(nstart[i] + frequency[i] - 1)] %*% t(states[state_index, nstart[i]:(nstart[i] + frequency[i] - 1)]))
            }
        }
        Seasonal <- xts(Seasonal, indx)
        if (simplify) S <- S + xts(rowSums(Seasonal), indx)
    } else {
        Seasonal <- S <- NULL
    }
    if (idx["AR","n"] > 0) {
        AR <- t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[state_index, idx["AR","start"]:idx["AR","end"]]))
        AR <- xts(AR, indx)
        colnames(AR) <- paste0("AR",idx["AR","n"])
        if (simplify) ARMA <- ARMA + AR
    } else {
        AR <- NULL
    }
    if (idx["MA","n"] > 0) {
        MA <- t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[state_index, idx["MA","start"]:idx["MA","end"]]))
        MA <- xts(MA, indx)
        colnames(MA) <- paste0("MA",idx["MA","n"])
        if (simplify) ARMA <- ARMA + MA
    } else {
        MA <- NULL
    }
    if (idx["AR","n"] == 0 & idx["MA","n"] == 0) {
        ARMA <- NULL
    }
    if (idx["X","n"] > 0) {
        beta <- matrix(object$parmatrix[which(grepl("kappa", object$parmatrix$parameters))]$value, ncol = 1)
        X <- object$spec$xreg$xreg
        X <- xts(X %*% beta, indx)
        colnames(X) <- "X"
    } else {
        X <- NULL
    }
    if (object$spec$variance$type == "dynamic") {
        Sigma <- sigma(object)
    } else {
        Sigma <- NULL
    }
    Irregular <- residuals(object, transformed = TRUE)
    colnames(Irregular) <- "Irregular"
    if (simplify) {
        decomposition <- cbind(Trend, S, ARMA, X, Irregular, Sigma)
    } else {
        decomposition <- cbind(Level, Slope, Seasonal, AR, MA, X, Irregular, Sigma)
    }
    return(decomposition)
}

#' @method tsdecompose tsissm.predict
#' @rdname tsdecompose
#' @export
tsdecompose.tsissm.predict <- function(object, simplify = FALSE, start = 1, ...)
{
    estimate <- NULL
    idx <- object$spec$idmatrix
    rnames <- rownames(idx)
    states <- object$states
    indx <- object$spec$target$index
    pars <- object$spec$parmatrix[estimate == 1]$value
    names(pars) <- object$spec$parmatrix[estimate == 1]$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    mdim <- object$spec$dims
    nsim <- nrow(object$distribution)
    fdates <- colnames(object$distribution)
    date_class <- attr(object$distribution, "date_class")
    L <- list()
    k <- 1
    w <- matrix(S[list("w")]$values, nrow = mdim[1], ncol = 1)
    w_t <- t(w)
    
    if (start == 1) {
        state_index <- 2:dim(states)[1]
    } else {
        state_index <- 1:(dim(states)[1] - 1)
    }
    
    d_names <- NULL
    s_names <- NULL
    if (simplify) {
        zero_matrix <- matrix(0, nrow = nsim, ncol = object$h)
        empty_list <- list(distribution = zero_matrix, original_series = xts(rep(0, NROW(object$original_series)), index(object$original_series)))
        Trend <- Seasonal <- ARMA <- Irregular <- X <- empty_list
    }
    if (idx["Level","n"] > 0) {
        d_names <- c(d_names, "Level")
        Level <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(states[state_index, idx["Level","start"]:idx["Level","end"], i],nrow = 1)
        }))
        # start at time zero since y[t] = x[t-1] so if we want to replicate we need this
        # Level <- Level[,1:(ncol(Level) - 1)]
        colnames(Level) <- fdates
        class(Level) <- "tsmodel.distribution"
        attr(Level, "date_class") <- date_class
        Level <- list(original_series = object$decomp$Level, distribution = Level)
        class(Level) <- "tsmodel.predict"
        L[[1]] <- Level
        k <- k + 1
        if (simplify) {
            s_names <- c(s_names, "Trend")
            Trend <- Level
        }
    }
    if (idx["Slope","n"] > 0) {
        d_names <- c(d_names, "Slope")
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(w_t[,nstart:nend] * states[state_index, idx["Slope","start"]:idx["Slope","end"], i], nrow = 1)
        }))
        colnames(Slope) <- fdates
        class(Slope) <- "tsmodel.distribution"
        attr(Slope, "date_class") <- date_class
        Slope <- list(original_series = object$decomp$Slope, distribution = Slope)
        class(Slope) <- "tsmodel.predict"
        L[[k]] <- Slope
        k <- k + 1
        if (simplify) {
            Trend$distribution <- Trend$distribution + Slope$distribution
            Trend$original_series <- Trend$original_series + Slope$original_series
        }
    }

    if (any(idx[grepl("Seasonal",rnames),"n"] > 0)) {
        ns <- length(idx[grepl("Seasonal",rnames),"n"])
        nstart <- idx[grepl("Seasonal",rnames),"start"]
        nend <- idx[grepl("Seasonal",rnames),"end"]
        frequency <- object$spec$seasonal$seasonal_frequency
        K <- object$spec$seasonal$seasonal_harmonics
        snames <- rnames[grepl("^Seasonal", rnames)]
        for (j in 1:ns) {
            d_names <- c(d_names, snames[j])
            tmp <- do.call(rbind, lapply(1:nsim, function(i){
                matrix(t(w_t[,nstart[j]:(nstart[j] + K[j] - 1)] %*% t(.retain_dimensions_array(states, i)[state_index, nstart[j]:(nstart[j] + K[j] - 1), drop = FALSE])), nrow = 1)
            }))
            colnames(tmp) <- fdates
            class(tmp) <- "tsmodel.distribution"
            attr(tmp, "date_class") <- date_class
            tmp <- list(original_series = object$decomp[,k], distribution = tmp)
            class(tmp) <- "tsmodel.predict"
            L[[k]] <- tmp
            if (simplify) {
                Seasonal$distribution <- Seasonal$distribution + tmp$distribution
                Seasonal$original_series <- Seasonal$original_series + tmp$original_series
            }
            k <- k + 1
        }
        if (simplify) s_names <- c(s_names, "Seasonal")
    } else {
        Seasonal <- NULL
    }
    if (idx["AR","n"] > 0) {
        d_names <- c(d_names, "AR")
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(.retain_dimensions_array(states, i)[state_index, idx["AR","start"]:idx["AR","end"], drop = FALSE])), nrow = 1)
        }))
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        tmp <- list(original_series = object$decomp[,k], distribution = tmp)
        class(tmp) <- "tsmodel.predict"
        L[[k]] <- tmp
        k <- k + 1
        if (simplify) {
            ARMA$distribution <- ARMA$distribution + tmp$distribution
            ARMA$original_series <- ARMA$original_series + tmp$original_series
        }
    }
    if (idx["MA","n"] > 0) {
        d_names <- c(d_names, "MA")
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(.retain_dimensions_array(states, i)[state_index, idx["MA","start"]:idx["MA","end"], drop = FALSE])), nrow = 1)
        }))
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        tmp <- list(original_series = object$decomp[,k], distribution = tmp)
        class(tmp) <- "tsmodel.predict"
        L[[k]] <- tmp
        k <- k + 1
        if (simplify) {
            ARMA$distribution <- ARMA$distribution + tmp$distribution
            ARMA$original_series <- ARMA$original_series + tmp$original_series
        }
    }
    if (idx["AR","n"] == 0 & idx["MA","n"] == 0) {
        ARMA <- NULL
    } else if (idx["AR","n"] > 0 | idx["MA","n"] > 0) {
        s_names <- c(s_names, "ARMA")
    }
    
    if (idx["X","n"] > 0) {
        d_names <- c(d_names, "X")
        s_names <- c(s_names, "X")
        beta <- matrix(object$spec$parmatrix[which(grepl("kappa", object$spec$parmatrix$parameters))]$value, ncol = 1)
        X <- object$spec$xreg$xreg
        xreg <- xts(X %*% beta, indx)
        colnames(xreg) <- "X"
        tmp <- coredata(object$newxreg) %*% beta
        tmp <- matrix(as.numeric(tmp), ncol = object$h, nrow = nsim, byrow = TRUE)
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        xtmp <- list(original_series = xreg, distribution = tmp)
        class(xtmp) <- "tsmodel.predict"
        L[[k]] <- xtmp
        k <- k + 1
    } else {
        X <- NULL
    }
    
    # Innovations
    Irregular <- object$innov
    L[[k]] <- Irregular
    s_names <- c(s_names, "Irregular")
    d_names <- c(d_names, "Irregular")
    
    names(L) <- d_names
    if (simplify) {
        if (!is.null(Seasonal)) class(Seasonal) <- "tsmodel.predict"
        if (!is.null(ARMA)) class(ARMA) <- "tsmodel.predict"
        class(Irregular) <- "tsmodel.predict"
        L <- list()
        L$Trend <- Trend
        if (!is.null(Seasonal)) L$Seasonal <- Seasonal
        if (!is.null(ARMA)) L$ARMA <- ARMA
        if (!is.null((X))) L$X <- X
        L$Irregular <- Irregular
    }
    return(L)
}


#' @method tsdecompose tsissm.simulate
#' @rdname tsdecompose
#' @export
tsdecompose.tsissm.simulate <- function(object, simplify = FALSE, start = 1, ...)
{
    estimate <- NULL
    idx <- object$spec$idmatrix
    rnames <- rownames(idx)
    states <- object$states
    indx <- object$spec$target$index
    pars <- object$spec$parmatrix[estimate == 1]$value
    names(pars) <- object$spec$parmatrix[estimate == 1]$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    mdim <- object$spec$dims
    nsim <- NROW(object$distribution)
    h <- NCOL(object$distribution)
    fdates <- colnames(object$distribution)
    date_class <- attr(object$distribution, "date_class")
    L <- list()
    k <- 1
    w <- matrix(S[list("w")]$values, nrow = mdim[1], ncol = 1)
    w_t <- t(w)
    
    if (start == 1) {
        state_index <- 2:NROW(states[[1]])
    } else {
        state_index <- 1:(NROW(states[[1]]) - 1)
    }
    d_names <- NULL
    if (simplify) {
        zero_matrix <- matrix(0, nrow = nsim, ncol = h)
        empty_list <- zero_matrix
        Trend <- Seasonal <- ARMA <- Irregular <- empty_list
    }
    if (idx["Level","n"] > 0) {
        d_names <- c(d_names, "Level")
        Level <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(states[[i]][state_index, idx["Level","start"]:idx["Level","end"]], nrow = 1)
        }))
        # start at time zero since y[t] = x[t-1] so if we want to replicate we need this
        # Level <- Level[,1:(ncol(Level) - 1)]
        colnames(Level) <- fdates
        class(Level) <- "tsmodel.distribution"
        attr(Level, "date_class") <- date_class
        L[[1]] <- Level
        k <- k + 1
        if (simplify) {
            Trend <- Level
        }
    }
    if (idx["Slope","n"] > 0) {
        d_names <- c(d_names, "Slope")
        nstart <- idx[grepl("Slope",rnames),"start"]
        nend <- idx[grepl("Slope",rnames),"end"]
        # w.t will be either 1 (not dampening) else the dampening parameter
        Slope <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(w_t[,nstart:nend] * states[[i]][state_index, idx["Slope","start"]:idx["Slope","end"]], nrow = 1)
        }))
        #  Slope <- Slope[,1:(ncol(Slope) - 1)]
        colnames(Slope) <- fdates
        class(Slope) <- "tsmodel.distribution"
        attr(Slope, "date_class") <- date_class
        L[[k]] <- Slope
        k <- k + 1
        if (simplify) {
            Trend <- Trend + Slope
        }
    }
    if (any(idx[grepl("Seasonal",rnames),"n"] > 0)) {
        ns <- length(idx[grepl("Seasonal",rnames),"n"])
        nstart <- idx[grepl("Seasonal",rnames),"start"]
        nend <- idx[grepl("Seasonal",rnames),"end"]
        frequency <- object$spec$seasonal$seasonal_frequency
        K <- object$spec$seasonal$seasonal_harmonics
        snames <- rnames[grepl("Seasonal",rnames)]
        for (j in 1:ns) {
            d_names <- c(d_names, snames[j])
            tmp <- do.call(rbind, lapply(1:nsim, function(i){
                matrix(t(w_t[,nstart[j]:(nstart[j] + K[j] - 1)] %*% t(states[[i]][state_index, nstart[j]:(nstart[j] + K[j] - 1), drop = FALSE])), nrow = 1)
            }))
            # tmp <- tmp[,1:(ncol(tmp) - 1)]
            colnames(tmp) <- fdates
            class(tmp) <- "tsmodel.distribution"
            attr(tmp, "date_class") <- date_class
            L[[k]] <- tmp
            if (simplify) {
                Seasonal <- Seasonal + tmp
            }
            k <- k + 1
        }
    } else {
        Seasonal <- NULL
    }
    if (idx["AR","n"] > 0) {
        d_names <- c(d_names, "AR")
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["AR","start"]:idx["AR","end"]] %*% t(states[[i]][state_index, idx["AR","start"]:idx["AR","end"], drop = FALSE])), nrow = 1)
        }))
        # tmp <- tmp[,1:(ncol(tmp) - 1)]
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        L[[k]] <- tmp
        k <- k + 1
        if (simplify) {
            ARMA <- ARMA + tmp
        }
    }
    if (idx["MA","n"] > 0) {
        d_names <- c(d_names, "MA")
        tmp <- do.call(rbind, lapply(1:nsim, function(i){
            matrix(t(w_t[,idx["MA","start"]:idx["MA","end"]] %*% t(states[[i]][state_index, idx["MA","start"]:idx["MA","end"], drop = FALSE])), nrow = 1)
        }))
        colnames(tmp) <- fdates
        class(tmp) <- "tsmodel.distribution"
        attr(tmp, "date_class") <- date_class
        L[[k]] <- tmp
        k <- k + 1
        if (simplify) {
            ARMA <- ARMA + tmp
        }
    }
    if (idx["AR","n"] == 0 & idx["MA","n"] == 0) {
        ARMA <- NULL
    }
    # Innovations
    d_names <- c(d_names, "Irregular")
    Irregular <- object$E[,-1, drop = FALSE]
    colnames(Irregular) <- fdates
    class(Irregular) <- "tsmodel.distribution"
    attr(Irregular, "date_class") <- date_class
    L[[k]] <- Irregular
    names(L) <- d_names
    if (simplify) {
        L <- list()
        L$Trend <- Trend
        if (!is.null(Seasonal)) L$Seasonal <- Seasonal
        if (!is.null(ARMA)) L$ARMA <- ARMA
        L$Irregular <- Irregular
    }
    return(L)
}


# logLik ---------------------------------------------------

#' Model Log-Likelihood
#'
#' @description Extract the log-likelihood from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @return Returns an object of class logLik. This is a number with at least one
#' attribute, "df" (degrees of freedom), giving the number of (estimated)
#' parameters in the model.
#' @aliases logLik
#' @method logLik tsissm.estimate
#' @rdname logLik
#' @export
#'
#'
logLik.tsissm.estimate <- function(object, ...)
{
    estimate <- NULL
    # estimated parameters + estimated states + unconditional variance (or sigma for constant model)
    np <- NROW(object$parmatrix[estimate == 1]) + length(object$model$xseed) + 1
    loglik <- -1.0 * object$model$loglik
    return(structure(loglik, df = np, class = "logLik"))
}

#' @method logLik tsissm.selection
#' @rdname logLik
#' @export
#'
#'
logLik.tsissm.selection <- function(object, ...)
{
    sapply(object$models, function(x) logLik(x))
}


# nobs ---------------------------------------------------

#' Extract the Number of Observations
#'
#' @description Extract the number of observations from an estimated model.
#' This is principally intended to be used in computing BIC and used in other
#' tidy methods
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases nobs
#' @method nobs tsissm.estimate
#' @rdname nobs
#' @export
#'
#'
nobs.tsissm.estimate <- function(object, ...)
{
    return(length(object$spec$target$y_orig))
}



# AIC/BIC ---------------------------------------------------

#' Akaike's An Information Criterion
#'
#' @description Extract the AIC from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @param k the penalty per parameter to be used; the default k = 2 is the
#' classical AIC.
#' @return a numeric value.
#' @aliases AIC
#' @method AIC tsissm.estimate
#' @rdname AIC
#' @export
#'
#'
AIC.tsissm.estimate <- function(object, ..., k = 2)
{
    estimate <- NULL
    np <- NROW(object$parmatrix[estimate == 1]) + length(object$model$xseed) + 1
    return(object$model$loglik + k * np)
}

#' @method AIC tsissm.selection
#' @rdname AIC
#' @export
#'
#'
AIC.tsissm.selection <- function(object, ..., k = 2)
{
    return(sapply(object$models, function(x) AIC(x, k = k)))
}

#' Bayesian Information Criterion
#'
#' @description Extract the BIC from an estimated model.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @return A numeric value.
#' @aliases BIC
#' @method BIC tsissm.estimate
#' @rdname BIC
#' @export
#'
#'
BIC.tsissm.estimate <- function(object, ...)
{
    nobs_good <- sum(object$spec$good)
    np <- NROW(object$parmatrix[estimate == 1]) + length(object$model$xseed) + 1
    out <- -2 * as.numeric(logLik(object)) + np * log(nobs_good)
    return(out)
}

#' @method BIC tsissm.selection
#' @rdname BIC
#' @export
#'
#'
BIC.tsissm.selection <- function(object, ...)
{
    return(sapply(object$models, function(x) BIC(x)))
}

# tsmetrics ---------------------------------------------------

#' Performance Metrics
#'
#' @description Performance metrics from an estimated or predicted tsissm model.
#' @param object an object of class \dQuote{tsissm.estimate} or \dQuote{tsissm.predict}
#' @param actual the actual data matched to the dates of the forecasts.
#' @param alpha the coverage level for distributional forecast metrics.
#' @param ... not currently used.
#' @returns a data.frame of performance metrics including MAPE, MSLRE, BIAS and AIC
#' for the estimate object and MAPE, MSLRE, BIAS, MASE, MIS and CRPS for predict object.
#' @aliases tsmetrics
#' @method tsmetrics tsissm.predict
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.tsissm.predict = function(object, actual, alpha = 0.1, ...)
{
    n <- NCOL(object$distribution)
    if (NROW(actual) != n) stop("\nactual length not equal to forecast length")
    m_mape <- mape(actual, colMeans(object$distribution))
    m_bias <- bias(actual, colMeans(object$distribution))
    m_mslre <- mslre(actual, colMeans(object$distribution))
    m_mase <- mase(actual, colMeans(object$distribution), object$original_series, frequency = object$frequency[1])
    if (!is.null(alpha)) {
        m_mis <- mis(actual, lower = apply(object$distribution, 2, quantile, alpha/2), upper = apply(object$distribution, 2, quantile, 1 - alpha/2), alpha = alpha)
    } else {
        m_mis <- as.numeric(NA)
    }
    m_crps <- crps(actual, object$distribution)
    return(data.frame("h" = n, "MAPE" = m_mape, "MASE" = m_mase, "MSLRE" = m_mslre, "BIAS" = m_bias, "MIS" = m_mis,"CRPS" = m_crps))
}

#' @method tsmetrics tsissm.estimate
#' @rdname tsmetrics
#' @export
#'
#'
tsmetrics.tsissm.estimate <- function(object, ...)
{
    parameters <- NULL
    estimate <- NULL
    lambda <- object$parmatrix[parameters == "lambda"]$value
    nr <- NROW(object$parmatrix[estimate == 1]) + NCOL(object$model$xseed)
    AIC <- object$model$loglik + 2 * nr
    MAPE <- mape(object$spec$target$y_orig, fitted(object))
    BIAS <- bias(object$spec$target$y_orig, fitted(object))
    MSLRE <- mslre(object$spec$target$y_orig, fitted(object))
    metrics <- data.frame("MAPE" = MAPE, "MSLRE" = MSLRE, "BIAS" = BIAS, "AIC" = AIC)
    return(metrics)
}

# vcov ---------------------------------------------------

#' The Covariance Matrix of the Estimated Parameters
#'
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param adjust logical. Should a finite sample adjustment be made? This amounts
#' to multiplication with n/(n-k) where n is the number of observations and k
#' the number of estimated parameters.
#' @param type valid choices are \dQuote{H} for using the analytic hessian
#' for the bread, \dQuote{OP} for the outer product of gradients, \dQuote{QMLE}
#' for the Quasi-ML sandwich estimator (Huber-White), and \dQuote{NW} for the Newey-West
#' adjusted sandwich estimator (a HAC estimator).
#' @param ... additional parameters passed to the Newey-West bandwidth function to
#' determine the optimal lags.
#' @returns The variance-covariance matrix of the estimated parameters.
#' @method vcov tsissm.estimate
#' @aliases vcov
#' @rdname vcov
#' @export
#'
vcov.tsissm.estimate <- function(object, adjust = FALSE, type = c("H","OP","QMLE","NW"), ...)
{
    estimate <- NULL
    type <- match.arg(type[1],c("H","OP","QMLE","NW"))
    if (type == "H") {
        V <- bread(object)
    } else if (type == "QMLE") {
        bread. <- bread(object)
        meat. <- meat_tsissm(object, adjust = adjust)
        V <- bread. %*% meat. %*% bread.
    } else if (type == "OP") {
        V <- vcovOPG(object, adjust = adjust)
    } else if (type == "NW") {
        bread. <- bread(object)
        meat. <- meatHAC_tsissm(object, adjust = adjust, ...)
        V <- bread. %*% meat. %*% bread.
    }
    par_names <- object$parmatrix[estimate == 1]$parameters
    colnames(V) <- rownames(V) <- par_names
    return(V)
}

# sigma ---------------------------------------------------

#' The Standard Deviation of the model
#'
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... none.
#' @returns If the model was estimated using dynamic variance (GARCH model), then
#' an xts matrix is returned with the estimated GARCH volatility else the estimated
#' standard deviation of the residuals in the case of constant variance.
#' @method sigma tsissm.estimate
#' @aliases sigma
#' @rdname sigma
#' @export
#'
sigma.tsissm.estimate <- function(object, ...) {
    if (object$spec$variance$type == "dynamic") {
        s <- xts(object$model$sigma, object$spec$target$index)
        colnames(s) <- "Sigma"
    } else {
        s <- object$model$sigma
    }
    return(s)
}


# equation ---------------------------------------------------

#' Model Equation (LaTeX)
#'
#' @description Generates a list of model equations in LaTeX.
#' @param object an object of class \dQuote{tsissm.estimate}.
#' @param ... not currently used.
#' @returns A list of equations in LaTeX which can be used in documents. This is
#' a list with 3 slots for the observation, state and variance equations,
#' @details This method is called in the summary when the format output option
#' chosen is \dQuote{flextable}.
#' @aliases tsequation
#' @method tsequation tsissm.estimate
#' @rdname tsequation
#' @export
#'
tsequation.tsissm.estimate <- function(object, ...)
{
    out <- .equation_issm(object$spec)   
    return(out)
}