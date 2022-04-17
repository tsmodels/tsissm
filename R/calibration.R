tscalibrate.tsissm.spec <- function(object, start = floor(length(object$target$y_orig))/2, end = length(object$target$y_orig),
                                  h = 1, nsim = 5000, solver = "nlminb", autodiff = TRUE,
                                  autoclean = FALSE, trace = FALSE, ...)
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
    if (is.null(object$target$frequency)) {
        if (is.null(object$seasonal$seasonal_frequency[1])) {
            frequency <- 1
        } else {
            frequency <- object$seasonal$seasonal_frequency[1]
        }
    } else {
        frequency <- object$target$frequency
    }
    start_date <- index(data)[start]
    end <- min(NROW(data), end)
    end_date <- index(data)[end - 1]
    if ((start/2) < max(object$seasonal$seasonal_frequency)) {
        warning(paste0("\nmaximum seasonal frequency > 1/2 initial data length (based on start)"))
    }
    seqdates <- index(data[paste0(start_date,"/", end_date)])
    elapsed_time <- function(idx, end_date, start_date) {
        min(h, which(end_date == idx) - which(start_date == idx))
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
                               transformation = object$transform$name, lower = object$transform$lower,
                               upper = object$transform$upper,
                               ar = object$arma$order[1], ma = object$arma$order[2], 
                               xreg = xreg_train, lambda = lambda, sampling = object$target$sampling)
        mod <- try(estimate(spec, solver = solver, autodiff = autodiff, use_hessian = use_hessian), silent = TRUE)
        if (inherits(mod, 'try-error') | is.null(mod)) {
           return(NULL)
        } else {
            p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(y_test))
            quantiles <- seq(0.1, 1, by = 0.1)
            qp <- apply(p$distribution, 2, quantile, quantiles)
            qp <- t(qp)
            colnames(qp) <- paste0("P", round(quantiles*100))
            if (spec$transform$include_lambda | spec$transform$name == "logit") {
                dist <- do.call(cbind, lapply(1:ncol(p$distribution), function(j) mod$spec$transform$transform(p$distribution[,j], lambda = mod$spec$transform$lambda)))
                ac <- mod$spec$transform$transform(y_test, lambda = mod$spec$transform$lambda)
                er <- unname(as.numeric(ac - colMeans(dist)))
                sigma_d <- sqrt(tsmoments.tsissm.estimate(mod, h = horizon[i])$var)
            } else {
                er <- as.numeric(y_test) - as.numeric(p$mean)
                sigma_d <- sqrt(tsmoments(mod, h = horizon[i])$var)
            }
            out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]),
                              "horizon" = 1:horizon[i],
                              "size" = rep(nrow(y_train), horizon[i]),
                              "forecast_dates" = as.character(index(y_test)),
                              "forecast" = as.numeric(p$mean),
                              "actual" = as.numeric(y_test),
                              "error"  = er,
                              "sigma" = sigma_d)
            out <- cbind(out, qp)
            return(out)
        }
    }, future.packages = c("tsmethods","tsaux","xts","tsissm","data.table"), future.seed = TRUE)
    b <- eval(b)
    b <- rbindlist(b)
    error <- NULL
    sigma <- NULL
    bx <- b[,list(scaling = abs(error)/sigma, error = error), by = c("estimation_date","horizon")]
    bs <- b[,list(rel_scale = sigma[1]/sigma, horizon = horizon), by = c("estimation_date")]
    bx <- merge(bx, bs, by = c("estimation_date","horizon"))
    bxe <- dcast(bx, estimation_date~horizon, value.var = "error")
    bxe <- as.matrix(bxe[,-1])
    bxe <- scale(bxe)
    if (NROW(na.omit(bxe)) < h) {
        V <- cov(bxe, use = "pairwise.complete")
    } else {
        V <- cov(na.omit(bxe))
    }
    V <- make.positive.definite(V)
    E <- eigen(V)
    W <- diag(1/sqrt(E$value)) %*% t(E$vectors)
    Z <- tcrossprod(bxe, W)
    Z <- sweep(Z, 2, colMeans(Z, na.rm = TRUE), "-")
    sigma_scale <- sapply(1:h, function(i) tail(jackknife(bx[horizon == i]$scaling, median)$jack.value, 1))
    sigma_scale_sigma <- sapply(1:h, function(i) tail(jackknife(bx[horizon == i]$scaling, sd)$jack.value, 1))
    samples <- do.call(cbind, lapply(1:h, function(i){
        kmod <- kde(as.numeric(na.omit(Z[,i])))
        s <- rkde(nsim, kmod)
        return(matrix(s, ncol = 1))
    }))
    return(list(samples = samples, sigma_scale = sigma_scale, sd_sigma_scale = sigma_scale_sigma))
}