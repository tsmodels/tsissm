tsbacktest.tsissm.spec <- function(object, start = floor(length(object$target$y_orig)/2), end = length(object$target$y_orig),
                                  h = 1, alpha = NULL, cores = 1, data_name = "y", save_output = FALSE,
                                  save_dir = "~/tmp/", solver = "optim", trace = FALSE, ...)
{
    if (save_output) {
        if (is.null(save_dir)) {
            stop("save_dir cannot be NULL when save.output is TRUE")
        }
        if (!dir.exists(save_dir)) {
            stop("save_dir does not exit. Create first and then resubmit")
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
    end_date <- index(data)[end - 1]
    if ((start/2) < max(object$seasonal$seasonal_frequency)) {
        warning(paste0("\nmaximum seasonal frequency > 1/2 initial data length (based on start)"))
    }
    seqdates <- index(data[paste0(start_date,"/", end_date)])
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
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (trace) {
        iterations <- length(seqdates)
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    b <- foreach(i = 1:length(seqdates), .packages = c("tsmethods","tsaux","xts","tsissm","data.table"), .options.snow = opts, .combine = rbind) %dopar% {
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
        spec <- issm_modelspec(y_train, slope = object$slope$include_slope, slope_damped = object$slope$include_damped,
                               seasonal = object$seasonal$include_seasonal,
                               seasonal_frequency = object$seasonal$seasonal_frequency,
                               seasonal_type = object$seasonal$seasonal_type,
                               seasonal_harmonics = object$seasonal$seasonal_harmonics,
                               ar = object$arma$order[1], ma = object$arma$order[2], xreg = xreg_train, lambda = lambda, sampling = object$target$sampling)
        mod <- estimate(spec, solver = solver)
        p <- predict(mod, h = horizon[i], newxreg = xreg_test, forc_dates = index(y_test))
        if (save_output) {
            saveRDS(mod, file = paste0(save_dir,"/model_", seqdates[i], ".rds"))
            saveRDS(p, file = paste0(save_dir,"/predict_", seqdates[i], ".rds"))
        }
        if (!is.null(quantiles)) {
            qp <- apply(p$distribution, 2, quantile, quantiles)
            if (length(quantiles) == 1) {
                qp <- matrix(qp, ncol = 1)
            } else{
                qp <- t(qp)
            }
            colnames(qp) <- paste0("P", round(quantiles*100))
        }
        out <- data.table("estimation_date" = rep(seqdates[i], horizon[i]),
                          "horizon" = 1:horizon[i],
                          "size" = rep(nrow(y_train), horizon[i]),
                          "forecast_dates" = as.character(index(y_test)),
                          "forecast" = as.numeric(p$mean), "actual" = as.numeric(y_test))
        if (!is.null(quantiles)) out <- cbind(out, qp)
        return(out)
    }
    stopCluster(cl)
    if (trace) {
        close(pb)
    }
    if (is.null(data_name)) data_name <- "y"
    actual <- NULL
    forecast <- NULL
    metrics <- b[,list(variable = data_name, MAPE = mape(actual, forecast), MSLRE = mslre(actual, forecast),
                       BIAS = bias(actual, forecast),
                       n = .N), by = "horizon"]
    if (!is.null(alpha)) {
        q_names <- matrix(paste0("P", round(quantiles*100)), ncol = 2, byrow = TRUE)
        q <- do.call(cbind, lapply(1:length(alpha), function(i){
            b[,list(mis = mis(actual, get(q_names[i,1]), get(q_names[i,2]), alpha[i])), by = "horizon"]
        }))
        q <- q[,which(grepl("mis",colnames(q))), with = FALSE]
        colnames(q) <- paste0("MIS[",alpha,"]")
        metrics <- cbind(metrics, q)
    }
    return(list(prediction = b, metrics = metrics))
}


