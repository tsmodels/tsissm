auto_issm <- function(y, slope = c(TRUE,FALSE), slope_damped = c(TRUE,FALSE), seasonal = c(TRUE,FALSE), 
                      seasonal_frequency = 1, seasonal_type = "trigonometric", seasonal_harmonics = list(),
                      ar = 0:2, ma = 0:2, xreg = NULL, transformation = "box-cox", lambda = 1, 
                      lower = 0, upper = 1, sampling = NULL, cores = 1, trace = FALSE, 
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
    # setup parallel infrastructure
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
    if (trace) {
        iterations <- n
        pb <- txtProgressBar(max = iterations, style = 3)
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
    } else {
        opts <- NULL
    }
    b <- foreach(i = 1:n, .packages = c("tsmethods","tsaux","xts","tsissm","data.table"), .options.snow = opts, .combine = rbind) %dopar% {
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
    }
    stopCluster(cl)
    if (trace) {
        close(pb)
    }
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
    return(mod)
}