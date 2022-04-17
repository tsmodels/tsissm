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
