hresiduals_issm <- function(object,  h = 12, nsim = 1000, start = 1, transformed = TRUE, ...)
{
    # truncate h to the available time periods
    N <- length(object$spec$target$y) + 1
    if ((start + h) > N) {
        h <- N - start
    }
    if (transformed) {
        y <- c(0, object$spec$target$y)
        actual <- y[(start + 1):(start + h)]
    } else {
        y <- c(0, object$spec$target$y_orig)
        actual <- y[(start + 1):(start + h)]
    }
    xseed <- object$model$states[start, ,drop = FALSE]
    if (!object$spec$xreg$include_xreg) {
        xreg = matrix(0, ncol = 1, nrow = h)
    } else {
        m <- NCOL(object$spec$xreg$xreg)
        xreg <- rbind(matrix(0, ncol = m, nrow = 1), object$spec$xreg$xreg)[(start + 1):(start + h), ,drop = FALSE]
    }
    forecast_mu <- tsmoments.tsissm.estimate(object, h = h, newxreg = xreg, init_states = xseed, transform = !transformed)$mean
    e <- actual - forecast_mu
    return(data.table(date = object$spec$target$index[start], error = e, forecast = forecast_mu, actual = actual, horizon = 1:h))
}
