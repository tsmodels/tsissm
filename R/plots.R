#' Object Plots
#'
#' @description Plots for objects generated from the tsissm functions.
#' @param x an object of class \dQuote{tsissm.estimate}, \dQuote{tsissm.simulate}
#' or \dQuote{tsissm.profile}.
#' @param y not used.
#' @param type type of profile plot for objects of class \dQuote{tsissm.profile}.
#' @param ... additional arguments passed to the underlying plot function.
#' @returns different plots depending on the input class.
#' @aliases plot
#' @method plot tsissm.estimate
#' @rdname plot
#' @export
#'
#'
plot.tsissm.estimate <- function(x, y = NULL, ...)
{
    opar <- par(mfrow = c(1,1))
    tsd <- tsdecompose(x)
    # fitted+actual and then states
    a <- x$spec$target$y_orig
    f <- as.numeric(fitted(x))
    dt <- x$spec$target$index
    n <- ncol(tsd) + 1
    colx <- .viridis_fun(option = "D", end = 0.8)(n - 1)
    par(bg = "white", mar = c(2,2,0.5,3))
    layout(mat = matrix(c(1:n), nrow = n))
    plot(dt, a, type = "l", main = "", ylab = "", xlab = "", cex.axis = 0.8, col = "black")
    lines(dt, f, col = "red", lty = 2)
    grid()
    mtext("Fitted", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    par(bg = "white",mar = c(0.5,2,0.5,3))
    for (i in 1:ncol(tsd)) {
        plot(dt, as.numeric(tsd[,i]), col = colx[i], type = "l", main = "", ylab = "", xlab = "", xaxt = "n", cex.axis = 0.8)
        mtext(colnames(tsd)[i], side = 4, adj = 0.5, padj = 0.5, cex = 0.7, font = 2, family = "mono")
        grid()
    }
    suppressWarnings(par(opar))
}

#' @method plot tsissm.simulate
#' @rdname plot
#' @export
plot.tsissm.simulate <- function(x, y = NULL, ...)
{
    opar <- par(mfrow = c(1,1))
    components <- tsdecompose(x)
    n <- length(components) + 1
    colx <- (.viridis_fun(option = "H", begin = 0.4, end = 0.9, alpha = 0.5)(n))
    par(bg = "white", mar = c(2,2,0.5,3))
    layout(mat = matrix(c(1:n), nrow = n))
    plot(x$distribution, gradient_color = "azure3", interval_color = "steelblue", median_width = 1, interval_type = 1, interval_width = 1)
    mtext("Series", side = 4, adj = 0.5, padj = 0.5, cex = 0.7, font = 2, family = "mono")
    component_names <- names(components)
    par(bg = "white",mar = c(0.5,2,0.5,3))
    for (i in 1:(n - 1)) {
        plot(components[[component_names[i]]], gradient_color = colx[i], main = "", ylab = "", xlab = "", x_axes = FALSE, cex.axis = 0.8, interval_color = "steelblue", median_width = 1, interval_type = 1, interval_width = 1)
        mtext(component_names[i], side = 4, adj = 0.5, padj = 0.5, cex = 0.7, font = 2, family = "mono")
    }
    suppressWarnings(par(opar))
}

#' @method plot tsissm.profile
#' @rdname plot
#' @export
plot.tsissm.profile <- function(x, y = NULL, type = c("coef","mape","mase","crps"), ...)
{
    Simulation <- Variable <- NULL
    opar <- par(mfrow = c(1,1))
    type <- match.arg(type[1], c("coef","mape","mase","crps"))
    if (type == "coef") {
        true_values <- data.table(Variable = names(x$true_coef),
                                  TrueValue = unname(x$true_coef))
        dt_merged <- merge(x$coef, true_values, by = "Variable", all.x = TRUE)
        cf <- true_values$Variable
        n <- length(cf)
        nf <- n2mfrow(n)
        colx <- .viridis_fun(alpha = 0.5, begin = 0, end = 0.7, option = "B")(n)
        par(mar = c(2.5,3,2,4), mfrow = nf)
        for (i in 1:n) {
            plot(density(dt_merged[Variable == cf[i]]$Value, adjust = 1.5), col = colx[i], xlab = "", ylab = "", main = cf[i])
            grid()
            abline(v = dt_merged[Variable == cf[i]]$TrueValue[1], col  = "steelblue", lty = 2, lwd = 1)
        }
    } else if (type == "mape") {
        tmp <- dcast(x$MAPE, Simulation~Horizon, value.var = "MAPE")
        tmp[,Simulation := NULL]
        colx <- rev(.viridis_fun(alpha = 0.5, begin = 0.2, end = 0.8, option = "H")(ncol(tmp)))
        boxplot(round(tmp * 100, 2), xlab = "Horizon", ylab = "MAPE (%)", col = colx, outline = FALSE, main = "MAPE by Horizon")
    } else if (type == "mase") {
        tmp <- dcast(x$MASE, Simulation~Horizon, value.var = "MASE")
        tmp[,Simulation := NULL]
        colx <- rev(.viridis_fun(alpha = 0.5, begin = 0, end = 0.7, option = "H")(ncol(tmp)))
        boxplot(round(tmp, 2), xlab = "Horizon", ylab = "MASE", col = colx, outline = FALSE, main = "MASE by Horizon")
        abline(h = 1, lty = 2, col = "grey")
    } else if (type == "crps") {
        tmp <- dcast(x$CRPS, Simulation~Horizon, value.var = "CRPS")
        tmp[,Simulation := NULL]
        colx <- rev(.viridis_fun(alpha = 0.5, begin = 0, end = 0.7, option = "H")(ncol(tmp)))
        boxplot(round(tmp, 2), xlab = "Horizon", ylab = "CRPS", col = colx, outline = FALSE, main = "CRPS by Horizon")
    }
    suppressWarnings(par(opar))
    return(invisible(x))
}
