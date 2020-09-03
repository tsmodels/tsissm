plot.tsissm.estimate <- function(x, y = NULL, ...)
{
    opar <- par()
    opar$cin <- NULL
    opar$cra <- NULL
    opar$csi <- NULL
    opar$cxy <- NULL
    opar$din <- NULL
    opar$page <- NULL
    tsd <- tsdecompose(x)
    # fitted+actual and then states
    a <- x$spec$target$y_orig
    f <- as.numeric(fitted(x))
    dt <- x$spec$target$index
    n <- ncol(tsd) + 1
    colx <- viridis_pal(option = "D", end = 0.8)(n - 1)
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

plot.tsissm.simulate <- function(x, y = NULL, ...)
{
    opar <- par()
    opar$cin <- NULL
    opar$cra <- NULL
    opar$csi <- NULL
    opar$cxy <- NULL
    opar$din <- NULL
    opar$page <- NULL
    n <- length(x$components) + 1
    colx <- viridis_pal(option = "D", end = 0.8)(n)
    par(bg = "white", mar = c(2,2,0.5,3))
    layout(mat = matrix(c(1:n), nrow = n))
    plot(x$Simulated, gradient_color = "whitesmoke")
    mtext("Series", side = 4, adj = 0.5, padj = 0.5, cex = 0.7, font = 2, family = "mono")
    component_names <- x$components
    par(bg = "white",mar = c(0.5,2,0.5,3))
    for (i in 1:(n - 1)) {
        plot(x[[component_names[i]]], gradient_color = colx[i], main = "", ylab = "", xlab = "", x_axes = FALSE, cex.axis = 0.8)
        mtext(component_names[i], side = 4, adj = 0.5, padj = 0.5, cex = 0.7, font = 2, family = "mono")
    }
    suppressWarnings(par(opar))
}


plot.tsissm.profile <- function(x, y = NULL, type = c("metrics", "coef"), ...)
{
    opar <- par()
    opar$cin <- NULL
    opar$cra <- NULL
    opar$csi <- NULL
    opar$cxy <- NULL
    opar$din <- NULL
    opar$page <- NULL
    type <- match.arg(type[1], c("metrics", "coef"))
    if (type == "metrics") {
        layout_matrix <- matrix(1:3, nrow = 3, ncol = 1)
        layout(mat = layout_matrix)
        par(mar = c(2,2,2,4))
        plot(x$MAPE*100, date_class = "integer", interval_quantiles = c(0.1, 0.9), main = "", xlab = "horizon")
        mtext("MAPE[%]", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
        plot(x$BIAS*100, date_class = "integer", interval_quantiles = c(0.1, 0.9), main = "", xlab = "horizon")
        mtext("BIAS[%]", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
        plot(x$MSLRE*100, date_class = "integer", interval_quantiles = c(0.1, 0.9), main = "", xlab = "horizon")
        mtext("MSLRE[%]", side = 4, adj = 0.5, padj = 0.5, cex = 0.8, font = 2, family = "mono")
    } else {
        Variable <- NULL
        cf <- unique(x$coef$Variable)
        n <- length(cf)
        nf <- n2mfrow(n)
        colx <- viridis_pal(alpha = 0.5)(10)
        par(mar = c(2.5,3,2,4), mfrow = nf)
        for (i in 1:n) {
            xlim_plot <- c(min(x$coef[Variable == cf[i]]$Value, x$true.coef[cf[i]]), max(x$coef[Variable == cf[i]]$Value, x$true.coef[cf[i]]))
            xlim_dist <- (xlim_plot[2] - xlim_plot[1])/10
            xlim_plot <- c(xlim_plot[1] - xlim_dist, xlim_plot[2] + xlim_dist)
            hist(x$coef[Variable == cf[i]]$Value, main = cf[i], xlab = "", col = colx[3],  ylab = "", prob = TRUE, freq = FALSE, xaxs = "i",yaxs = "i", xlim = xlim_plot)
            abline(v = x$true.coef[cf[i]], col  = "tomato2", lty = 2, lwd = 2)
            box()
            grid()
        }
    }
    suppressWarnings(par(opar))
}
