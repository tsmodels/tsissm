tsfilter.tsissm.estimate <- function(object, y = NULL, newxreg = NULL, ...)
{
    parameters <- NULL
    yold <- xts(object$spec$target$y_orig, object$spec$target$index)
    ynew <- y
    exc <- which(index(ynew) %in% index(yold))
    if (length(exc) == 0) {
        y <- ynew
    } else{
        y <- ynew[-exc]
        if (NROW(y) == 0) {
            warning("\nno new data in y which is not already in object!")
            return(object)
        }
    }
    if (object$spec$xreg$include_xreg) {
        nx <- NCOL(object$spec$xreg$xreg)
        if (!is.null(newxreg)) {
            if (ncol(newxreg) != nx) stop(paste0("\nExpected ", nx, " columns in newxreg but got ", NCOL(newxreg)))
            newcnames <- colnames(newxreg)
            oldcnames <- colnames(object$spec$xreg$xreg)
            if (!is.null(newcnames) & !is.null(oldcnames)) {
                if (!all(sort(oldcnames) %in% sort(newcnames))) {
                    stop("\ncolnames of newxreg do not match those of original xreg")
                }
                newxreg <- newxreg[, oldcnames]
            }
            X <- coredata(newxreg)
            if (length(exc) > 0) {
                X <- X[-exc,,drop = FALSE]
            }
        } else {
            X <- matrix(0, ncol = nx, nrow = nrow(y))
        }
    } else {
        X <- matrix(0, ncol = 1, nrow = NROW(y))
    }
    newindex <- index(y)
    yneworig <- y
    xseed <- tail(object$model$states, 1)
    pars <- object$parmatrix$optimal
    names(pars) <- object$parmatrix$parameters
    Mnames <- na.omit(object$spec$S$pars)
    S <- object$spec$S
    S[which(!is.na(pars)),"values"] <- pars[Mnames]
    ##################################################################
    parnames <- names(pars)
    mdim = object$spec$dims
    mdim[2] <- NROW(y)
    f <- issfilter(f0_ = S[list("F0")]$values,
                      f1_ = S[list("F1")]$values,
                      f2_ = S[list("F2")]$values,
                      w_ = as.numeric(S[list("w")]$values),
                      g_ = as.numeric(S[list("g")]$values),
                      y_ = as.numeric(y),
                      lambda_ = as.numeric(pars["lambda"]),
                      xreg_ = as.vector(X),
                      kappa_ = S[list("kappa")]$values,
                      mdim = mdim, xseed_ = as.vector(xseed))
    # update all vectors with the y
    object$spec$target$y_orig <- c(object$spec$target$y_orig, as.numeric(yneworig))
    object$spec$target$index <- c(object$spec$target$index, newindex)
    object$spec$target$y <- c(object$spec$target$y, as.numeric(y))
    mdim[2] <- length(object$spec$target$y_orig)
    if (!is.null(newxreg)) {
        object$spec$xreg$xreg <- rbind(object$spec$xreg$xreg, X)
    } else {
        object$spec$xreg$xreg <- cbind(object$spec$xreg$xreg, matrix(0, ncol = ncol(object$spec$xreg$xreg), nrow = NROW(yneworig)))
    }
    if (object$spec$transform$lambda != 1) {
        filt <- object$spec$transform$inverse(f$fitted, lambda = as.numeric(pars["lambda"]))
    } else {
        filt <- f$fitted
    }
    filt <- filt[-1]
    err <- as.numeric(y) - filt
    object$model$fitted <- c(object$model$fitted, filt)
    object$model$states <- rbind(object$model$states, f$states[-1,,drop = FALSE])
    object$model$residuals <- c(object$model$residuals, err)
    return(object)
}