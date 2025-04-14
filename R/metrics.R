accuracy <- function(actual, predicted)
{
    actual <- as.numeric(actual)
    predicted <- as.numeric(predicted)
    n <- length(actual)
    if (length(predicted) != n) stop("\nactual and predicted must be of the same length.")
    b <- cbind(actual, predicted)
    b <- na.omit(b)
    n <- NROW(b)
    if (n == 0) return(NA)
    s <- sign(b[,1] * b[,2])
    a <- which(s == 1)
    metric <- (length(a)/n) * 100
    return(metric)
}


.rosner_lambda <- function(n, k = 10, alpha = 0.05) 
{
    aruments_matrix <- cbind(n = as.vector(n), k = as.vector(k), alpha = as.vector(alpha))
    for (i in c("n", "k", "alpha")) assign(i, aruments_matrix[, i])
    l <- k - 1
    p <- 1 - ((alpha/2)/(n - l))
    t_critical <- qt(p = p, df = n - l - 2)
    lambda <- t_critical * (n - l - 1)/sqrt((n - l - 2 + t_critical^2) * (n - l))
    names(lambda) <- NULL
    return(lambda)
}


.rosner_test <- function(x, k = 3, alpha = 0.05, warn = TRUE)
{
    data.name <- deparse(substitute(x))
    if (!is.numeric(x)) stop("'x' must be a numeric vector")
    obs.num <- 1:length(x)
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        x <- x[x.ok]
        obs.num <- obs.num[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    if (n < 3) stop("There must be at least 3 non-missing finite observations in 'x'")
    if (length(k) != 1 || !is.numeric(k) || !is.finite(k) || 
        k != round(k) || k < 1 || k > (n - 2)) 
        stop(paste("'k' must be a positive integer less than or equal to n-2,", 
                   "where 'n' denotes the number of finite, non-missing observations in 'x'"))
    if (length(alpha) != 1 || !is.numeric(alpha) || !is.finite(alpha) || 
        any(alpha <= 0) || any(alpha >= 1)) 
        stop("'alpha' must be a numeric scalar greater than 0 and less than 1")
    if (warn) {
        if (k > 10 | k > floor(n/2)) {
            warning(paste("The true Type I error may be larger than assumed.", 
                          "Although the help file for 'rosnerTest' has a table with information", 
                          "on the estimated Type I error level,", "simulations were not run for k > 10 or k > floor(n/2).", 
                          sep = "\n"))
        }
        else {
            warn.conds <- (alpha > 0.01 & ((n >= 15 & n < 25 & k > 2) | (n < 15 & k > 1))) | (alpha <= 0.01 & (n < 15 & k > 1))
            if (warn.conds) 
                warning(paste("The true Type I error may be larger than assumed.", 
                              "See the help file for 'rosnerTest' for a table with information", 
                              "on the estimated Type I error level."))
        }
    }
    R <- rep(as.numeric(NA), k)
    mean.vec <- rep(as.numeric(NA), k)
    sd.vec <- rep(as.numeric(NA), k)
    x.vec <- rep(as.numeric(NA), k)
    obs.num.vec <- rep(as.numeric(NA), k)
    new.x <- x
    new.obs.num <- obs.num
    for (i in 1:k) {
        mean.vec[i] <- mean(new.x)
        sd.vec[i] <- sd(new.x)
        if (sd.vec[i] > 0) {
            abs.z = abs(new.x - mean.vec[i])/sd.vec[i]
            R[i] <- max(abs.z)
            index <- which(abs.z == R[i])[1]
            x.vec[i] <- new.x[index]
            obs.num.vec[i] <- new.obs.num[index]
            new.x <- new.x[-index]
            new.obs.num <- new.obs.num[-index]
        }
        else {
            R[i:k] <- NA
            break
        }
    }
    num.outlier.vec <- 1:k
    lambda <- .rosner_lambda(n = n, k = 1:k, alpha = alpha)
    outlier <- R > lambda
    if (any(outlier)) {
        index <- max(num.outlier.vec[outlier], na.rm = TRUE)
        outlier[1:index] <- TRUE
    }
    stats_table <- data.frame(num.outlier.vec - 1, mean.vec, sd.vec, x.vec, obs.num.vec, R, lambda, outlier)
    names(stats_table) <- c("i", "Mean.i", "SD.i", "Value", "Obs.Num", "R.i+1", "lambda.i+1", "Outlier")
    return(stats_table)
}