ss_matrix_A <- function(frequency, K, type = "trigonometric")
{
    if (type == "regular") {
        n <- length(frequency)
        x <- lapply(1L:n, function(i){
            rbind(
                cbind(matrix(0, 1L, frequency[i] - 1L),
                      matrix(1, 1L, 1L)),
                cbind(diag(1,   frequency[i] - 1L, frequency[i] - 1L),
                      matrix(0, frequency[i] - 1L, 1L))
        )})
    } else {
        n <- length(frequency)
        r <- NULL
        r <- lapply(1L:n, function(i) 2 * pi * (1L:K[i]) / frequency[i])
        x <- lapply(1L:n, function(i){
            rbind(
                cbind(diag( cos(r[[i]]), K[i], K[i]),
                      diag( sin(r[[i]]), K[i], K[i])),
                cbind(diag(-sin(r[[i]]), K[i], K[i]),
                      diag( cos(r[[i]]), K[i], K[i])))
        })
    }
    A <- bdiag(x)
    return(A)
}

ss_vector_gamma <- function(gamma, frequency, type = "trigonometric", K)
{
    if (type == "trigonometric") {
        k <- rep(K, 1, each = 2L)
        gamma <- rep(gamma, 1, each = 2L)
        gamma_name <- rep(paste0("gamma", frequency), 1, each = 2L)
        gamma_name[seq(1L, length(gamma_name), by = 2)] <- paste0(gamma_name[seq(1L, length(gamma_name), by = 2)],".1")
        gamma_name[seq(2L, length(gamma_name), by = 2)] <- paste0(gamma_name[seq(2L, length(gamma_name), by = 2)],".2")
        g <-  c(unlist(sapply(1L:length(k), function(i) rep(gamma[i], k[i]))))
        gnames <- c(unlist(sapply(1L:length(k), function(i) rep(gamma_name[i], k[i]))))
        names(g) <- gnames
    } else {
        g <- c(unlist(sapply(1L:length(frequency), function(i) c(gamma[i], rep(0, frequency[i] - 1)))))
        names(g)[which(is.na(g))] <- paste0("gamma",frequency)
    }
    return(g)
}

ss_vector_a <- function(frequency, type = "trigonometric", K)
{
    if (type == "trigonometric") {
        a <- matrix(unlist(sapply(1L:length(K), function(i) c(rep(1, K[i]), rep(0, K[i])))), nrow = 1)
    } else {
        a <- matrix(unlist(sapply(1L:length(frequency), function(i) c(rep(0, frequency[i] - 1), 1))), nrow = 1)
    }
    return(a)
}

ss_matrix_x <- function(N, slope = TRUE, frequency, type = "trigonometric", K, ar = 0, ma = 0)
{
    if (frequency[1] > 1) {
        if (type == "trigonometric") {
            s <- matrix(unlist(sapply(1L:length(K), function(i) rep(0, 2*K[i]))), nrow = 1L)
        } else {
            s <- matrix(unlist(sapply(1L:length(frequency), function(i) rep(0, frequency[i] - 1))), nrow = 1L)
        }
    } else {
        s <- NULL
    }
    if (ar > 0) {
        p <- matrix(rep(0, ar), nrow = 1L)
    } else {
        ar <- NULL
    }
    if (ma > 0) {
        q <- matrix(rep(0 ,ma), nrow = 1L)
    } else {
        q <- NULL
    }
    if (slope) {
        b <- matrix(0,1,1)
    } else {
        b <- NULL
    }
    x <- t(replicate(n = N + 1L, cbind(matrix(0,1,1), b, s, p, q), simplify = "matrix"))
    return(x)
}

ss_state_names <- function(level = 1, slope = 1, frequency, type = "trigonometric", K, ar = 0, ma = 0)
{
    ss_names <- c()
    if (level == 1) {
        ss_names <- c(ss_names, "Level")
    }
    if (slope == 1) {
        ss_names <- c(ss_names, "Slope")
    }
    if (frequency[1L] > 1) {
        if (type == "trigonometric") {
            k <- rep(K, 1, each = 2L)
            frequency <- rep(frequency, 1, each = 2L)
            frequency[seq(2L,length(frequency), by = 2L)] <- paste0(frequency[seq(2L,length(frequency), by = 2L)],".")
            tmp <- c(unlist(sapply(1L:length(frequency), function(i) paste0("S",frequency[i],".",1L:k[i]))))
            ss_names <- c(ss_names, tmp)
        } else {
            tmp <- c(unlist(sapply(1L:length(frequency), function(i) paste0("S",frequency[i],".",1L:(frequency[i] - 1L)))))
            ss_names <- c(ss_names, tmp)
        }
    }
    if (ar > 0) {
        ss_names <- c(ss_names, paste0("ar", 1L:ar))
    }
    if (ma > 0) {
        ss_names <- c(ss_names, paste0("ma", 1L:ma))
    }
    return(ss_names)
}


ss_seasonal_dimensions <- function(frequency, type, K)
{
    if (frequency[1] > 1) {
        if (type == "trigonometric") {
            n <- 2*sum(K)
        } else {
            n <- sum(frequency)
        }
    } else {
        n <- 0
    }
    return(n)
}



ss_level <- function(slope = TRUE, frequency = 1, type = "trigonometric", K = NULL, ar = 0, ma = 0)
{
    if (frequency[1] > 1) {
        tau <- ifelse(type == "trigonometric", 2 * sum(K), sum(frequency))
    } else {
        tau <- 0
    }
    f0 <- matrix(1, 1L, 1L)
    if (slope) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = 1L))
    if (frequency[1] > 0) f0 <- rbind(f0, matrix(0, nrow = tau, ncol = 1))
    if (ar > 0) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = 1L))
    if (ar > 1) f0 <- rbind(f0, matrix(0, nrow = ar - 1, ncol = 1))
    if (ma > 0) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = 1L))
    if (ma > 1) f0 <- rbind(f0, matrix(0, nrow = ma - 1, ncol = 1))
    w <- c(1)
    w_names <- NULL
    g <- c(NA)
    g_names <- c("alpha")
    f1 <- matrix(1, nrow = nrow(f0), ncol = ncol(f0))
    f2 <- f1
    f1_names <- NULL
    f2_names <- NULL
    return(list(f0 = f0, f1 = f1, f2 = f2, f1_names = f1_names, f2_names = f2_names, g = g, g_names = g_names, w = w, w_names = w_names))
}

ss_slope <- function(slope = TRUE, damped_slope = FALSE, frequency = 1, type = "trigonometric", K = NULL, ar = 0, ma = 0)
{
    if (frequency[1] > 1) {
        tau <- ifelse(type == "trigonometric", 2 * sum(K), sum(frequency))
    } else {
        tau <- 0
    }
    f0 <- matrix(1, 1L, 1L)
    f0 <- rbind(f0, matrix(1, nrow = 1L, ncol = 1L))
    if (frequency[1] > 0) f0 <- rbind(f0, matrix(0, nrow = tau, ncol = 1))
    if (ar > 0) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = 1L))
    if (ar > 1) f0 <- rbind(f0, matrix(0, nrow = ar - 1, ncol = 1))
    if (ma > 0) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = 1L))
    if (ma > 1) f0 <- rbind(f0, matrix(0, nrow = ma - 1, ncol = 1))
    if (damped_slope){
        phi <- NA
        w_names <- "phi"
    } else {
        phi <- 1
        w_names <- NULL
    }
    w <- c(phi)
    g <- c(NA)
    g_names <- c("beta")
    f1 <- matrix(1, nrow = nrow(f0), ncol = ncol(f0))
    f2 <- f1
    f2[1:2,1] <- phi
    if (damped_slope) f2_names <- c("phi","phi") else f2_names = NULL
    f1_names <- NULL
    return(list(f0 = f0, f1 = f1, f2 = f2, f1_names = f1_names, f2_names = f2_names, g = g, g_names = g_names, w = w, w_names = w_names))
}

ss_seasonal <- function(slope = TRUE, frequency = 1, type = "trigonometric", K = NULL, ar = 0, ma = 0)
{
    if (frequency[1] > 1) {
        A <- ss_matrix_A(frequency = frequency, type = type, K = K)
        tau <- ifelse(type == "trigonometric", 2 * sum(K), sum(frequency))
        a <- ss_vector_a(frequency, type = type, K)
    } else {
        a <- NULL
        tau <- 0
        A <- NULL
    }
    f0 <- matrix(0, nrow = 1L, ncol = tau)
    if (slope) f0 <- rbind(f0, f0)
    f0 <- rbind(f0, A)
    if (ar > 0) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = tau))
    if (ar > 1) f0 <- rbind(f0, matrix(0, nrow = ar - 1, ncol = tau))
    if (ma > 0) f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = tau))
    if (ma > 1) f0 <- rbind(f0, matrix(0, nrow = ma - 1, ncol = tau))
    f2 <- f1 <- matrix(1, nrow = nrow(f0), ncol = ncol(f0))
    f1_names = f2_names = c()
    w <- a
    w_names <- NULL
    g <- ss_vector_gamma(rep(NA, length(frequency)), frequency, type = type, K)
    g_names <- na.omit(names(g))
    return(list(f0 = f0, f1 = f1, f2 = f2, f1_names = f1_names, f2_names = f2_names, g = g, g_names = g_names, w = w, w_names = w_names))
}

ss_ar <- function(slope = TRUE, frequency = 1, type = "trigonometric", K = NULL, ar = 0, ma = 0)
{
    if (frequency[1] > 1) {
        tau <- ifelse(type == "trigonometric", 2 * sum(K), sum(frequency))
    } else {
        tau <- 0
    }
    f0 <- matrix(1, nrow = 1L, ncol = ar)
    f1 <- matrix(NA, nrow = 1L, ncol = ar)
    f2 <- matrix(NA, nrow = 1L, ncol = ar)
    if (slope) {
        f0 <- rbind(f0, f0)
        f1 <- rbind(f1, f1)
        f2 <- rbind(f2, f2)
    }
    if (frequency[1L] > 1) {
        gamma <- ss_vector_gamma(rep(NA, length(frequency)), frequency, type = type, K)
        f0 <- rbind(f0, matrix(1, ncol = ar, nrow = tau))
        if (type == "trigonometric") {
            f1 <- rbind(f1, matrix(NA, ncol = ar, nrow = tau))
            f2 <- rbind(f2, matrix(NA, ncol = ar, nrow = tau))
        } else {
            tmp <- gamma %*% t(c(1L:ar))
            f1 <- rbind(f1, tmp)
            f2 <- rbind(f2, tmp)
        }
    }
    f0 <- rbind(f0, matrix(1, ncol = ar, nrow = 1L))
    f1 <- rbind(f1, matrix(1, ncol = ar, nrow = 1L))
    f2 <- rbind(f2, matrix(NA, ncol = ar, nrow = 1L))

    if (ar > 1) {
        f0 <- rbind(f0, diag(1L, ar - 1L, ar))
        f1 <- rbind(f1, matrix(1, nrow = ar - 1L, ncol = ar))
        f2 <- rbind(f2, matrix(1, nrow = ar - 1L, ncol = ar))
    }
    if (ma > 0) {
        f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = ar))
        f1 <- rbind(f1, matrix(1, nrow = 1L, ncol = ar))
        f2 <- rbind(f2, matrix(1, nrow = 1L, ncol = ar))
    }
    if (ma > 1) {
        f0 <- rbind(f0, matrix(0, nrow = ma - 1L, ncol = ar))
        f1 <- rbind(f1, matrix(1, nrow = ma - 1L, ncol = ar))
        f2 <- rbind(f2, matrix(1, nrow = ma - 1L, ncol = ar))
    }
    f1_names = f2_names = c()
    for (i in 1:ncol(f1)) {
        f1_names <- c(f1_names, "alpha")
        f2_names <- c(f2_names, c(paste0("theta",i)))
        if (slope) {
            f1_names <- c(f1_names, "beta")
            f2_names <- c(f2_names, c(paste0("theta",i)))
        }
        if (frequency[1] > 1) {
            if (type == "trigonometric") {
                f1_names <- c(f1_names, names(gamma))
                f2_names <- c(f2_names, rep(paste0("theta",i), tau))
            } else {
                f1_names <- c(f1_names, na.omit(names(gamma)))
                f2_names <- c(f2_names, rep(paste0("theta",i), length(frequency)))
            }
        }
        f2_names <- c(f2_names, paste0("theta",i))
    }
    w <- rep(NA, ar)
    w_names <- paste0("theta",1:ar)
    g <- c(1)
    if (ar > 1) g <- c(g, rep(0, ar - 1))
    g_names <- NULL
    return(list(f0 = f0, f1 = f1, f2 = f2, f1_names = f1_names, f2_names = f2_names, g = g, g_names = g_names, w = w, w_names = w_names))
}

ss_ma <- function(slope = TRUE, frequency = 1, type = "trigonometric", K = NULL, ar = 0, ma = 0)
{
    if (frequency[1] > 1) {
        tau <- ifelse(type == "trigonometric", 2 * sum(K), sum(frequency))
    } else {
        tau <- 0
    }
    f0 <- matrix(1, nrow = 1L, ncol = ma)
    f1 <- matrix(NA, nrow = 1L, ncol = ma)
    f2 <- matrix(NA, nrow = 1L, ncol = ma)
    if (slope) {
        f0 <- rbind(f0, f0)
        f1 <- rbind(f1, f1)
        f2 <- rbind(f2, f2)
    }
    if (frequency[1L] > 1) {
        gamma <- ss_vector_gamma(rep(NA, length(frequency)), frequency = frequency, type = type, K = K)
        f0 <- rbind(f0, matrix(1, ncol = ma, nrow = tau))
        if (type == "trigonometric") {
            f1 <- rbind(f1, matrix(NA, ncol = ma, nrow = tau))
            f2 <- rbind(f2, matrix(NA, ncol = ma, nrow = tau))
        } else {
            tmp <- gamma %*% t(c(1L:ma))
            f1 <- rbind(f1, tmp)
            f2 <- rbind(f2, tmp)
        }
    }
    if (ar > 0) {
        f0 <- rbind(f0, matrix(1, ncol = ma, nrow = 1L))
        f1 <- rbind(f1, matrix(1, ncol = ma, nrow = 1L))
        f2 <- rbind(f2, matrix(NA, ncol = ma, nrow = 1L))
    }
    if (ar > 1) {
        f0 <- rbind(f0, diag(0, ar - 1L, ma))
        f1 <- rbind(f1, matrix(1, nrow = ar - 1L, ncol = ma))
        f2 <- rbind(f2, matrix(1, nrow = ar - 1L, ncol = ma))
    }
    f0 <- rbind(f0, matrix(0, nrow = 1L, ncol = ma))
    f1 <- rbind(f1, matrix(1, nrow = 1L, ncol = ma))
    f2 <- rbind(f2, matrix(1, nrow = 1L, ncol = ma))

    if (ma > 1) {
        f0 <- rbind(f0, diag(1, ma - 1L, ma))
        f1 <- rbind(f1, matrix(1, nrow = ma - 1L, ncol = ma))
        f2 <- rbind(f2, matrix(1, nrow = ma - 1L, ncol = ma))
    }
    f1_names = f2_names = c()
    for (i in 1:ncol(f1)) {
        f1_names <- c(f1_names, "alpha")
        f2_names <- c(f2_names, c(paste0("psi",i)))
        if (slope) {
            f1_names <- c(f1_names, "beta")
            f2_names <- c(f2_names, c(paste0("psi",i)))
        }
        if (frequency[1] > 1) {
            if (type == "trigonometric") {
                f1_names <- c(f1_names, names(gamma))
                f2_names <- c(f2_names, rep(paste0("psi",i), tau))
            } else {
                f1_names <- c(f1_names, na.omit(names(gamma)))
                f2_names <- c(f2_names, rep(paste0("psi",i), length(frequency)))
            }
        }
        if (ar > 0 ) f2_names <- c(f2_names, paste0("psi",i))
    }
    w <- rep(NA, ma)
    w_names <- paste0("psi",1:ma)
    g <- c(1)
    if (ma > 1) g <- c(g, rep(0, ma - 1))
    g_names <- NULL
    return(list(f0 = f0, f1 = f1, f2 = f2, f1_names = f1_names, f2_names = f2_names, g = g, g_names = g_names, w = w, w_names = w_names))
}


ss_matrices <- function(y, slope = TRUE, damped_slope = FALSE, frequency = 1, type = "trigonometric", K = NULL, ar = 0, ma = 0, include_xreg = FALSE, xreg = NULL)
{
    if (is.null(frequency)) frequency <- 1
    M <- ss_level(slope = slope, frequency = frequency, type = type, K = K, ar = ar, ma = ma)
    F0 <- M$f0
    F1 <- M$f1
    F2 <- M$f2
    F1_names <- M$f1_names
    F2_names <- M$f2_names
    g <- M$g
    g_names <- M$g_names
    w <- M$w
    w_names <- M$w_names
    pars <- c("alpha")
    lower_bounds <- c(1e-12)
    upper_bounds <- c(0.99)
    if (slope) {
        tmp <- ss_slope(slope = slope, damped_slope = damped_slope, frequency = frequency, type = type, K = K, ar = ar, ma = ma)
        pars <- c(pars, "beta")
        lower_bounds <- c(lower_bounds, 1e-12)
        upper_bounds <- c(upper_bounds, 1 - 1e-2)
        if (damped_slope) {
            pars <- c(pars, "phi")
            lower_bounds <- c(lower_bounds, 0.5)
            upper_bounds <- c(upper_bounds, 1)
        }
        F0 <- cbind(F0, tmp$f0)
        F1 <- cbind(F1, tmp$f1)
        F2 <- cbind(F2, tmp$f2)
        F1_names <- c(F1_names, tmp$f1_names)
        F2_names <- c(F2_names, tmp$f2_names)
        g <- c(g, tmp$g)
        g_names <- c(g_names, tmp$g_names)
        w <- c(w, tmp$w)
        w_names <- c(w_names, tmp$w_names)
    }
    if (frequency[1] > 1) {
        tmp <- ss_seasonal(slope = slope, frequency = frequency, type = type, K = K, ar = ar, ma = ma)
        pars <- c(pars, unique(tmp$g_names))
        lower_bounds <- c(lower_bounds, rep( -0.01, length(unique(tmp$g_names))))
        upper_bounds <- c(upper_bounds, rep( 0.99, length(unique(tmp$g_names))))
        F0 <- cbind(F0, tmp$f0)
        F1 <- cbind(F1, tmp$f1)
        F2 <- cbind(F2, tmp$f2)
        F1_names <- c(F1_names, tmp$f1_names)
        F2_names <- c(F2_names, tmp$f2_names)
        g <- c(g, tmp$g)
        g_names <- c(g_names, tmp$g_names)
        w <- c(w, tmp$w)
        w_names <- c(w_names, tmp$w_names)
    }
    if (ar > 0) {
        tmp <- ss_ar(slope = slope, frequency = frequency, type = type, K = K, ar = ar, ma = ma)
        pars <- c(pars, paste0("theta",1L:ar))
        lower_bounds <- c(lower_bounds, rep(-1 + 1e-2, ar))
        upper_bounds <- c(upper_bounds, rep( 1 - 1e-2, ar))
        F0 <- cbind(F0, tmp$f0)
        F1 <- cbind(F1, tmp$f1)
        F2 <- cbind(F2, tmp$f2)
        F1_names <- c(F1_names, tmp$f1_names)
        F2_names <- c(F2_names, tmp$f2_names)
        g <- c(g, tmp$g)
        g_names <- c(g_names, tmp$g_names)
        w <- c(w, tmp$w)
        w_names <- c(w_names, tmp$w_names)
    }
    if (ma > 0) {
        tmp <- ss_ma(slope = slope, frequency = frequency, type = type, K = K, ar = ar, ma = ma)
        pars <- c(pars, paste0("psi",1L:ma))
        lower_bounds <- c(lower_bounds, rep(-1 + 1e-2, ma))
        upper_bounds <- c(upper_bounds, rep( 1 - 1e-2, ma))
        F0 <- cbind(F0, tmp$f0)
        F1 <- cbind(F1, tmp$f1)
        F2 <- cbind(F2, tmp$f2)
        F1_names <- c(F1_names, tmp$f1_names)
        F2_names <- c(F2_names, tmp$f2_names)
        g <- c(g, tmp$g)
        g_names <- c(g_names, tmp$g_names)
        w <- c(w, tmp$w)
        w_names <- c(w_names, tmp$w_names)
    }

    p1 <- rep(NA,length(F1))
    p1[which(is.na(F1))] <- F1_names
    p2 <- rep(NA,length(F2))
    p2[which(is.na(F2))] <- F2_names
    p3 <- rep(NA, length(w))
    p3[which(is.na(w))] <- w_names
    p4 <- rep(NA, length(g))
    p4[which(is.na(g))] <- g_names
    # Level, Slope, Seasonal, AR, MA
    Smatrix <- matrix(0, ncol = 3, nrow = 2 + length(frequency) + 2 + 1)
    colnames(Smatrix) <- c("n","start","end")
    rownames(Smatrix) <- c("Level","Slope", paste0("Seasonal",frequency), paste0("AR"),paste0("MA"),"X")
    k <- 1
    Smatrix[k, 1:2] <- 1
    Smatrix[k, 3]  <- 1
    k <- k + 1
    if (slope) {
        Smatrix[k, 1] <- 1
        Smatrix[k, 2] <- max(Smatrix[1:(k - 1), 3]) + 1
        Smatrix[k, 3] <- max(Smatrix[1:(k - 1), 3]) + 1
    }
    if (frequency[1] > 1) {
        for (i in 1:length(frequency)) {
            k <- k + 1
            if (type == "trigonometric") {
                Smatrix[k, 1] <- K[i]
                Smatrix[k, 2] <- max(Smatrix[1:(k - 1), 3]) + 1
                Smatrix[k, 3] <- max(Smatrix[1:(k - 1), 3]) + 2 * K[i]
            } else {
                Smatrix[k, 1] <- frequency[i]
                Smatrix[k, 2] <- max(Smatrix[1:(k - 1), 3]) + 1
                Smatrix[k, 3] <- max(Smatrix[1:(k - 1), 3]) + frequency[i]
            }
        }
    } else {
        k <- k + 1
    }
    k <- k + 1
    if (ar > 0) {
        Smatrix[k, 1] <- ar
        Smatrix[k, 2] <- max(Smatrix[1:(k - 1), 3]) + 1
        Smatrix[k, 3] <- max(Smatrix[1:(k - 1), 3]) + ar
    }
    k <- k + 1
    if (ma > 0) {
        Smatrix[k, 1] <- ma
        Smatrix[k, 2] <- max(Smatrix[1:(k - 1), 3]) + 1
        Smatrix[k, 3] <- max(Smatrix[1:(k - 1), 3]) + ma
    }
    k <- k + 1
    Smatrix[k ,1] <- ifelse(include_xreg, NCOL(xreg), 0)
    Smatrix[k ,2] <- 1
    Smatrix[k ,3] <- NCOL(xreg)
    D <- c("m" = NROW(F0), "n" = length(y), "s" = ifelse(frequency[1] > 1, ifelse(type == "trigonometric", 1, 2), 0))
    if (frequency[1] > 1) {
        seasonal_start <- Smatrix[which(grepl("Seasonal",rownames(Smatrix)))[1],2]
        seasonal_end <- Smatrix[max(which(grepl("Seasonal",rownames(Smatrix)))),3]
    } else {
        seasonal_start = 1
        seasonal_end = max(Smatrix[1:2,3])
    }
    arma_n <- ar + ma
    D <- c(D, seasonal_start, seasonal_end, arma_n, NCOL(xreg))
    setup <- rbind(
        data.table(matrix = "F0", values = as.numeric(F0), pars = NA),
        data.table(matrix = "F1", values = as.numeric(F1), pars = p1),
        data.table(matrix = "F2", values = as.numeric(F2), pars = p2),
        data.table(matrix = "w", values = w, pars = p3),
        data.table(matrix = "g", values = g, pars = p4))
    # create dims for seasonal type and indices of parameters
    L <- list(setup = setup, dims = D, idmatrix = Smatrix, pars = pars, lower = lower_bounds, upper = upper_bounds)
    return(L)
}

bdiag <- function(x = NULL, ...)
{
    if (is.null(x)) {
        if (nargs() == 1L) x <- as.list(...) else x <- list(...)
    }
    n <- length(x)
    if (n == 0L) return(NULL)
    x <- lapply(x, function(y) if (length(y)) as.matrix(y) else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1L, ]
    cc <- d[2L, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1L, -1L] <- rcum[-n]
    ind[2L, ] <- rcum
    ind[3L, -1L] <- ccum[-n]
    ind[4L, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1L] + 1):y[2L], (y[3L] + 1):y[4L]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}
