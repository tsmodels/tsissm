ar_jconstraint <- function(pars, fun, issmenv)
{
    arconstraint <- function(parms) {
        # Overload assignment for AD-aware indexing
        "[<-" <- ADoverload("[<-")
        pars <- parms$issmpars[issmenv$arindex]
        n <- length(pars)
        # Construct the companion matrix.
        # For a polynomial 1 - a1*x - ... - an*x^n = 0, I want the monic form:
        #   x^n + (a_{n-1}/a_n)*x^(n-1) + ... + (a_1/a_n)*x - 1/a_n = 0.
        # The standard companion matrix for
        #   x^n + c_{n-1}*x^(n-1) + ... + c_0 = 0
        # is:
        #   [ -c_{n-1}   -c_{n-2}  ...  -c_0 ]
        #   [     1          0     ...   0  ]
        #   [     0          1     ...   0  ]
        #   [    ...       ...    ...  ... ]
        #
        # Here, c_{n-1} = a_{n-1}/a_n, c_{n-2} = a_{n-2}/a_n, ..., c_0 = -1/a_n.
        # Because the input vector is ordered as (a1, a2, ..., an), need to reverse
        # the first n-1 coefficients.
        if (n > 1) {
            first_row <- c(-rev(pars[-n]) / pars[n], 1 / pars[n])
        } else {
            first_row <- 1 / pars[n]
        }
        # Build the companion matrix
        C <- matrix(0, n, n)
        C[1, ] <- first_row
        if (n > 1) {
            C[cbind(2:n, 1:(n - 1))] <- 1
        }
        # Compute the eigenvalues of the companion matrix.
        # RTMB overloads eigen() so that it is AD‐compatible.
        proots <- eigen(C, symmetric = FALSE)$values
        cons <- 1.0 - Mod(proots)
        ADREPORT(cons)
        return(1)
    }
    
    parameters <- issmenv$parameters
    obj <- makeADFun2(func = arconstraint, parameters = parameters, map = issmenv$map, ADreport = TRUE, silent = TRUE)
    return(obj)
}

ma_jconstraint <- function(pars, fun, issmenv)
{
    maconstraint <- function(parms) {
        # Overload assignment for AD-aware indexing
        "[<-" <- ADoverload("[<-")
        pars <- -parms$issmpars[issmenv$maindex]
        n <- length(pars)
        # Construct the companion matrix.
        # For a polynomial 1 - a1*x - ... - an*x^n = 0, I want the monic form:
        #   x^n + (a_{n-1}/a_n)*x^(n-1) + ... + (a_1/a_n)*x - 1/a_n = 0.
        # The standard companion matrix for
        #   x^n + c_{n-1}*x^(n-1) + ... + c_0 = 0
        # is:
        #   [ -c_{n-1}   -c_{n-2}  ...  -c_0 ]
        #   [     1          0     ...   0  ]
        #   [     0          1     ...   0  ]
        #   [    ...       ...    ...  ... ]
        #
        # Here, c_{n-1} = a_{n-1}/a_n, c_{n-2} = a_{n-2}/a_n, ..., c_0 = -1/a_n.
        # Because the input vector is ordered as (a1, a2, ..., an), need to reverse
        # the first n-1 coefficients.
        if (n > 1) {
            first_row <- c(-rev(pars[-n]) / pars[n], 1 / pars[n])
        } else {
            first_row <- 1 / pars[n]
        }
        # Build the companion matrix
        C <- matrix(0, n, n)
        C[1, ] <- first_row
        if (n > 1) {
            C[cbind(2:n, 1:(n - 1))] <- 1
        }
        # Compute the eigenvalues of the companion matrix.
        # RTMB overloads eigen() so that it is AD‐compatible.
        proots <- eigen(C, symmetric = FALSE)$values
        cons <- 1.0 - Mod(proots)
        ADREPORT(cons)
        return(1)
    }
    parameters <- issmenv$parameters
    obj <- makeADFun2(func = maconstraint, parameters = parameters, map = issmenv$map, ADreport = TRUE, silent = TRUE)
    return(obj)
}

slope_constraint <- function(pars, fun, issmenv) {
    pars[2] - pars[1] - 0.01
}

slope_jacobian <- function(pars, fun, issmenv) {
    J <- matrix(0, ncol = length(pars), nrow = 1)
    J[1,1:2] <- c(-1, 1)
    return(J)
}


garch_constraint <- function(pars, fun, issmenv) {
    estimate <- group <- NULL
    tmp <- copy(issmenv$parmatrix)
    tmp[estimate == 1, value := pars]
    sum(tmp[group == "eta"]$value) + sum(tmp[group == "delta"]$value) - 0.99
}

garch_jacobian <- function(pars, fun, issmenv) {
    J <- matrix(0, ncol = length(pars), nrow = 1)
    J[1,issmenv$garchindex] <- 1
    return(J)
}


issm_constraint <- function(pars, fun, issmenv)
{
    allpars <- fun$env$data$allpars
    ppindex <- fun$env$data$ppindex + 1
    findex <- fun$env$data$findex + 1
    fpindex <- fun$env$data$fpindex + 1
    V <- fun$env$data$V
    fshape <- fun$env$data$fshape
    fshape[1] <- fshape[1] + 1
    fshape[3] <- fshape[3] + 1
    fshape[5] <- fshape[5] + 1
    fshape[7] <- fshape[7] + 1
    fshape[9] <- fshape[9] + 1
    issmpars <- issmenv$issmparsinx
    # add armaorder and eliminate from D
    modeli <- fun$env$data$modeli
    f <- function(parms) {
        "[<-" <- ADoverload("[<-")
        n <- length(ppindex)
        pars <- parms$issmpars
        for (i in 1:n) allpars[ppindex[i]] <- pars[i]
        m <- length(findex)
        for (i in 1:m) V[findex[i]] <- allpars[fpindex[i]]
        F0tmp <- V[fshape[1]:(fshape[1] + fshape[2] - 1)]
        F0 <- matrix(F0tmp, modeli[1], modeli[1])
        F1tmp <- V[fshape[3]:(fshape[3] + fshape[4] - 1)]
        F1 <- matrix(F1tmp, modeli[1], modeli[1])
        F2tmp <- V[fshape[5]:(fshape[5] + fshape[6] - 1)]
        F2 <- matrix(F2tmp, modeli[1], modeli[1])
        F <- matrix(as.vector(F0) * as.vector(F1) * as.vector(F2), modeli[1], modeli[1])
        W <- V[fshape[7]:(fshape[7] + fshape[8] - 1)]
        G <- V[fshape[9]:(fshape[9] + fshape[10] - 1)]
        M <- t(matrix(W, modeli[1], 1))
        B <- matrix(G, modeli[1], 1) %*% M
        D <- matrix(as.vector(F) - as.vector(B), modeli[1], modeli[1])
        D <- as.matrix(D[1:modeli[4],1:modeli[4]])
        d <- eigen(D, symmetric = FALSE, only.values = TRUE)
        tmp <- Mod(d$values)
        e <- tmp - 1.01
        REPORT(tmp)
        ADREPORT(e)
        return(sum(e))
    }
    parameters <- issmenv$parameters
    obj <- makeADFun2(func = f, parameters = parameters, map = issmenv$map, ADreport = TRUE, silent = TRUE)
    return(obj)
}

model_case <- function(spec) {
    # id [level slope|seasonal ARMA GARCH]
    id <- c(1, 0, 0, 0)
    if (spec$slope$include_slope | spec$seasonal$include_seasonal) {
        id[2] <- 1        
    }
    if (sum(spec$arma$order) > 1) {
        id[3] <- 1
    }
    if (sum(spec$variance$order) > 0) {
        id[4] <- 1
    }
    return(id)
}

make_constraint <- function(spec, fun, issmenv) {
    # case 1 only level
    # case 2 only level and GARCH
    # case 3 only level and arma
    # case 4 only level and arma and GARCH
    # case 5 level +
    # case 6 level + and GARCH
    # case 7 level + and arma
    # case 8 level + and arma and GARCH
    case_id <- model_case(spec)
    use_case <- c(0, 0, 0)
    if (case_id[2] == 0) {
        issm_c <- function(pars, fun, issmenv) {
            NULL
        }
        issm_j <- function(pars, fun, issmenv) {
            NULL
        }
    } else {
        issm_c <- function(pars, fun, issmenv) {
            issmenv$constraint$fn(pars)
        }
        issm_j <- function(pars, fun, issmenv) {
            if (all(abs(issmenv$constraint$report(pars)$tmp - 1) <= 1e-12)) {
                J <- matrix(1, ncol = length(pars), nrow = 1)
            } else {
                J <- issmenv$constraint$gr(pars)
            }
            return(J)
        }
        use_case[1] <- 1
    }
    if (case_id[3] > 0) {
        arma_order <- spec$arma$order
        if (arma_order[1] > 1 & arma_order[2] > 1) {
            arma_c <- function(pars, fun, issmenv) {
                c(issmenv$arcons$fn(pars),issmenv$macons$fn(pars))
            }
            arma_j <- function(pars, fun, issmenv) {
                rbind(issmenv$arcons$gr(pars), issmenv$macons$gr(pars))
            }
        } else if (arma_order[1] > 1 & arma_order[2] <= 1) {
            arma_c <- function(pars, fun, issmenv) {
                issmenv$arcons$fn(pars)
            }
            arma_j <- function(pars, fun, issmenv) {
                issmenv$arcons$gr(pars)
            }
        } else if (arma_order[1] <= 1 & arma_order[2] > 1) {
            arma_c <- function(pars, fun, issmenv) {
                issmenv$macons$fn(pars)
            }
            arma_j <- function(pars, fun, issmenv) {
                issmenv$macons$gr(pars)
            }
        } else {
            arma_c <- function(pars, fun, issmenv) {
                NULL
            }
            arma_j <- function(pars, fun, issmenv) {
                NULL
            }
        }
        use_case[2] <- 1
    } else {
        arma_c <- function(pars, fun, issmenv) {
            NULL
        }
        arma_j <- function(pars, fun, issmenv) {
            NULL
        }
    }
    if (case_id[4] > 0) {
        garch_c <- function(pars, fun, issmenv) {
            garch_constraint(pars, fun, issmenv)
        }
        garch_j <- function(pars, fun, issmenv){
            garch_jacobian(pars, fun, issmenv)
        }
        use_case[3] <- 1
    } else {
        garch_c <- function(pars, fun, issmenv) {
            NULL
        }
        garch_j <- function(pars, fun, issmenv){
            NULL
        }
    }
    
    C <- list()
    J <- list()
    k <- 0
    if (use_case[1] == 1) {
        k <- k + 1
        C[[k]] <- issm_c
        J[[k]] <- issm_j
    }
    if (use_case[2] == 1) {
        k <- k + 1
        C[[k]] <- arma_c
        J[[k]] <- arma_j
    }
    if (use_case[3] == 1) {
        k <- k + 1
        C[[k]] <- garch_c
        J[[k]] <- garch_j
    }
    
    if (sum(use_case) == 0) {
        constraint_fun <- NULL
    } else {
        constraint_fun <- function(pars, fun, issmenv) {
            con <- unlist(sapply(C, do.call, list(pars, fun, issmenv)))
            jac <- do.call(rbind, lapply(J, function(x) x(pars, fun, issmenv)))
            return(list("constraints" = con, "jacobian" = jac))
        }
    }
    return(constraint_fun)
}

companion_matrix <- function(params) {
    n <- length(params)
    if (params[n] == 0) stop("The last parameter must not be zero to form a monic polynomial.")
    # Reverse the order of params[-n] so that the coefficient for x^(n-1) comes first.
    if(n > 1){
        first_row <- c(-rev(params[-n]) / params[n], 1 / params[n])
    } else {
        first_row <- 1 / params[n]
    }
    # Build an n x n matrix and fill the first row
    C <- matrix(0, n, n)
    C[1, ] <- first_row
    # Fill in the sub-diagonal with ones
    if(n > 1){
        C[cbind(2:n, 1:(n-1))] <- 1
    }
    return(C)
}

