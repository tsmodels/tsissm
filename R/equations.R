.distribution_abb <- function(distribution, garch = FALSE, standardized = FALSE)
{
    out <- switch(distribution,
                  "norm" = "N\\left(0,1\\right)",
                  "snorm" = "\\text{SN}\\left(0,1,\\zeta\\right)",
                  "std" = "\\text{T}\\left(0,1,\\nu\\right)",
                  "sstd" = "\\text{ST}\\left(0,1,\\zeta, \\nu\\right)",
                  "ged" = "\\text{GED}\\left(0,1,\\nu\\right)",
                  "sged" = "\\text{SGED}\\left(0,1,\\zeta, \\nu\\right)",
                  "jsu" = "\\text{JSU}\\left(0,1,\\zeta, \\nu\\right)",
                  "nig" = "\\text{NIG}\\left(0,1,\\zeta, \\nu\\right)",
                  "ghst" = "\\text{GHST}\\left(0,1,\\zeta, \\nu\\right)",
                  "gh" = "\\text{GH}\\left(0,1,\\zeta, \\nu,\\lambda\\right)")
    if (!standardized) {
        if (!garch) {
            out <- gsub("0,1","0,\\\\sigma",out)
        } else {
            out <- gsub("0,1","0,\\\\sigma_t",out)
            
        }
    }
    return(out)
}

.equation_regressors <- function(xreg, include_xreg)
{
    eq1 <- NULL
    if (include_xreg) {
        n <- NCOL(xreg)
        if (n == 1) {
            eq1 <- "\\kappa_1 \\xi_{1,t}"
        } else {
            eq1 <- paste0("\\sum_{j=1}^{",n,"} \\kappa_j \\xi_{j,t}")
        }
    }
    return(list(regressor = eq1))
}

.equation_garch <- function(garch_order) {
    eq_constant <- "\\hat\\sigma^2 (1 - P)"
    v_equation <- paste0("\\sigma^2_{t} = ", eq_constant)
        
    if (garch_order[1] > 0) {
        if (garch_order[1] > 1) {
            eq_alpha <- paste0("\\sum_{j=1}^",garch_order[1],"\\eta_j \\varepsilon^2_{t-j}")
        } else {
            eq_alpha <- paste0("\\eta_1\\varepsilon^2_{t-1}")
        }
    } else {
        eq_alpha <- NULL
    }
    if (garch_order[2] > 0) {
        if (garch_order[2] > 1) {
            eq_beta <- paste0("\\sum_{j=1}^",garch_order[2],"\\deltq_j \\sigma^2_{t-j}")
        } else {
            eq_beta <- paste0("\\delta_1\\sigma^2_{t-1}")
        }
    } else {
        eq_beta <- NULL
    }
    v_equation <- paste(v_equation, eq_alpha, eq_beta, sep = "+")
    return(v_equation)
}

.equation_issm <- function(spec)
{
    if (spec$variance$type == "constant") {
        dist <- .distribution_abb(spec$distribution$type, garch = FALSE, standardized = FALSE)
        eq_variance <- NULL
    } else {
        dist <- .distribution_abb(spec$distribution$type, garch = TRUE, standardized = FALSE)
        eq_variance <- .equation_garch(spec$variance$order)
    }
    if (spec$xreg$include_xreg) {
        eq_reg <- .equation_regressors(spec$xreg$xreg, spec$xreg$include_xreg)
        eq_obs <- paste0("y^{\\omega}_t=w'x_{t-1} + ",eq_reg," + \\varepsilon, \\varepsilon \\sim ",dist)
    } else {
        eq_obs <- paste0("y^{\\omega}_t=w'x_{t-1} + \\varepsilon, \\varepsilon \\sim ",dist)
    }
    eq_state <- paste("x_t = Fx_{t-1}+g\\varepsilon_{t}")
    return(list(observation = eq_obs, state = eq_state, variance = eq_variance))
}
