#' ISSM solver control parameters
#'
#' @param algorithm (nloptr) the algorithm to use
#' @param trace (integer) controls print level information
#' @param xtol_rel (nloptr) relative tolerance on optimization parameters
#' @param maxeval (nloptr) number of function evaluations to stop on
#' @param xtol_abs (nloptr) absolute tolerances on optimization parameters
#' @param max_iter (solnp) maximum number of major (outer) iterations
#' @param min_iter (solnp) maximum number of minor (inner) iterations per outer iteration
#' @param rho (solnp) onitial penalty parameter for the augmented Lagrangian
#' @param tol (solnp) convergence tolerance
#' @details
#' The two function provide defaults for use when using either the \code{\link{nloptr}[nloptr]} or
#' \code{\link{csolnp}[Rsolnp]} functions. For the former, additional control parameters may be
#' appended to the list if the user so wishes (\sQuote{nloptr} has many more options).
#' @returns a list with the options which is then passed to the nloptr or solnp solver.
#' @rdname solver_control
#' @export
issm_control <- function(solver = "nloptr", algorithm = c("SLSQP","AUGLAG/MMA","AUGLAG/CCSAQ"), 
                         trace = 0, xtol_rel = 1e-14, maxeval = 1000, xtol_abs = 1e-12,
                         rho = 1, max_iter = 400, min_iter = 800, tol = 1e-8)
{
    solver <- match.arg(solver[1], c("nloptr","solnp"))
    if (solver == "nloptr") {
        algorithm <- match.arg(algorithm[1], c("SLSQP","AUGLAG/MMA","AUGLAG/CCSAQ"))
        maxeval <- max(1, maxeval)
        xtol_abs <- abs(xtol_abs)
        trace <- max(0, trace)
        xtol_rel <- abs(xtol_rel)
        ctrl <- switch(algorithm,
               "SLSQP" = list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_SLSQP", xtol_rel = xtol_rel, maxeval = maxeval, xtol_abs = xtol_abs, check_derivatives = FALSE),
               "AUGLAG/MMA" = list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_AUGLAG", xtol_rel = xtol_rel, xtol_abs = xtol_abs, maxeval = maxeval, check_derivatives = FALSE,
                                  local_opts = list(algorithm = "NLOPT_LD_MMA", maxeval = maxeval, xtol_rel = xtol_rel)),
               "AUGLAG/CCSAQ" = list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_AUGLAG", xtol_rel = xtol_rel, xtol_abs = xtol_abs, maxeval = maxeval, check_derivatives = FALSE,
                                     local_opts = list(algorithm = "NLOPT_LD_CCSAQ", maxeval = maxeval, xtol_rel = xtol_rel))
        )
        return(ctrl)
    } else if (solver == "solnp") {
        max_iter <- max(1, max_iter)
        min_iter <- min(1, min_iter)
        trace <- max(0, trace)
        tol <- abs(tol)
        ctrl <- list(trace = trace, max_iter = max_iter, min_iter = min_iter, tol = tol)
        return(ctrl)
    }
}
