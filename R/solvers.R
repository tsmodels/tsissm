#' ISSM solver control parameters
#'
#' @param solver only \sQuote{nloptr} at present
#' @param algorithm (nloptr) the algorithm to use
#' @param trace (integer) controls print level information
#' @param xtol_rel (nloptr) relative tolerance on optimization parameters
#' @param maxeval (nloptr) number of function evaluations to stop on
#' @param xtol_abs (nloptr) absolute tolerances on optimization parameters
#' @details
#' The function provides defaults for use, but additional control parameters may 
#' be appended to the list if the user so wishes (\sQuote{nloptr} has many more options).
#' @returns a list with the options which is then passed to the solver.
#' @rdname solver_control
#' @export
issm_control <- function(solver = "nloptr", algorithm = c("SLSQP","AUGLAG/MMA","AUGLAG/CCSAQ"), 
                         trace = 0, xtol_rel = 1e-14, maxeval = 1000, xtol_abs = 1e-12)
{
    solver <- match.arg(solver[1], c("nloptr"))
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
    }
}
