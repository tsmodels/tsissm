issm_control <- function(solver = "nloptr", algorithm = c("SLSQP","AUGLAG/MMA","AUGLAG/CCSAQ"), trace = 0, xtol_rel = 1e-14, maxeval = 1000, xtol_abs = 1e-12)
{
    algorithm <- match.arg(algorithm[1], c("SLSQP","AUGLAG/MMA","AUGLAG/CCSAQ"))
    switch(algorithm,
           "SLSQP" = list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_SLSQP", xtol_rel = xtol_rel, maxeval = maxeval, xtol_abs = xtol_abs, check_derivatives = FALSE),
           "AUGLAG/MMA" = list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_AUGLAG", xtol_rel = xtol_rel, xtol_abs = xtol_abs, maxeval = maxeval, check_derivatives = FALSE,
                              local_opts = list(algorithm = "NLOPT_LD_MMA", maxeval = 500, xtol_rel = 1e-12)),
           "AUGLAG/CCSAQ" = list(print_level = ifelse(trace, 1, 0), algorithm = "NLOPT_LD_AUGLAG", xtol_rel = xtol_rel, xtol_abs = xtol_abs, maxeval = maxeval, check_derivatives = FALSE,
                                 local_opts = list(algorithm = "NLOPT_LD_CCSAQ", maxeval = 500, xtol_rel = 1e-12))
    )
}