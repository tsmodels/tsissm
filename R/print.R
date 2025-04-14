pvalue_format <- function(x) {
    z <- cut(x, breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf), labels = c("***", "**", "*", ".", ""))
    as.character(z)
}

.print_screen <- function(x, digits = max(3L, getOption("digits") - 3L), signif.stars = getOption("show.signif.stars"), ...)
{
    term <- NULL
    df <- x$df
    rdf <- x$n_parameters + 1
    coefs <- copy(x$coefficients)
    coef_names <- coefs$term
    coefs <- coefs[,term := NULL]
    coefs <- as.data.frame(coefs)
    rownames(coefs) <- coef_names
    cat("ISSM Model:", x$model)
    cat("\nCoefficients:\n")
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    if (x$n_missing > 0) {
        cat("\nN:", as.integer(x$n_obs), " | Missing :", x$n_missing)
    } else {
        cat("\nN:", as.integer(x$n_obs))
    }
    cat("\nsigma^2:", format(signif(x$init_variance, digits = digits)))
    #cat("\nStates [0]:", format(signif(x$init_states, digits = digits)))
    cat("\nLogLik:", format(signif(x$loglikelihood, digits = 2 + digits)))
    cat(",  ")
    cat("AIC: ", format(signif(x$AIC, digits = 2 + digits)))
    cat(",  ")
    cat("BIC:", format(signif(x$BIC, digits = 2 + digits)))
    cat("\n")
    cat("DAC (%):", format(signif(x$DAC, digits = 2 + digits)))
    cat("\n")
    cat("MAPE (%):", format(signif(x$MAPE, digits = 2 + digits)))
    cat("\n")
    invisible(x)
}


.print_flextable <- function(x, digits = max(3L, getOption("digits") - 3L),
                             signif.stars = getOption("show.signif.stars"),
                             include.symbols = TRUE, include.equation = TRUE,
                             include.statistics = TRUE, table.caption = paste0("ISSM Model: ", x$model), ...)
{
    signif <- `Pr(>|t|)` <- NULL
    n <- nrow(x$coefficients)
    cf <- copy(x$coefficients)
    k <- 0
    if (signif.stars) {
        cf[,signif := pvalue_format(`Pr(>|t|)`)]
        out <- flextable(cf) |> set_caption(caption = table.caption) |>
            align(j = "term", align = "left") |>
            align(j = "signif", align = "left") |>
            padding(padding.right = 0, j = "`Pr(>|t|)`", part  = "all") |>
            bold(j = "signif", bold = TRUE) |>
            padding(padding.left = 0, j = "signif", part  = "all") |>
            set_header_labels(term = "", Estimate = "Estimate",
                              `Std. Error` = "Std. Error", `t value` = "t value",
                              `Pr(>|t|)` = "Pr(>|t|)", signif = "" )
        out <- out |> add_footer_lines(values = c("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1"))
        k <- 1
    } else {
        out <- flextable(cf) |> set_caption(caption = table.caption) |>
            align(j = "term", align = "left") |>
            set_header_labels(term = "", Estimate = "Estimate",
                              `Std. Error` = "Std. Error", `t value` = "t value",
                              `Pr(>|t|)` = "Pr(>|t|)")
    }
    if (include.symbols) {
        for (i in 1:n) {
            out <- compose(out, i = i, j = 1, as_paragraph(as_chunk(' '))) |> append_chunks(i = i,j = 1, as_equation(x$symbol[i]))
        }
    }
    if (include.statistics) {
        out <- out |> add_footer_lines(values = c(paste0(sprintf("sigma^2: %s",formatC(x$init_variance))),
                                                  paste0("LogLik: ", format(x$loglikelihood)),
                                                  sprintf("AIC: %s | BIC: %s", formatC(x$AIC), formatC(x$BIC)),
                                                  sprintf("DAC : %s | MAPE : %s", formatC(x$DAC, digits = 2), formatC(x$MAPE, digits = 2))))
        k <- k + 5
    }
    if (include.equation) {
        if (k == 0) flag <- TRUE else flag <- FALSE
        out <- out |> add_footer_lines(top = FALSE, values = "Model Equation")
        out <- out |> add_footer_row(top = FALSE, values = " ", colwidths = length(out$col_keys))
        k <- k + 1
        out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k, j = "term", part = "footer", as_equation(paste0(x$equation$observation)))
        out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k + 1, j = "term", part = "footer", as_equation(paste0(x$equation$state)))
        if (x$variance_type == "dynamic") {
            out <- out |> add_footer_lines(values = " ", top = FALSE) |> append_chunks(i = k + 2, j = "term", part = "footer", as_equation(paste0(x$equation$variance)))
            k <- k + 1
        }
        if (!flag) out <- out |> hline(part = "footer",i = k - 3, j = 1)
        out <- out |> hline(part = "footer",i = k - 3, j = 1)
    }
    out <- colformat_double(out, digits = digits) |> autofit()
    return(out)
}