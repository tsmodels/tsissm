#' @rawNamespace useDynLib(tsissm,.registration = TRUE); useDynLib(tsissm_TMBExports)
#' @keywords internal
#' @import methods
#' @import tsmethods
#' @import data.table
#' @importFrom RTMB AD ADoverload ADREPORT REPORT matrix eigen as.vector.advector
#' @importFrom TMB MakeADFun
#' @importFrom utils head tail data txtProgressBar setTxtProgressBar
#' @importFrom stats sd var acf na.pass pchisq pgamma qqline qqnorm simulate na.omit median fitted coef quantile residuals predict logLik cov cor qt pnorm AIC BIC nobs sigma ar arima dnorm printCoefmat vcov density
#' @importFrom graphics grid legend lines par plot points abline axis axis.Date axis.POSIXct box polygon layout mtext title hist boxplot
#' @importFrom grDevices gray n2mfrow
#' @importFrom zoo index as.zoo zoo coredata na.locf na.fill is.zoo
#' @importFrom sandwich estfun bwNeweyWest vcovHAC vcovOPG bread
#' @importFrom nloptr nloptr
#' @importFrom xts xts as.xts is.xts
#' @importFrom flextable flextable as_flextable set_caption italic fontsize separate_header add_footer_row add_footer_lines append_chunks as_chunk as_equation as_paragraph compose colformat_double set_header_labels padding bold align autofit hline width
#' @importFrom future.apply future_lapply
#' @importFrom future %<-% future value nbrOfWorkers
#' @importFrom progressr handlers progressor
#' @importFrom viridisLite viridis
#' @importFrom Rcpp evalCpp loadModule
#' @importFrom copula rCopula normalCopula P2p
#' @importFrom tsaux future_dates sampling_frequency box_cox crps auto_clean tstransform tslinear mape mase crps bias mslre mis
#' @importFrom tsdistributions distribution_modelspec ddist pdist rdist qdist
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
