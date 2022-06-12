#' @rawNamespace useDynLib(tsissm,.registration = TRUE)
#' @keywords internal
#' @import methods
#' @import tsmethods
#' @import data.table
#' @importFrom utils head tail data txtProgressBar setTxtProgressBar
#' @importFrom stats lm optimize qnorm rnorm runif sd var Box.test acf na.pass pchisq pgamma qqline qqnorm simulate
#' @importFrom tsaux future_dates sampling_frequency box_cox mape mase mslre mis bias crps auto_clean tstransform
#' @importFrom stats median na.omit fitted coef quantile residuals predict optim nlminb logLik cov cor cov2cor pnorm AIC
#' @importFrom graphics grid legend lines par plot points abline axis axis.Date axis.POSIXct box polygon layout mtext title hist
#' @importFrom grDevices gray colorRampPalette n2mfrow
#' @importFrom zoo index as.zoo zoo coredata na.locf
#' @importFrom Rsolnp solnp
#' @importFrom knitr kable
#' @importFrom tsissmad estimate_ad.tsissm.spec
#' @importFrom xts xts as.xts is.xts endpoints
#' @importFrom EnvStats rosnerTest
#' @importFrom future.apply future_lapply
#' @importFrom future %<-%
#' @importFrom progressr handlers progressor
#' @importFrom viridis viridis_pal
#' @importFrom corpcor make.positive.definite
#' @importFrom ks rkde kde
#' @importFrom bootstrap jackknife
#' @importFrom Rcpp evalCpp loadModule
"_PACKAGE"
