\name{issm_modelspec}
\alias{issm_modelspec}
\title{
Model Specification
}
\description{
Specifies an issm model prior to estimation.
}
\usage{
issm_modelspec(y, slope = TRUE, slope_damped = FALSE, seasonal = FALSE,
seasonal_frequency = 1, seasonal_type = c("trigonometric", "regular"),
seasonal_harmonics = NULL, ar = 0, ma = 0, xreg = NULL, lambda = 1,
sampling = NULL, ...)
}
\arguments{
\item{y}{
An xts vector.
}
\item{slope}{
(Logical) Slope component.
}
\item{slope_damped}{
(Logical) Slope Dampening component.
}
\item{seasonal}{
(Logical) Seasonal component(s).
}
\item{seasonal_frequency}{
Vector of numeric seasonal frequency. For trigonometric this can be fractional, but
must be integer for regular seasonality.
}
\item{seasonal_type}{
Either trigonometric or regular. The latter currently does not allow multiple
seasonality.
}
\item{seasonal_harmonics}{
The number of harmonics per seasonal frequency for the trigonometric type.
}
\item{ar}{
The AR order.
}
\item{ma}{
The MA order.
}
\item{xreg}{
An xts matrix of external regressors.
}
\item{lambda}{
The Box Cox lambda. If not NULL, then either a numeric value or NA denoting automatic calculation.
}
\item{sampling}{
(Optional) sampling frequency of the dataset. If NULL, will try to identify from the timestamps of
y. This is useful for plotting and extending the timestamps in the prediction horizon.
}
\item{\dots}{
No other arguments at present.
}
}
\details{
The specification performs some sanity checks on the arguments provided and sets up the
required state space matrices and parameters which are used in the estimation stage.
}
\value{
An object of class \dQuote{tsissm.spec} with the following slots:\cr
\item{target}{A list with original data series, the data series index and the sampling
frequency}
\item{slope}{A list with details about the slope state}
\item{seasonal}{A list with details about the seasonal state}
\item{xreg}{A list with details on the external regressors}
\item{transform}{A list with details on the Box Cox transformation}
\item{arma}{A list with details on the ARMA state}
\item{S}{A data.table with the vectorized state matrices}
\item{dims}{A vector with dimensions and flags used in the estimation code}
\item{parmatrix}{A data.table of the model parameters}
\item{idmatrix}{A matrix with index information on the parameters}
}
\references{
De Livera, Alysha M and Hyndman, Rob J and Snyder, Ralph D, 2011,
Forecasting time series with complex seasonal patterns using exponential smoothing,
\emph{Journal of the American Statistical Association}, \bold{106(496)}, 1513--1527.
}
\author{
Alexios Galanos
}
\keyword{"ts"}
