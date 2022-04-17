\name{auto_issm}
\alias{auto_issm}
\title{
Automatic Model Selection
}
\description{
Automatic model selection based on enumeration of options space
}
\usage{
auto_issm(y, slope = c(TRUE, FALSE), slope_damped = c(TRUE, FALSE), 
seasonal = c(TRUE, FALSE), seasonal_frequency = 1, 
seasonal_type = "trigonometric", seasonal_harmonics = list(), 
ar = 0:2, ma = 0:2, xreg = NULL, transformation = "box-cox", lambda = 1, 
lower = 0, upper = 1, sampling = NULL, trace = FALSE, 
return_table = FALSE, solver = "nlminb", autodiff = TRUE, ...)
}
\arguments{
\item{y}{
An xts vector
}
\item{slope}{
Whether to include or not the slope. Having both TRUE and FALSE will estimate
all possible combination with and without a slope.
}
\item{slope_damped}{
Whether to include or not the slope dampening. Having both TRUE and FALSE 
will estimate all possible combination with and without a slope dampening.
}
\item{seasonal}{
Whether to include or not a seasonal component. Having both TRUE and FALSE 
will estimate all possible combination with and without a seasonal component.
}
\item{seasonal_frequency}{
A vector of seasaonal frequencies.
}
\item{seasonal_type}{
The type of seasonality to include. Trigonometric is the preferred approach
and the only one accepting multiple seasonal frequencies.
}
\item{seasonal_harmonics}{
A list with slots for each seasonal frequency listing the sequence of harmonics
to test for each (see details).
}
\item{ar}{
A vector of the ar terms to test for (see details).
}
\item{ma}{
A vector of the ma terms to test for (see details).
}
\item{xreg}{
An xts matrix of exogenous regressors to include in the model.
}
\item{transformation}{
The transformation to use (defaults to box-cox and is effectively NULL if lambda
is NULL or 1).
}
\item{lambda}{
The box-cox lambda parameter.
}
\item{lower}{
The lower bound for the transformation.
}
\item{upper}{
The upper bound for the transformation.
}
\item{sampling}{
The sampling frequency of the data.
}
\item{trace}{
Whether to show the progress bar. The user is expected to have set up 
appropriate handlers for this using the \dQuote{progressr} package.
}
\item{return_table}{
Whether to return the table with the enumerated options and the AIC and MAPE
for each combination of those options used.
}
\item{solver}{
The solver to use.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. 
This makes use of the tsissmad package.
}
\item{\dots}{
Not currently used.
}
}
\details{
The user is responsible for passing reasonable options for the options. For 
instance, the harmonics must be strictly positive and less than one half the
seasonal frequency. The \code{\link{expand.grid}} function is used to enumerate
all possible combinations of the options with some sanity checks and 
eliminations if testing for both seasonal and non-seasonal models, or for the
case of slope and dampening.
If the user prefers to ensemble models or use instead the MAPE criterion, the
\sQuote{return_table} should be set to TRUE and the returned table used to
estimate an alternative to the top AIC model.
}
\value{
An object of class \dQuote{tsissm.estimate} returning the best model based on
AIC (minimum).
}
\note{
The function can use parallel functionality as long as the user has set up a
\code{\link[future]{plan}} using the future package.
}
\author{
Alexios Galanos
}
