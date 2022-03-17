\name{predict}
\alias{predict}
\alias{predict.tsissm.estimate}
\title{
Model Prediction
}
\description{
Prediction function for class \dQuote{tsissm.estimate}.
}
\usage{
\method{predict}{tsissm.estimate}(object, h = 12, newxreg = NULL,
nsim = 1000, forc_dates = NULL, innov = NULL,
init_states = NULL, exact_moments = TRUE, innov_type = "q", 
sigma_scale = NULL, ...)
}
\arguments{
\item{object}{
An object of class tsissm.estimate
}
\item{h}{
The forecast horizon.
}
\item{newxreg}{
The matrix of external regressors in the forecast horizon.
}
\item{nsim}{
The number of simulations to use for generating the simulated
predictive distribution.
}
\item{forc_dates}{
An optional vector of forecast dates equal to h. If NULL will
use the implied periodicity of the data to generate a regular
sequence of dates after the last available date in the data.
}
\item{init_states}{
An optional vector of states to initialize the forecast. If NULL,
will use the last available state from the estimated model.
}
\item{innov}{
This is an optional vector of uniform innovations which will be
translated to regular innovations using the appropriate
distribution quantile function and model standard deviation.
The length of this vector should be equal to nsim x horizon.
}
\item{exact_moments}{
Whether to rescale the mean and variance of the simulated distribution by their 
exact (analytic) moments. This is performed on the transformed data.
}
\item{innov_type}{
If \sQuote{innov} is not NULL, then this denotes the type of values passed,
with \dQuote{q} denoting quantile probabilities (default and backwards
compatible) and \dQuote{z} for standardized errors.
}
\item{sigma_scale}{
A vector of length h denoting a scaling factor which is applied to rescale the
standard deviation of each simulated horizon's distribution.
}
\item{\dots}{
Not currently used
}
}
\value{
An object of class \dQuote{tsissm.predict} which also inherits \dQuote{tsmodel.predict}, with slots
for the simulated prediction distribution, the original series (as a zoo object), the original
specification object and the mean forecast. The predictive distribution is back transformed if lambda was
not set to NULL in the specification.
}
\author{
Alexios Galanos
}
