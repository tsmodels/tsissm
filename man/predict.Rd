% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.R
\name{predict.tsissm.estimate}
\alias{predict.tsissm.estimate}
\alias{predict}
\title{Model Prediction}
\usage{
\method{predict}{tsissm.estimate}(
  object,
  h = 12,
  newxreg = NULL,
  nsim = 1000,
  forc_dates = NULL,
  innov = NULL,
  innov_type = "q",
  init_states = NULL,
  exact_moments = TRUE,
  sigma_scale = NULL,
  ...
)
}
\arguments{
\item{object}{an object of class \dQuote{tsissm.estimate}.}

\item{h}{the forecast horizon.}

\item{newxreg}{a matrix of external regressors in the forecast horizon.}

\item{nsim}{the number of simulations to use for generating the simulated
predictive distribution.}

\item{forc_dates}{an optional vector of forecast dates equal to h. If NULL will
use the implied periodicity of the data to generate a regular sequence of
dates after the last available date in the data.}

\item{innov}{an optional vector of uniform innovations which will be translated
to regular innovations using the appropriate distribution quantile function
and model standard deviation. The length of this vector should be equal to
nsim x horizon.}

\item{innov_type}{if \sQuote{innov} is not NULL, then this denotes the type of values
passed, with \dQuote{q} denoting quantile probabilities (default and
backwards compatible) and \dQuote{z} for standardized errors.}

\item{init_states}{an optional vector of states to initialize the forecast.
If NULL, will use the last available state from the estimated model.}

\item{exact_moments}{whether to rescale the mean and variance of the simulated 
distribution by their exact (analytic) moments. This is performed on the 
transformed data.}

\item{sigma_scale}{a vector of length h denoting a scaling factor which is
applied to rescale the standard deviation of each simulated horizon's
distribution.}

\item{...}{not currently used.}
}
\value{
An object of class \dQuote{tsissm.predict} which also inherits
\dQuote{tsmodel.predict}, with slots for the simulated prediction distribution,
the original series (as a zoo object), the original specification object and
the mean forecast. The predictive distribution is back transformed if lambda was
not set to NULL in the specification.
}
\description{
Prediction function for class \dQuote{tsissm.estimate}.
}
\details{
Like all models in the ts framework, prediction is done by
simulating h-steps ahead in order to build a predictive distribution.
}
