\name{residuals}
\alias{residuals}
\alias{residuals.tsissm.estimate}
\title{
Model Residuals
}
\description{
Extract the residual values from an estimated model.
}
\usage{
\method{residuals}{tsissm.estimate}(object, raw = FALSE,
h = 1, seed = NULL, trace = FALSE, index_start = 1, ...)
}
\arguments{
\item{object}{
An object of class tsissm.estimate
}
\item{raw}{
Whether to return residuals on the scale of the Box Cox transformation.
}
\item{h}{
The horizon (steps) ahead residuals required. The default represents the
standard residuals whilst for h>1 these are the (1:h)-step ahead in-sample
predicted residuals for each time point under fixed coefficients.
}
\item{seed}{
A seed value which initializes the simulated predictive distribution from
which the h-step ahead forecasts are made in order to calculate the residuals.
}
\item{trace}{
Whether to show the progress bar. The user is expected to have set up 
appropriate handlers for this using the \dQuote{progressr} package.
}
\item{index_start}{
The numerix index of the series from which to start the evaluation (defaults
to the first data point). For very large series, one may be interested in
discarding earlier periods.
}
\item{\dots}{
Not currently used
}
}
\details{
For h>1, this is like performing an in-sample backtest starting at
time 1 with fixed coefficients. The purpose of having the matrix of
h-step ahead residuals is in order to calculate the 1:h covariance
matrix as well as the cross 1:h covariance matrix when ensembling
series at multiple horizons.
}
\value{
An xts vector of the model residuals for h = 1, else a data.table
with rows representing the first prediction date and columns the
h-ahead forecast residuals.
}
\note{
The function can use parallel functionality (for h > 1) as long as the user has 
set up a \code{\link[future]{plan}} using the future package.
}
\author{
Alexios Galanos
}
