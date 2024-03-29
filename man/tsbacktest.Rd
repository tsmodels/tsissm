% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/backtest.R
\name{tsbacktest.tsissm.spec}
\alias{tsbacktest.tsissm.spec}
\alias{tsbacktest}
\title{Walk Forward Model Backtest}
\usage{
\method{tsbacktest}{tsissm.spec}(
  object,
  start = floor(length(object$target$y_orig)/2),
  end = length(object$target$y_orig),
  h = 1,
  estimate_every = 1,
  FUN = NULL,
  alpha = NULL,
  solver = "nlminb",
  autodiff = TRUE,
  autoclean = FALSE,
  trace = FALSE,
  ...
)
}
\arguments{
\item{object}{an object of class \dQuote{tsissm.spec}.}

\item{start}{numeric data index from which to start the backtest.}

\item{end}{numeric data index on which to end the backtest. The backtest will
end 1 period before that date in order to have at least 1 out of sample value
to compare against.}

\item{h}{forecast horizon. As the expanding window approaches the \dQuote{end},
the horizon will automatically shrink to the number of available out of sample
periods.}

\item{estimate_every}{number of periods at which the model is re-estimated
and new predictions are generated (defaults to 1).}

\item{FUN}{optional function which is applied across all horizons for each
draw (i.e. operating on each row of the distribution which represents a
single path). For example, using the max function will return the distribution
of the maximum across all horizons, whereas sum (for flow variables) would
represent the aggregate value distribution. The P50 of this distribution is
returned and aligned with the terminal horizon for each re-estimation period,
and if alpha is not NULL, then the quantiles of this distributions with
respect to the coverage (alpha) chosen.}

\item{alpha}{optional numeric vector of coverage rates for which to calculate
the quantiles.}

\item{solver}{solver to use.}

\item{autodiff}{whether to use automatic differentiation for estimation.
This makes use of the tsissmad package.}

\item{autoclean}{whether to perform automatic cleaning on the training data
prior to prediction as per the \sQuote{auto_clean} function in the tsaux package.}

\item{trace}{whether to show the progress bar. The user is expected to have
set up appropriate handlers for this using the \dQuote{progressr} package.}

\item{...}{additional arguments passed to the \dQuote{auto_clean} function.}
}
\value{
A list with the following data.tables:
\itemize{
\item prediction : the backtest table with forecasts and actuals
\item metrics: a summary performance table showing metrics by
forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
}
}
\description{
Generates an expanding window walk forward backtest.
}
\note{
The function can use parallel functionality as long as the user has
set up a \code{\link[future]{plan}} using the future package.
}
