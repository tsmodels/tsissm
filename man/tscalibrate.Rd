% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calibration.R
\name{tscalibrate.tsissm.spec}
\alias{tscalibrate.tsissm.spec}
\alias{tscalibrate}
\title{Walk Forward Calibration of Distribution Variance (Experimental)}
\usage{
\method{tscalibrate}{tsissm.spec}(
  object,
  start = floor(length(object$target$y_orig))/2,
  end = length(object$target$y_orig),
  h = 1,
  nsim = 5000,
  solver = "nlminb",
  autodiff = TRUE,
  autoclean = FALSE,
  trace = FALSE,
  ...
)
}
\arguments{
\item{object}{an object of class \dQuote{tsissm.spec}.}

\item{start}{numeric data index from which to start the backtest calibration.}

\item{end}{numeric data index on which to end the backtest. The backtest will
end 1 period before that date in order to have at least 1 out of sample value
to compare against.}

\item{h}{forecast horizon. As the expanding window approaches the \dQuote{end},
the horizon will automatically shrink to the number of available out of sample
periods.}

\item{nsim}{number of samples to draw for the simulated prediction distribution.}

\item{solver}{solver to use.}

\item{autodiff}{whether to use automatic differentiation for estimation. This
makes use of the tsissmad package.}

\item{autoclean}{whether to perform automatic cleaning on the training data
prior to prediction as per the \sQuote{auto_clean} function in the tsaux package.}

\item{trace}{whether to show the progress bar. The user is expected to have set
up appropriate handlers for this using the \dQuote{progressr} package.}

\item{...}{Additional arguments passed to the \dQuote{auto_clean} function.}
}
\value{
A list with the following data.tables:
\itemize{
\item sample : multi-horizon white noise sample
\item sigma_scale: The scaling values per horizon which can be used to scale
the predictive distribution standard deviation.
\item sd_sigma_scale: The standard deviation of the sigma scaling factor
}
}
\description{
Using an expanding window walk forward backtest, provides a
calibrated set of adjustment values which can be used in the predict function
for scaling the multi-horizon distributional standard deviation.
}
\note{
The walk forward predictive errors per horizon are used to estimate the
realized variance against the expected variance (sigma_scale). A jackknife
procedure is used for estimating the mean of this value. The output also
includes samples from a kernel fitted on the whitened errors per horizon
which can then be used as inputs to the prediction function. The function can
use parallel functionality as long as the user has set up a
\code{\link[future]{plan}} using the future package.
}
