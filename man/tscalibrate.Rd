\name{tscalibrate}
\alias{tscalibrate}
\alias{tscalibrate.tsissm.spec}
\title{
Walk Forward Calibration of Distribution Variance (Experimental)
}
\description{
Using an expanding window walk forward backtest, provides a calibrated
set of adjustment values which can be used in the predict function for
scaling the multi-horizon distributional standard deviation.
}
\usage{
\method{tscalibrate}{tsissm.spec}(object,
start = floor(length(object$target$y_orig))/2,
end = length(object$target$y_orig), h = 1, nsim = 5000, cores = 1,
solver = "nlminb", autodiff = FALSE, autoclean = FALSE, trace = FALSE,  ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsissm.spec}.
}
\item{start}{
The numeric data index from which to start the backtest calibration.
}
\item{end}{
The numeric data index on which to end the backtest. The backtest
will end 1 period before that date in order to have at least 1 out
of sample value to compare against.
}
\item{h}{
The forecast horizon. As the expanding window approaches the \dQuote{end},
the horizon will automatically shrink to the number of available out of
sample periods.
}
\item{nsim}{
Number of sample to draw for the simulated prediction distribution.
}
\item{cores}{
The number of cores to use for parallel estimation. The package
makes use of the doSNOW package.
}
\item{solver}{
The solver to use.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsetsad package.
}
\item{autoclean}{
Whether to perform automatic cleaning on the training data prior to prediction
as per the \sQuote{auto_clean} function in the tsaux package.
}
\item{trace}{
Whether to show the progress bar.
}
\item{\dots}{
Additional arguments passed to the \dQuote{auto_clean} function.
}
}
\value{
A list with the following data.tables:
\itemize{
\item sample : multihorizon white noise sample
\item sigma_scale: The scaling values per horizon which can be used to
scale the predictive distribution standard deviation.
\item sd_sigma_scale: The standard deviation of the sigma scaling factor
}
}
\note{
The walk forward predictive errors per horizon are used to estimate the realized
variance against the expected variance (sigma_scale). A jackknife procedure is
used for estimating the mean of this value. The output also includes samples from
a kernal fitted on the whitened errors per horizon which can then be used as
inputs to the prediction function.
}
\author{
Alexios Galanos
}