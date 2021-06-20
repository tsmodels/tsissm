\name{tsbacktest}
\alias{tsbacktest}
\alias{tsbacktest.tsissm.spec}
\title{
Walk Forward Model Backtest
}
\description{
Generates an expanding window walk forward backtest.
}
\usage{
\method{tsbacktest}{tsissm.spec}(object, start = floor(length(object$target$y_orig)/2),
end = length(object$target$y_orig), h = 1, alpha = NULL,
cores = 1, data_name = "y", save_output = FALSE, save_dir = "~/tmp/",
solver = "optim", autodiff = FALSE, trace = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsissm.spec}.
}
\item{start}{
The numeric data index from which to start the backtest.
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
\item{alpha}{
An optional numeric vector of coverage rates for which to calculate
the quantiles.
}
\item{cores}{
The number of cores to use for parallel estimation. The package
makes use of the doSNOW package.
}
\item{data_name}{
The name of the dataset. Useful when combining mulitple results
together.
}
\item{save_output}{
Whether to save each estimated and predicted object.
}
\item{save_dir}{
The directory to save the estimated and predicted objects.
These are saved as rds files using \dQuote{model_estimationDate}
and \dQuote{predict_estimationDate} as names with estimationDate
the date on which the model was estimated.
}
\item{solver}{
The solver to use.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsissmad package.
}
\item{trace}{
Whether to show the progress bar.
}
\item{\dots}{
Not currently used
}
}
\value{
A list with the following data.tables:
\itemize{
\item prediction : the backtest table with forecasts and actuals
\item metrics: a summary performance table showing metrics by
forecast horizon (MAPE, MSLRE, BIAS and MIS if alpha was not NULL).
}
}
\author{
Alexios Galanos
}