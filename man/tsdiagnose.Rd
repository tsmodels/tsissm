\name{tsdiagnose}
\alias{tsdiagnose.tsissm.estimate}
\alias{tsdiagnose}
\title{
Model Diagnostics
}
\description{
Creates a short summary of model based diagnosticts.
}
\usage{
\method{tsdiagnose}{tsissm.estimate}(object, plot = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsissm.estimate}.
}
\item{plot}{
Whether to generate diagnostic plots to accompany summary.
}
\item{\dots}{
Not currently used.
}
}
\value{
A list of tables (printed out and returned insivibly) with
diagnostics for ARMA (if used) invertibility, forecastability (roots of D matrix),
Ljung-Box test for residual autocorrelation and outlier dates using the Rosner
Test from the \code{EnvStats} package.
}
\author{
Alexios Galanos.
}
