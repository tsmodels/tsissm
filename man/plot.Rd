% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot.tsissm.estimate}
\alias{plot.tsissm.estimate}
\alias{plot}
\alias{plot.tsissm.simulate}
\alias{plot.tsissm.profile}
\title{Object Plots}
\usage{
\method{plot}{tsissm.estimate}(x, y = NULL, ...)

\method{plot}{tsissm.simulate}(x, y = NULL, ...)

\method{plot}{tsissm.profile}(x, y = NULL, type = c("metrics", "coef"), ...)
}
\arguments{
\item{x}{an object of class \dQuote{tsissm.estimate}, \dQuote{tsissm.simulate}
or \dQuote{tsissm.profile}.}

\item{y}{not used.}

\item{...}{additional arguments passed to the underlying plot function.}

\item{type}{type of profile plot for objects of class \dQuote{tsissm.profile}.}
}
\description{
Plots for objects generated from the tsissm functions.
}
