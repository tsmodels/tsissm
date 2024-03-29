% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{logLik.tsissm.estimate}
\alias{logLik.tsissm.estimate}
\alias{logLik}
\title{Model Log-Likelihood}
\usage{
\method{logLik}{tsissm.estimate}(object, ...)
}
\arguments{
\item{object}{an object of class \dQuote{tsissm.estimate}.}

\item{...}{not currently used.}
}
\value{
Returns an object of class logLik. This is a number with at least one
attribute, "df" (degrees of freedom), giving the number of (estimated)
parameters in the model.
}
\description{
Extract the log-likelihood from an estimated model.
}
