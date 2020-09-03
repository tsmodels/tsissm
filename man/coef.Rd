\name{coef.tsissm.estimate}
\alias{coef}
\alias{coef.tsissm.estimate}
\alias{logLik.tsissm.estimate}
\alias{AIC.tsissm.estimate}
\alias{logLik}
\title{
Estimated Model Extractor Functions
}
\description{
Coefficient, log-likelihood and AIC from an estimated tsissm model.
}
\usage{
\method{coef}{tsissm.estimate}(object, ...)
\method{logLik}{tsissm.estimate}(object, ...)
\method{AIC}{tsissm.estimate}(object, ..., k = 2)
}
\arguments{
\item{object}{
An object of class \dQuote{tsissm.estimate}
}
\item{k}{
(numeric) The penalty per parameter to be used;
the default k = 2 is the classical AIC.
}
\item{\dots}{
Not currently used.
}
}
\author{
Alexios Galanos
}
