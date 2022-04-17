\name{tsmoments}
\alias{tsmoments}
\alias{tsmoments.tsissm.estimate}
\title{Analytic Forecast Moments}
\description{
Estimates a model given a specification object using maximum likelihood.
}
\usage{
\method{tsmoments}{tsissm.estimate}(object, h = 1, newxreg = NULL, 
init_states = NULL, transform = FALSE, ...)
}
\arguments{
\item{object}{
An object of class tsissm.estimate.
}
\item{h}{
The forecast horizon.
}
\item{newxreg}{
The matrix of external regressors in the forecast horizon.
}
\item{init_states}{
An optional vector of states to initialize the forecast. If NULL,
will use the last available state from the estimated model.
}
\item{transform}{
Whether to backtransform the mean forecast. For the Box-Cox and logit 
transformations this uses a Taylor Series expansion to adjust for the
bias.
}
\item{\dots}{
Not currently used.
}
}
\author{
Alexios Galanos
}