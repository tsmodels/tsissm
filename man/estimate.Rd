\name{estimate}
\alias{estimate}
\alias{estimate.tsissm.spec}
\title{Model Estimation}
\description{
Estimates a model given a specification object using maximum likelihood.
}
\usage{
\method{estimate}{tsissm.spec}(object, solver = "nlminb",
control = list(trace = 0), autodiff = TRUE, ...)
}
\arguments{
\item{object}{
An object of class tsissm.spec.
}
\item{solver}{
One of either \dQuote{solnp}, \dQuote{nlminb} or \dQuote{optim}. The
latter uses the Nelder-Mead algorithm. For the nlminb algorithm,
an initial run of 200 iterations using the Nelder-Mead algorithm is
used prior to starting the main solver.
}
\item{control}{
Solver control parameters.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsissmad package.
}
\item{\dots}{
Not currently used.
}
}
\author{
Alexios Galanos
}