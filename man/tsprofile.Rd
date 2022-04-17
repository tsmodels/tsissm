\name{tsprofile}
\alias{tsprofile}
\alias{tsprofile.tsissm.estimate}
\title{
Model Simulation Based Profiling
}
\description{
Profiling of model dynamics using simulation/estimation/prediction.
}
\usage{
\method{tsprofile}{tsissm.estimate}(object, h = 1, nsim = 100,
seed = NULL, trace = FALSE, sigma_scale = 1,
solver = "nlminb", autodiff = TRUE, ...)
}
\arguments{
\item{object}{
An object of class tsissm.estimate
}
\item{h}{
The forecast horizon on which to evaluate performance metrics.
}
\item{nsim}{
The number of paths to generate.
}
\item{seed}{
An object specifying if and how the random number generator should be initialized.
See the simulate documentation for more details.
}
\item{trace}{
Whether to show the progress bar. The user is expected to have set up 
appropriate handlers for this using the \dQuote{progressr} package.
}
\item{sigma_scale}{
A scaling factor for the innovations sd.
}
\item{solver}{
Choice of solver to use for the estimation of the paths.
}
\item{autodiff}{
Whether to use automatic differentiation for estimation. This makes use
of the tsissmad package.
}
\item{\dots}{
Not currently used
}
}
\details{
The function profiles an estimated model by simulating and then estimating multiple paths from the
assumed DGP while leaving h values out for prediction evaluation. Each simulated path is equal to
the size of the original dataset plus h additional values, and initialized with the initial state
vector from the model. The resulting output contains the distribution of the MAPE, percent bias (BIAS)
and root mean squared log error (RMSLE) per horizon h. Since these matrices are of class
\dQuote{tsmodel.distribution} they can be readily plotted with the special purpose \dQuote{plot} function
for this class from the \dQuote{tsmethods} package. Additionally, a data.table matrix is also 
returned with the distribution of the coefficients from each path estimation.
}
\value{
See details above.
}
\note{
The function can use parallel functionality as long as the user has set up a
\code{\link[future]{plan}} using the future package.
}
\author{
Alexios Galanos
}