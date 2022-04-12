\name{tsdecompose}
\alias{tsdecompose}
\alias{tsdecompose.tsissm.estimate}
\alias{tsdecompose.tsissm.predict}
\title{
Model Decomposition
}
\description{
Decomposes the estimated model or prediction into its component parts.
}
\usage{
\method{tsdecompose}{tsissm.estimate}(object, simplify = FALSE, ...)
\method{tsdecompose}{tsissm.predict}(object, simplify = FALSE, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsissm.estimate} or \dQuote{tsissm.predict}.
}
\item{simplify}{
Simplification of the returned states aggregates the level and slope 
(if present) into a Trend component, all Seasonal components, all ARMA components 
and the error terms into an Irregular component. This may be useful when bagging 
is carried out (assuming equal lambda in the box-cox transform or the logit 
transform is used). This also simplifies the ability to created custom 
overwrites of the Trend and rebuilt the predictive distribution.
}
\item{\dots}{
Not currently used
}
}
\value{
For the estimated object, returns an xts matrix of the state components 
(including Irregular term). For the predicted object, a list with the simulated 
state components of class \dQuote{tsmodel.predict} which includes the 
predictive distribution and original (estimated) series components.
}
\details{
The 1-step ahead prediction is given by the following equation:

\deqn{y_{t} = x_{t-1} w + \varespsilon_{t}}

Because the decomposition prelags the states so that the seed state is aligned
with the error term, then summing the state distribution of each component with 
the returned error distribution will ensure that the exact same
predicted value matched to the correct date is returned.
}
\author{
Alexios Galanos
}
