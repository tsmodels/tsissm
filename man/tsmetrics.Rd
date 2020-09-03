\name{tsmetrics.tsissm.estimate}
\alias{tsmetrics}
\alias{tsmetrics.tsissm.estimate}
\alias{tsmetrics.tsissm.predict}
\title{
Performance Metrics
}
\description{
Performance metrics from an estimated or predicted tsissm model.
}
\usage{
\method{tsmetrics}{tsissm.estimate}(object, ...)
\method{tsmetrics}{tsissm.predict}(object, actual, alpha = 0.1, ...)
}
\arguments{
\item{object}{
An object of class \dQuote{tsissm.estimate}
}
\item{actual}{
The actual data matched to the dates of the forecasts.
}
\item{alpha}{
The coverage level for distributional forecast metrics.
}
\item{\dots}{
Not currently used.
}
}
\author{
Alexios Galanos
}
