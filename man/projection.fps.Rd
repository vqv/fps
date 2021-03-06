% Generated by roxygen2 (4.0.1): do not edit by hand
\name{projection.fps}
\alias{projection.fps}
\title{Extract projection matrix estimates}
\usage{
\method{projection}{fps}(object, lambda, fixrank = FALSE, ...)
}
\arguments{
\item{object}{fps object}

\item{lambda}{lambda value to extract}

\item{fixrank}{should the rank of the projection matrix be fixed?}

\item{...}{other arguments}
}
\description{
Returns a fitted projection matrix
}
\examples{
data(wine)
out <- fps(cor(wine), ndim = 2)
projection(out, lambda = 0.5)
}

