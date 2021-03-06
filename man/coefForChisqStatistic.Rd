% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trunc_inf.R
\name{coefForChisqStatistic}
\alias{coefForChisqStatistic}
\title{Compute the coefficients related to truncated chi-squared test}
\usage{
coefForChisqStatistic(y, X, sigma, outlier.det, g)
}
\arguments{
\item{y, }{the response.}

\item{X, }{the design matrix (including intercepts).}

\item{sigma, }{the noise level \eqn{\sigma}.}

\item{outlier.det, }{indexes of outliers detected.}

\item{g, }{a subset of \eqn{[p]} indicating which group structure to test.}
}
\value{
This function returns a list of \code{(M, Mc, g, gc, P, chi, w, z, df)}.
    See the manuscript for details of those quantities.
}
\description{
This function computes the coefficients related to the truncated chi-squared test.
}
\keyword{internal}
