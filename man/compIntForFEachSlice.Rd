% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trunc_inf.R
\name{compIntForFEachSlice}
\alias{compIntForFEachSlice}
\title{Solve the complement of one slice of the truncation set for the \eqn{F} statistic}
\usage{
compIntForFEachSlice(C, coeffs, tol = 1e-08, tol2 = 1e-06)
}
\arguments{
\item{C, coeffs, }{the coefficients related to the \eqn{F} statistic.}

\item{tol, tol2}{the tolerance of roots.}
}
\value{
This function returns an "Intervals" object.
}
\description{
This function computes one slice of the truncation set for the \eqn{F} statistic, then
    returns the complement of the truncation set (w.r.t. to \eqn{\mathbb{R}}).
}
\details{
This function is essentially a wrapper function for
    \code{selectiveInference:::TF_roots}, with minor changes so that
    the function does not stop when the constraint is infeasible (instead, this
    function returns an empty "Intervals" object).
    The reason we compute the complement instead of the original solution set
    is that when taking intersection of multiple sets, the strategy of
    "taking union of complements, then taking comlement" is substantially faster
    than taking intersections directly based on the \code{intervals} package.
}
\references{
Tibshirani, R., et al. "selectiveInference: Tools for Post-Selection Inference."
    R package version 1.3 (2016).
}
\keyword{internal}
