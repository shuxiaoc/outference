% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outference.R
\name{coeftest}
\alias{coeftest}
\title{Test for a single regression coefficient from an \code{"outference"} object}
\usage{
coeftest(object, index)
}
\arguments{
\item{object, }{an object of class \code{"outference"}.}

\item{index, }{the index of parameter under consideration. The intercept is labelled as \code{index = 1}.}
}
\value{
If \code{sigma} is known or estimated, then this function returns a list with the following
    components:
    \item{Z, }{the value of the \eqn{Z}-statistic.}
    \item{sd, }{the standard error of the regression coefficient.}
    \item{truncation, }{the truncation set for the \eqn{Z}-statistic.}
    \item{pval, }{the corrected p-value.}
    If \code{sigma = NULL}, then this function returns a list with the following components:
    \item{truncation, }{the truncation set for the \eqn{F}-statistic.}
    \item{pval, }{the corrected p-value.}
}
\description{
This function does a two-sided selective test for a single regression coefficient
    (test whether it is zero or not) from an \code{"outference"} object.
}
\examples{
## Brownlee’s Stack Loss Plant Data
data("stackloss")
## fit the model
## detect outlier using Cook's distance with cutoff = 4
fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
coeftest(fit, 2)

}
\seealso{
\code{\link{grptest}} for testing group structures from an \code{"outference"} object.
}
\author{
Shuxiao Chen <sc2667@cornell.edu>
}
