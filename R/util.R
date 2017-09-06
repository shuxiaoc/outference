# ----- general purpose helper functions -----

# check the kkt condition of the lasso program
# min 1/2n ||y - X beta||_2^2 + lambda ||beta||_1
# return TRUE if the KKT condition is satisfied
#' Check the KKT condition of lasso solutions
#'
#' This function returns \code{TRUE} if and only if the lasso solution
#' satisfy the KKT condition.
#'
#' This function assumes the quadratic loss is scaled by \eqn{1/2n}. That is,
#'     \deqn{1/(2n) * ||y-X\beta||_2^2 + \lambda * ||\beta||_1}
#'
#' @keywords internal
#'
#' @param y, the response vector.
#' @param X, the design matrix.
#' @param lambda, the tuning parameter.
#' @param beta.hat, the lasso solution.
#' @param tol, the tolerance.
#'
#' @return This function returns the desired logical result.
checkKKT <- function(y, X, lambda, beta.hat, tol = 1e-5) {

  n <- length(y)
  p <- ncol(X)
  kkt <- 1/n * t(X) %*% (X %*% beta.hat - y)
  active <- which(beta.hat != 0)

  if (length(active) == 0) {
    # all beta.hat = 0
    # each entry of kkt should be between -lambda and lambda
    if (any(abs(kkt) > lambda + tol)) {
      return(FALSE)
    }
    return(TRUE)
  }

  # we know beta.hat is not all 0
  sign.active <- sign(beta.hat[active])
  # check kkt for active variables
  if (any(abs(kkt[active] + lambda*sign.active) > tol)) {
    return(FALSE)
  }


  if (length(active) < p) {
    # not all variables are active
    # check kkt for inactive variables
    if (any(abs(kkt[-active]) > lambda + tol)) {
      return(FALSE)
    }
  }
  return(TRUE)

}


#' Estimate the tuning parameter of lasso regression
#'
#' This function returns \eqn{\lambda = E[||X^T \epsilon||_\infty] / n},
#'     where \eqn{\epsilon ~ N(0, \sigma^2 I)}.
#'
#' @keywords internal
#'
#' @param X, the design matrix.
#' @param sigma, the noise level.
#' @param nsim, the number of Monte Carlo simulations to run.
#'
#' @return This function returns the value of tuning parameter.
estimateLambda <- function(X, sigma, nsim = 1000) {
  n <- nrow(X)
  eps <- matrix(stats::rnorm(n * nsim, 0, sigma), n, nsim)
  tepsX <- abs(crossprod(eps, X))
  imax <- max.col(tepsX) # max index of each row
  lambda = tepsX[cbind(1:nsim, imax)]
  return(mean(lambda)/n)
}


#' Estimate the noise level in linear regression with presence of outliers
#'
#' This function estimates the noise level \eqn{\sigma} in linear regression
#'     with possible presence of outliers, based on "lasso-refitting" strategy.
#'
#' Assume the mean-shift model \deqn{y = X \beta + u + \epsilon,} where
#'     \eqn{\epsilon ~ N(0, \sigma^2 I)}. This is equivalent to
#'     \deqn{y = X.enlarged \beta + \epsilon,} where \eqn{X.enlarged = (X : I_n)}.
#'     This function fits a lasso regression based on \eqn{(y, X.enlarged)}
#'     with the cross-validated tuning parameter,
#'     and then computes the residual sum of square, scaled by \eqn{1/(n-s)},
#'     where \eqn{s} is the number of active variables estimated by lasso.
#'
#' @keywords internal
#'
#' @param y, the response.
#' @param X, the design matrix with the first column being 1 (the column of intercepts).
#'
#' @return This function returns an estimate of the noise level \eqn{\sigma}.
#'
estimateSigma = function(y, X) {
  # y = X beta + u + eps is equivalent to
  # y = X.enlarged %*% c(beta, u) + eps
  n <- length(y)
  p <- ncol(X)
  X.enlarged <- cbind(X, diag(n))
  # we first choose lambda by cross-validation.
  # notice that we are NOT penalizing the intercept
  suppressWarnings(fit.lasso <- glmnet::cv.glmnet(x = X.enlarged[,-1], y = y, standardize = FALSE, intercept = T))
  lambda <- fit.lasso$lambda.min
  # utils::capture.output(
  #   # again notice that we are NOT penalizing the intercept
  #   fit.lasso <- penalized::penalized(response = y, penalized = X.enlarged[,-1], standardize = F,
  #                          lambda1 = lambda*n)
  # )
  fit.lasso <- glmnet::glmnet(x = X.enlarged[, -1], y = y, family = "gaussian", alpha = 1, lambda = lambda,
                              standardize = F, intercept = T)
  param.hat <- as.numeric(glmnet::coef.glmnet(fit.lasso))
  sigma <- sqrt(sum((y-X.enlarged%*%param.hat)^2)/(n-sum(param.hat!=0)))
  return(sigma)
}



#' Give two-sided p-values
#'
#' This function returns \code{2 * min(surv, 1-surv)}
#'
#' @keywords internal
#'
#' @param surv, usually the one-sided p-value.
#'
#' @return This function returns the desired two-sided p-value.
giveTwoSidePval <- function(surv){
  return(2 * min(surv, 1-surv))
}

# assume X has full column rank,
# return the projection matrix onto the column space of X.
# return 0 matrix is ncol(X) = 0
#' Projection matrix onto the column space
#'
#' This function assumes \eqn{X} has full column rank, and
#'     returns the projection matrix onto the column space of \eqn{X}.
#'
#' The projection matrix onto the column space of \eqn{X} can be written as
#'     \deqn{X(X^TX)^{-1}X^T.}
#'     When \eqn{X} is an empty matrix, then this function returns the zero matrix.
#'
#' @keywords internal
#'
#' @param X, a matrix with full column rank, can be an empty matrix with \code{ncol(X) == 0}.
#'
#' @return This function returns the desired square matrix with dimension \code{nrow(X)}.
proj <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) == 0) return(matrix(0, nrow = nrow(X), ncol = nrow(X)))
  return(X %*% MASS::ginv(X))
}





