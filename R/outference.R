#' Fit a linear model with outliers detected and removed
#'
#' This function detects outliers using a user-specified method, and
#'     fits a linear regression model with outliers removed. The object
#'     returned by this function can be used for valid inference corrected
#'     for outlier removal through generic functions like \code{\link{summary}},
#'     \code{\link{confint}}, \code{\link{predict}}.
#'
#' This function uses the same syntax as \code{\link{lm}} for the \code{formula} and \code{data} arguments.
#'     Users can access the original \code{"lm"} objects through \code{$fit.full} and \code{$fit.rm}.
#'     Common generic functions for \code{lm}, including \code{\link{coef}}, \code{\link{confint}},
#'     \code{\link{plot}}, \code{\link{predict}} and \code{\link{summary}} are re-written so that
#'     they can be used to extract useful features of the object returned by this function.
#'
#'     Currently, this function supports three outlier detection methods. For \code{"cook"}, the \eqn{i}-th
#'     observation is considered as an outlier when its Cook's distance is greater than \code{cutoff/n},
#'     where \code{n} is the number of observations. For \code{"dffits"}, the \eqn{i}-th observation is
#'     considered as an outlier when the square of its DFFITS measure is greater than \code{cutoff*p/(n-p)},
#'     where \code{p} is the number of variables (including the intercept). The rule of thumb of \code{cutoff}
#'     for both methods are 4, which is the default value when one sets \code{cutoff = NULL}.
#'     The outlier detection event of both methods can be characterized as a set of quadratic constraints
#'     in the response \eqn{y}:
#'     \deqn{\bigcap_{i \in [n]} {y^T Q_i y \ge 0},}
#'     and the constraint returned by this function is the list of \eqn{Q_i} matrices.
#'     For \code{"lasso"}, we assume the \emph{mean-shift model}
#'     \eqn{y = X \beta + u + \epsilon}, where \eqn{u} is the "outlying coefficients" and
#'     \eqn{\epsilon ~ N(0, \sigma^2 I)} is the noise. We solve the following program:
#'     \deqn{(\hat \beta, \hat u) = argmin ||y-X\beta-u||_2^2 + cutoff*||u||_1.} The \eqn{i}-th observation
#'     is considered as an outlier when \eqn{\hat u_i} differs from \eqn{0}. The default cutoff for
#'     \code{"lasso"} is \eqn{0.75*E[||X^T \epsilon||_\infty]/n}, which is a less conservative choice
#'     than the prediction-optimal cutoff \eqn{2*E[||X^T \epsilon||_\infty]/n}. This cutoff is computed
#'     by Monte Carlo simulation and \eqn{\sigma} is replaced by an estimate when the true noise level
#'     is unknown. The outlier detection event of \code{"lasso"} can be characterized as a
#'     set of affine constraints in the response \eqn{y}:
#'     \deqn{A y \ge b, }
#'     where the \eqn{"\ge"} is interpreted as element-wise. The constraint returned by this function is
#'     then a list of \code{(A, b)}.
#'
#' @export
#'
#' @aliases print.outference
#'
#' @param formula, an object of class \code{"\link{formula}"}, the same syntax as in \code{\link{lm}}.
#' @param data, an optional data frame, list or environment containing the variables in the model, the same
#'     syntax as in \code{\link{lm}}.
#' @param method, the outlier detection method, must be one of \code{"cook", "dffits", "lasso"}. See also 'details'.
#' @param cutoff, the cutoff of the outlier detection method. If \code{cutoff = NULL}, then this function
#'     uses the default values. For \code{"cook"} or \code{"dffits"}, the default cutoff is \eqn{4};
#'     for \code{"lasso"}, the default cutoff is \eqn{0.75*E[||X^T \epsilon||_\infty]/n}.
#'     See also 'details'.
#' @param sigma, the noise level. Must be one of \code{NULL, "estimate"}, or a positive scaler value. If
#'     \code{sigma = NULL}, then the inference will assume the noise level is unknown;
#'     if \code{sigma = "estimate"}, then the inference will base on an estimated noise level.
#' @param x, an object of class \code{"outference"}.
#' @param digits, the number of significant digits to use when printing.
#' @param ..., other arguments.
#'
#' @return This function returns an object of \code{\link{class}} \code{"outference"}.
#'
#'     The function \code{\link{summary}} is used to obtain and print a summary (including p-values)
#'     of the results. The generic functions \code{\link{coef}}, \code{\link{confint}}, \code{\link{plot}},
#'     \code{\link{predict}} are used to extract useful features of the object returned by this function.
#'
#'     An object of class \code{"outference"} is a list containing the following components:
#'
#'
#'     \item{fit.full, }{an \code{"lm"} object representing the fit using the full data (no outliers are removed).}
#'     \item{fit.rm, }{an \code{"lm"} object representing the fit using the data after outlier removal.}
#'     \item{method, }{the method used for outlier detection.}
#'     \item{cutoff, }{the cutoff of the method.}
#'     \item{outlier.det, }{indexes of detected outliers.}
#'     \item{magnitude, }{a measure of "outlying-ness". For \code{"cook"} and \code{"dffits"}, this is
#'         the vector of the Cook's distance or DFFITS for all observations; for \code{"lasso"}, this is
#'         the vector of "outlying coefficients" estimated by lasso. See also 'details'.}
#'     \item{constraint, }{the constraint in the response that characterizes the outlier detection event.
#'         For \code{"cook"} and \code{"dffits"}, this is a list of \eqn{n} by \eqn{n} matrices;
#'         for \code{"lasso"}, this is a list of \code{(A, b)}, where \code{A} is a matrix and
#'         \code{b} is a vector. See also 'details'.}
#'     \item{sigma, }{the noise level used in the fit.}
#'     \item{call, }{the function call.}
#'
#' @seealso \code{\link{summary.outference}} for summaries;
#'
#'     \code{\link{coef.outference}} for extracting coefficients;
#'
#'     \code{\link{confint.outference}} for confidence intervals of regression coefficients;
#'
#'     \code{\link{plot.outference}} for plotting the outlying measure;
#'
#'     \code{\link{predict.outference}} for making predictions.
#'
#' @references Lee, Jason D., et al. "Exact post-selection inference, with application to the lasso."
#'     The Annals of Statistics 44.3 (2016): 907-927.
#' @references S. Chen and J. Bien. “Valid Inference Corrected for Outlier Removal”. arXiv preprint arXiv:1711.10635 (2017).
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' head("stackloss")     # look at the dataset
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' plot(fit)             # plot the Cook's distance of each observation
#' ## observation 21 is considered as an outlier with cutoff = 4
#' summary(fit$fit.full) # look at the fit with all the data
#' summary(fit$fit.rm)   # look at the fit with observation 21 deleted
#' summary(fit)          # extract the corrected p-values after outlier removal
outference <- function(formula, data, method = c("cook", "dffits", "lasso"),
                      cutoff = NULL, sigma = NULL) {
  this.call <- match.call()
  method <- match.arg(method)
  lm.call.index<- match(c("formula", "data"), names(this.call), 0)
  lm.call <- this.call[c(1, lm.call.index)]
  lm.call[[1]] <- quote(stats::lm)

  fit.full <- eval(lm.call, parent.frame())
  y <- stats::model.response(fit.full$model)
  X <- stats::model.matrix(fit.full)
  n <- length(y)
  p <- ncol(X)


  if (!is.null(sigma)) {
    if (sigma == "estimate") sigma = estimateSigma(y = y, X = X)
    else if (!(is.numeric(sigma) && length(sigma) == 1 && sigma > 0))
      stop("sigma can only be NULL, \"estimate\" or a user-specified positive scaler")
  }

  if (p == 0) stop("empty model")
  # some constants
  XtX <- crossprod(X)
  Xplus <- chol2inv(chol(XtX)) %*% t(X)
  PX <- X %*% Xplus
  PXperp <- diag(n) - PX


  # affine constraints
  if (method == "lasso") {
    PXperpY <- PXperp %*% y
    # we set cutoff =  constant * E[||PXperp %*% eps||_\infty]
    if (is.null(sigma)) {
      sigma.lasso <- estimateSigma(y, X)
    }
    else {
      sigma.lasso <- sigma
    }

    if (is.null(cutoff)) {
      # using a less conservative choice of 0.75 rather than 2
      cutoff <- 0.75*estimateLambda(PXperp, sigma.lasso)
    }
    else if (is.character(cutoff)) {
      times <- as.numeric(substr(cutoff, start = 1, stop = nchar(cutoff)-1))
      cutoff <- times * estimateLambda(PXperp, sigma.lasso)
    }
    else if (is.numeric(cutoff) && (length(cutoff) == 1)) {
      cutoff <- as.numeric(cutoff)
    }
    else {
      stop("for lasso, cutoff must be null, a scaler, or a string like \'0.5x\'")
    }

    # # the calculation of polyhedron is highly sensitive to the lasso solution,
    # # so we use the penalized package for high-precision calculation of lasso solutions.
    # # our preliminary simulation shows that the solution by glmnet occationally, though very rare,
    # # is not accurate enough for polyhedron calculations.
    # utils::capture.output(
    #   fit.lasso <- penalized::penalized(response = PXperpY, penalized = PXperp, unpenalized = ~0, standardize = F,
    #                          lambda1 = cutoff*n)
    # )
    # u.hat = penalized::coef(fit.lasso, "all")
    fit.lasso <- glmnet::glmnet(x = PXperp, y = PXperpY, family = "gaussian", alpha = 1,
                                lambda = cutoff, standardize = F, intercept = F,
                                thresh = 1e-10, maxit = 1e7)

    u.hat <- as.numeric(glmnet::coef.glmnet(fit.lasso)[-1])

    # ad-hoc check of kkt condition to make sure the lasso problem is what we actually want
    if (checkKKT(y = PXperpY, X = PXperp, lambda = cutoff, beta.hat = u.hat) == FALSE) {
      warning("this lasso solution does not satisfy kkt conditions!")
    }

    outlier.det <- which(u.hat != 0)
    if (n - length(outlier.det) <= p) stop("number of remaining observations less than number of variables, the model is singular")
    if (length(outlier.det) == 0) fit.rm <- fit.full
    else { # at least one outlier detected
      lm.call$subset <- -outlier.det
      fit.rm <- eval(lm.call, parent.frame())
    }

    constr <- constrInResponseLasso(n, p, PXperp, outlier.det, sign(u.hat[outlier.det]), cutoff)

    # ad-hoc check of the polyhedron: Ay >= b should always hold
    if (any(constr$A %*% y - constr$b < -1e-5)) {
      stop("constraint for lasso is problematic")
    }

    out <- list(fit.full = fit.full, fit.rm = fit.rm, method = method, cutoff = cutoff, outlier.det = outlier.det,
               magnitude = u.hat, constraint = constr, sigma = sigma, call = this.call)
    class(out) <- "outference"
    return(out)
  }

  # quadratic constraints
  if (method == "cook") {

    if (is.null(cutoff)) cutoff <- 4
    magnitude <- as.numeric(stats::cooks.distance(fit.full))
    outlier.det <- which(magnitude >= cutoff/n)
    if (n - length(outlier.det) <= p) stop("number of remaining observations less than number of variables, the model is singular")
    if (length(outlier.det) == 0) fit.rm <- fit.full
    else { # at least one outlier detected
      lm.call$subset <- -outlier.det
      fit.rm <- eval(lm.call, parent.frame())
    }
    constr <- constrInResponseCook(n, p, PX, PXperp, outlier.det, cutoff)
    out <- list(fit.full = fit.full, fit.rm = fit.rm, method = method, cutoff = cutoff, outlier.det = outlier.det,
               magnitude = magnitude, constraint = constr, sigma = sigma, call = this.call)
    class(out) <- "outference"
    return(out)
  }


  # quadratic constraints
  if (method == "dffits") {
    if (is.null(cutoff)) cutoff <- 4
    magnitude <- as.numeric(stats::dffits(fit.full))
    outlier.det <- which((magnitude)^2 >= cutoff * p / (n-p))
    if (n - length(outlier.det) <= p) stop("number of remaining observations less than number of variables, the model is singular")
    if (length(outlier.det) == 0) fit.rm <- fit.full
    else { # at least one outlier detected
      lm.call$subset <- -outlier.det
      fit.rm <- eval(lm.call, parent.frame())
    }

    constr <- constrInResponseDffits(n, p, PX, PXperp, outlier.det, cutoff)
    out <- list(fit.full = fit.full, fit.rm = fit.rm, method = method, cutoff = cutoff, outlier.det = outlier.det,
               magnitude = magnitude, constraint = constr, sigma = sigma, call = this.call)
    class(out) = "outference"
    return(out)
  }
}



#' Fit a linear model with outliers detected SEQUENTIALLY
#'
#' This function detects outliers by using Cook's distance sequentially, and
#'     fits a linear regression model with outliers removed. The object
#'     returned by this function can be used for valid inference corrected
#'     for outlier removal through generic functions like \code{\link{summary}},
#'     \code{\link{confint}}, \code{\link{predict}}.
#'
#' This function uses the same syntax as \code{\link{lm}} for the \code{formula} and \code{data} arguments.
#'     Users can access the original \code{"lm"} objects through \code{$fit.full} and \code{$fit.rm}.
#'     Common generic functions for \code{lm}, including \code{\link{coef}}, \code{\link{confint}},
#'     \code{\link{plot}}, \code{\link{predict}} and \code{\link{summary}} are re-written so that
#'     they can be used to extract useful features of the object returned by this function.
#'
#'     The \eqn{i}-th observation is considered as an outlier when its Cook's distance
#'     rank among top \eqn{k}, where \eqn{k} is the user-specified number of outliers to be detected.
#'     The outlier detection event can be characterized as a set of quadratic constraints
#'     in the response \eqn{y}:
#'     \deqn{\bigcap_{i \in I} {y^T Q_i y \ge 0},}
#'     where \eqn{I} is a finite index set, and the constraint returned by this function
#'     is the list of \eqn{Q_i} matrices.
#'
#' @keywords internal
#'
#' @param formula, an object of class \code{"\link{formula}"}, the same syntax as in \code{\link{lm}}.
#' @param data, an optional data frame, list or environment containing the variables in the model, the same
#'     syntax as in \code{\link{lm}}.
#' @param sigma, the noise level. Must be one of \code{NULL, "estimate"}, or a positive scaler value. If
#'     \code{sigma = NULL}, then the inference will assume the noise level is unknown;
#'     if \code{sigma = "estimate"}, then the inference will base on an estimated noise level.
#' @param numOfOutlier, the number of outliers to be detected.
#'
#'
#' @return This function returns an object of \code{\link{class}} \code{c("outference_seq", "outference")}.
#'
#'     The function \code{\link{summary}} is used to obtain and print a summary (including p-values)
#'     of the results. The generic functions \code{\link{coef}}, \code{\link{confint}}, \code{\link{plot}},
#'     \code{\link{predict}} are used to extract useful features of the object returned by this function.
#'
#'     An object of class \code{c("outference_seq", "outference")} is a list containing the following components:
#'
#'     \item{fit.full, }{an \code{"lm"} object representing the fit using the full data (no outliers are removed).}
#'     \item{fit.rm, }{an \code{"lm"} object representing the fit using the data after outlier removal.}
#'     \item{method, }{"cook".}
#'     \item{cutoff, }{\code{NULL}.}
#'     \item{numOfOutlier, }{the number of outliers to be detected.}
#'     \item{outlier.det, }{indexes of detected outliers.}
#'     \item{magnitude, }{ the vector of the Cook's distance for all observations}
#'     \item{constraint, }{the constraint in the response that characterizes the outlier detection event.
#'         A list of \eqn{n} by \eqn{n} matrices.}
#'     \item{sigma, }{the noise level used in the fit.}
#'     \item{call, }{the function call.}
#'
#' @seealso \code{\link{summary.outference}} for summaries;
#'
#'     \code{\link{coef.outference}} for extracting coefficients;
#'
#'     \code{\link{confint.outference}} for confidence intervals of regression coefficients;
#'
#'     \code{\link{plot.outference}} for plotting the outlying measure;
#'
#'     \code{\link{predict.outference}} for making predictions.
#'
#'
#' @references S. Chen and J. Bien. “Valid Inference Corrected for Outlier Removal”. arXiv preprint arXiv:1711.10635 (2017).
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
outference_seq <- function(formula, data, sigma = NULL, numOfOutlier) {

  this.call <- match.call()
  lm.call.index <- match(c("formula", "data"), names(this.call), 0)
  lm.call <- this.call[c(1, lm.call.index)]
  lm.call[[1]] <- quote(stats::lm)

  fit.full <- eval(lm.call, parent.frame())
  y <- stats::model.response(fit.full$model)
  X <- stats::model.matrix(fit.full)
  n <- length(y)
  p <- ncol(X)


  if (!is.null(sigma)) {
    if (sigma == "estimate") sigma = estimateSigma(y = y, X = X)
    else if (is.numeric(sigma) && length(sigma) == 1 && sigma > 0) sigma = sigma
    else stop("sigma can only be NULL, \"estimate\" or a user-specified positive scaler")
  }

  if (p == 0) stop("empty model")
  # some constants
  XtX <- crossprod(X)
  Xplus <- chol2inv(chol(XtX)) %*% t(X)
  PX <- X %*% Xplus
  PXperp <- diag(n) - PX

  # quadratic constraints
  magnitude <- as.numeric(stats::cooks.distance(fit.full))
  obs.ordered <- order(magnitude, decreasing = T)
  outlier.det <- obs.ordered[1:numOfOutlier]
  if (n - length(outlier.det) <= p) stop("number of remaining observations less than number of variables, the model is singular")
  if (length(outlier.det) == 0) fit.rm <- fit.full
  else { # at least one outlier detected
    lm.call$subset <- -outlier.det
    fit.rm <- eval(lm.call, parent.frame())
  }
  constr <- constrInResponseCookSeq(n = n, p = p, PX = PX, PXperp = PXperp, obs.ordered = obs.ordered, numOfOutlier = numOfOutlier)
  out <- list(fit.full = fit.full, fit.rm = fit.rm, method = "cook", cutoff = NULL, numOfOutlier = numOfOutlier, outlier.det = outlier.det,
              magnitude = magnitude, constraint = constr, sigma = sigma, call = this.call)
  class(out) <- c("outference_seq", "outference")
  return(out)

}


#' @rdname outference
#' @export
print.outference <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  # first we essentially do print.lm(x$fit.rm)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Outlier detection method: ", x$method, " with cutoff = ", x$cutoff, sep = "")
  if (length(x$outlier.det) == 0) {
    cat("\nOutlier detected: none")
  }
  else {
    cat("\nOutlier detected:", x$outlier.det)
  }
  cat("\n\nCoefficients:\n")
  print.default(format(stats::coef(x$fit.rm), digits = digits), print.gap = 2, quote = FALSE)
  #cat("")
  invisible(x)
}


#' Plot the outlying measure from an \code{"outference"} object
#'
#' This function plots the Cook's distance, or DFFITS, or the estimated
#' outlying coefficients, as well as the cutoff based on the outlier detection method.
#'
#' @export
#'
#' @param x, an object of class \code{"outference"}.
#' @param ..., other arguments.
#'
#' @seealso \code{\link{outference}} for model fitting;
#'
#'     \code{\link{summary.outference}} for summaries;
#'
#'     \code{\link{coef.outference}} for extracting coefficients;
#'
#'     \code{\link{confint.outference}} for confidence intervals of regression coefficients;
#'
#'     \code{\link{predict.outference}} for making predictions.
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' plot(fit)             # plot the Cook's distance of each observation
#'
plot.outference <- function(x, ...) {
  y <- stats::model.response(x$fit.full$model)
  X <- stats::model.matrix(x$fit.full)
  n <- length(y)
  p <- ncol(X)

  method <- x$method
  cutoff <- x$cutoff
  magnitude <- x$magnitude
  if (method == "lasso") {
    cutoff <- 0
    ylab.info <- "Estimated magnitude of outliers"
  }
  else if (method == "cook") {
    cutoff <- cutoff/n
    ylab.info <- "Cook's distance"
  }
  else if (method == "dffits") {
    cutoff <- sqrt(cutoff * p / (n-p))
    magnitude <- abs(magnitude)
    ylab.info <- "|DFFITS|"
  }
  else {
    stop("method should be one of \"cook\", \"dffits\", \"lasso\"")
  }

  graphics::plot(magnitude, ylab = ylab.info)
  graphics::abline(h = cutoff, lty = 2)
}



#' Summarize from an \code{"outference"} object
#'
#' This function produces a summary from an \code{"outference"} object, with a similar fasion
#'     as \code{\link{summary.lm}}.
#'
#' This function is written in a similar fasion as \code{\link{summary.lm}}. Users can get access
#'     to the \code{"summary.lm"} objects through \code{$summary.full} and \code{$summary.rm}.
#'
#' @export
#'
#' @aliases print.summary.outference
#'
#' @param object, an object of class \code{"outference"}.
#' @param x, an object of class \code{"summary.outference"}.
#' @param digits, the number of significant digits to use when printing.
#' @param signif.stars, should the 'significance starts' be printed?
#' @param ..., other arguments.
#'
#' @return This function returns an object of class \code{"summary.outference"}, which is a list containing
#'     the following components:
#'     \item{call, }{the function call.}
#'     \item{summary.full, }{an object of class \code{"summary.lm"}, representing the summary
#'         from the fit using the full data.}
#'     \item{summary.rm, }{an object of class \code{"summary.lm"}, representing the summary
#'         from the fit after outlier removal.}
#'     \item{method, }{the method used for outlier detection.}
#'     \item{cutoff, }{the cutoff of the method.}
#'     \item{outlier.det, }{indexes of detected outliers.}
#'     \item{magnitude, }{a measure of "outlying-ness". For \code{"cook"} and \code{"dffits"}, this is
#'         the vector of the Cook's distance or DFFITS for all observations; for \code{"lasso"}, this is
#'         the vector of "outlying coefficients" estimated by lasso.}
#'     \item{sigma, }{the noise level used in the fit.}
#'     \item{coefficients, }{a data frame summarizing the estimates, standard errors, values of the test
#'         statistics and corrected p-values of regression coefficients.}
#'     \item{truncation.coef, }{a list of the truncation sets of the test statistics for each regression
#'         coefficient.}
#'     \item{chisqstatistic, fstatistic, }{a list containing the value, the degree(s) of freedom, the
#'         truncation set, and the corrected p-value for testing the global null.}
#'
#' @seealso \code{\link{outference}} for model fitting;
#'
#'     \code{\link{coef.outference}} for extracting coefficients;
#'
#'     \code{\link{confint.outference}} for confidence intervals of regression coefficients;
#'
#'     \code{\link{plot.outference}} for plotting the outlying measure;
#'
#'     \code{\link{predict.outference}} for making predictions.
#'
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' summary(fit)          # extract the corrected p-values after outlier removal
summary.outference <- function(object, ...) {
  out <- list()
  call <- object$call
  summary.full <- stats::summary.lm(object$fit.full)
  summary.rm <- stats::summary.lm(object$fit.rm)

  out$call <- call
  out$summary.full <- summary.full
  out$summary.rm <- summary.rm
  out$method <- object$method
  outlier.det <- object$outlier.det
  out$outlier.det <- outlier.det
  out$cutoff <- object$cutoff
  out$magnitude <- object$magnitude
  sigma <- object$sigma
  out$sigma <- sigma
  # extract stuff to be used for inference
  y <- stats::model.response(object$fit.full$model)
  X <- stats::model.matrix(object$fit.full)
  n <- length(y)
  p <- ncol(X)


  test.global <- (p > 1)
  if (attr(X, "assign")[1] == 0) { # we know intercept is included automatically
    test.global <- (p > 1)
    g.global <- 2:p # if p = 1, then test.global = FALSE, and g.global will never be used
  }
  else { # we know the user does not want to include intercept
    test.global <- (p >= 1)
    g.global <- 1:p
  }

  if (is.null(sigma)) { # do selective t and f test

    # inference for each beta^M_j
    temp <- lapply(1:p, function(index) {
      return(coeftest(object, index))
    })

    truncation.coef <- lapply(temp, function(temp.each) {
      temp.each$truncation
    })
    p.coef <- vapply(temp, function(temp.each) {
      temp.each$pval
    }, FUN.VALUE = 0.1)
    # format the result
    coef.rm <- summary.rm$coef
    coef.rm[, 4] <- p.coef
    colnames(coef.rm)[4] <- "Corrected p-value"
    out$coefficients <- coef.rm
    out$truncation.coef <- truncation.coef

    # do the global test
    if (test.global) {
      temp <- grptest(object, g.global)
      out$fstatistic <- list(values = summary.rm$fstatistic, truncation.global = temp$truncation, p.global = temp$pval)
    }
  }
  else { # do selective z and chisq test
    temp <- lapply(1:p, function(index) {
      coeftest(object, index)
    })
    truncation.coef <- lapply(temp, function(temp.each) {
      return(temp.each$truncation)
    })
    # the Z-statistic
    Z.coef <- vapply(temp, function(temp.each) {
      return(temp.each$Z)
    }, FUN.VALUE = 0.1)
    # standard error for Z statistic
    se.coef <- vapply(temp, function(temp.each) {
      return(temp.each$se)
    }, FUN.VALUE = 0.1)
    # p-value
    p.coef <- vapply(temp, function(temp.each) {
      return(temp.each$pval)
    }, FUN.VALUE = 0.1)
    coef.rm <- summary.rm$coef
    coef.rm[, -1] <- cbind(se.coef, Z.coef, p.coef)
    colnames(coef.rm)[3:4] <- c("z value", "Corrected p-value")
    out$coefficients <- coef.rm
    out$truncation.coef <- truncation.coef

    # do the global test
    if (test.global) {
      temp <- grptest(object, g.global)
      chisqstatistic <- c(temp$statistic, temp$df)
      names(chisqstatistic) <- c("value", "df")
      out$chisqstatistic <- list(values = chisqstatistic, truncation.global = temp$truncation, p.global = temp$pval)
    }


  }
  class(out) <- "summary.outference"
  return(out)
}

#' @rdname summary.outference
#' @export
print.summary.outference <- function(x, digits = max(3, getOption("digits") - 3),
                                     signif.stars = getOption("show.signif.stars"), ...) {

  cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
  cat("Outlier detection method: ", x$method, " with cutoff = ", x$cutoff, sep = "")

  if (length(x$outlier.det) == 0) {
    cat("\nOutlier detected: none")
  }
  else {
    cat("\nOutlier detected:", x$outlier.det)
  }

  cat("\n\n")
  r <- x$summary.rm$residuals
  r <- structure(zapsmall(stats::quantile(r), digits = digits + 1), names = c("Min", "1Q", "Median", "3Q", "Max"))
  print(r)
  cat("\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE)
  cat("\nResidual standard error:", format(signif(x$summary.rm$sigma, digits)), "on", x$summary.rm$df[2], "degrees of freedom")
  cat("\nMultiple R-squared: ", formatC(x$summary.rm$r.squared, digits = digits), "   Adjusted R-squared: ",
      format(x$summary.rm$adj.r.squared, digits = digits))
  if (!is.null(x$fstatistic)) {
    temp <- x$fstatistic
    cat("\nF-statistic: ", formatC(temp$values[1], digits = digits), " on ", temp$values[2],  " and ", temp$values[3],
        " DF,  Corrected p-value: ", format.pval(temp$p.global, digits = digits), sep = '')
  }
  else if (!is.null(x$chisqstatistic)) {
    temp = x$chisqstatistic
    cat("\nChisq-statistic: ", formatC(temp$values[1], digits = digits), " on ", temp$values[2]," DF,  Corrected p-value: ",
        format.pval(temp$p.global, digits = digits), sep = '')
  }
  cat("\n")
  invisible(x)
}



#' Extract coefficients from an \code{"outference"} object
#'
#' This function extracts the estimated regression coefficients from an \code{"outference"} object.
#'
#' @export
#'
#' @param object, an object of class \code{"outference"}.
#' @param model, if \code{model = "removed"}, then the coefficients are from the fit with outliers removed;
#'     if \code{model = "original"}, then the coefficients are from the fit with all the data.
#' @param ..., other arguments.
#'
#' @return This function returns the desired vector of regression coefficients.
#'
#' @seealso \code{\link{outference}} for model fitting;
#'
#'     \code{\link{summary.outference}} for summaries;
#'
#'     \code{\link{confint.outference}} for confidence intervals of regression coefficients;
#'
#'     \code{\link{plot.outference}} for plotting the outlying measure;
#'
#'     \code{\link{predict.outference}} for making predictions.
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' coef(fit, model = "original")    # the coefficients from the fit using all the data
#' coef(fit, model = "removed")     # the coefficients from the fit after outlier removal
coef.outference <- function(object, model = c("removed", "original"), ...) {
  arg <- match.arg(model)
  if (arg == "removed") return(stats::coefficients(object$fit.rm))
  return(stats::coefficients(object$fit.full))
}



#' Make predictions from an \code{"outference"} object
#'
#' This function gives predictions as well as confidence intervals for the regression surfaces
#'     and the prediction intervals from an \code{"outference"} object. The syntax is the same
#'     as \code{\link{predict.lm}}.
#'
#' If \code{alpha.tilde = NULL}, then this function iterates over
#'     \code{alpha.tilde in seq(0, 1-level, length.out = 100)[c(5, 25, 50, 75, 95)]} and returns
#'     the results with the shortest prediction intervals.
#'
#' @export
#'
#' @param object an object of class \code{"outference"}.
#' @param newdata, an optional data frame in which to look for variables with which to predict. If
#'     omitted, the fitted values are used. WARNING: making predictions for many new data points with
#'     \code{interval = "confidence"} or \code{"prediction"} can be time-consuming,
#'     since for each data point, the function needs to compute the truncation set and
#'     solve the roots for a truncated survival function.
#' @param interval, type of interval calculation. If set to \code{"none"}, then only
#'     point predictions are made; if set to \code{"confidence"}, then this function returns
#'     confidence intervals for the regression surface; if set to \code{"prediction"}, then this
#'     function returns the prediction intervals.
#' @param level, confidence level, default to \eqn{0.95}.
#' @param alpha.tilde, an extra parameter between \code{0} and \code{1-level}, which is  used when
#'     computing prediction intervals. If left \code{NULL}, then this function searches the
#'     \code{alpha.tilde} that gives shortest prediction intervals. See also 'details'.
#' @param ..., other arguments.
#'
#' @return This function gives a vector of predictions or a matrix of predictions and intervals
#'     with column names \code{fit, lwr, upr} if \code{interval} is set.
#'
#'
#' @seealso \code{\link{outference}} for model fitting;
#'
#'     \code{\link{summary.outference}} for summaries;
#'
#'     \code{\link{coef.outference}} for extracting coefficients;
#'
#'     \code{\link{confint.outference}} for confidence intervals of regression coefficients;
#'
#'     \code{\link{plot.outference}} for plotting the outlying measure;
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' ## predictions at the first two observations
#' predict(fit, newdata = stackloss[1:2, ], interval = "none")
#' predict(fit, newdata = stackloss[1:2, ], interval = "confidence")
#' predict(fit, newdata = stackloss[1:2, ], interval = "prediction")
#'

predict.outference <- function(object, newdata, interval = c("none", "confidence", "prediction"),
                               level = 0.95, alpha.tilde = NULL, ...) {

  interval <- match.arg(interval)
  if (interval == "none") return(stats::predict.lm(object$fit.rm, newdata = newdata, level = level))
  # we need to induce selective inference
  alpha <- 1-level
  # extract stuff to be used for inference
  y <- stats::model.response(object$fit.full$model)
  X <- stats::model.matrix(object$fit.full)
  n <- length(y)
  p <- ncol(X)
  outlier.det <- object$outlier.det
  sigma <- object$sigma
  if (is.null(sigma)) sigma <- estimateSigma(y, X)

  if(missing(newdata) || is.null(newdata)) {
    print("Warning: form intervals for all fitted values can be slow!")
    X.new <- stats::model.matrix(object$fit.rm)
  }
  else {
    X.new <- stats::model.matrix(stats::delete.response(stats::terms(object$fit.full)), data = newdata)
  }

  if (interval == "confidence") {

    out <- stats::predict.lm(object$fit.rm, newdata = newdata, interval = "confidence", level = level)
    int <- vapply(1:nrow(X.new), function(i) {
      x0 <- X.new[i, ]
      v <- contrastForSurf(x0 = x0, X = X, outlier.det = outlier.det)
      z = y - v * sum(v*y) / sum(v*v)
      truncation <- intForZAll(method = object$method, constraint = object$constraint,
                              v = v, z = z, sigma = sigma, outlier.det = outlier.det)
      return(computeCI(v, y, sigma, truncation, alpha))

    }, FUN.VALUE = c(0.1, 0.2))
    out[, 2:3] <- t(int)
    return(out)
  }
  else { # interval = "prediction"

    out <- stats::predict.lm(object$fit.rm, newdata = newdata, interval = "prediction", level = level)

    if (is.null(alpha.tilde)) alpha.tilde = seq(0, alpha, length.out = 100)[c(5, 25, 50, 75, 95)]

    int <- vapply(1:nrow(X.new), function(i) {
      x0 <- X.new[i, ]
      v <- contrastForSurf(x0 = x0, X = X, outlier.det = outlier.det)
      z <- y - v * sum(v*y) / sum(v*v)
      truncation <- intForZAll(method = object$method, constraint = object$constraint,
                              v = v, z = z, sigma = sigma, outlier.det = outlier.det)

      int.candidate <- vapply(alpha.tilde, function(alpha.each) {
        computeCI(v, y, sigma, truncation, alpha.each)
      }, FUN.VALUE = c(0.1, 0.2))

      int.candidate <- t(int.candidate)
      int.error <- stats::qnorm(p = 1 - (alpha-alpha.tilde)/2, mean = 0, sd = 1, lower.tail = TRUE) * sigma
      int.candidate[, 1] <- int.candidate[, 1] - int.error
      int.candidate[, 2] <- int.candidate[, 2] + int.error
      int.length <- int.candidate[, 2] - int.candidate[, 1]
      return(int.candidate[which.min(int.length), ])
    }, FUN.VALUE = c(0.1, 0.2))

    out[, 2:3] <- t(int)
    return(out)

  }
}



#' Form confidence intervals for regression coeffcients from an \code{"outference"} object
#'
#' This function constructs confidence intervals for the regression coefficients from an
#'     \code{"outference"} object. The syntax is the same as \code{\link{confint.lm}}.
#'
#' @export
#'
#' @param object, an object of class \code{"outference"}.
#' @param parm, indexes of which parameter to consider. If set to \code{NULL}, then
#'     all parameters are considered.
#' @param level, the confidence level.
#' @param ..., other arguments.
#'
#' @return A matrix with columns being lower and upper confidence limits for each parameter.
#'
#' @seealso \code{\link{outference}} for model fitting;
#'
#'     \code{\link{summary.outference}} for summaries;
#'
#'     \code{\link{coef.outference}} for extracting coefficients;
#'
#'     \code{\link{plot.outference}} for plotting the outlying measure;
#'
#'     \code{\link{predict.outference}} for making predictions.
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' confint(fit)
#'

confint.outference <- function(object, parm = NULL, level = 0.95, ...) {
  alpha <- 1-level
  # extract stuff to be used for inference
  y <- stats::model.response(object$fit.full$model)
  X <- stats::model.matrix(object$fit.full)
  n <- length(y)
  p <- ncol(X)
  outlier.det <- object$outlier.det
  sigma <- object$sigma
  if (is.null(sigma)) sigma <- estimateSigma(y, X)

  # which parameter to consider?
  label <- names(coef.outference(object, model = "removed"))
  if (is.null(parm)) parm <- 1:p

  res <- vapply(parm, function(j) {

    v <- contrastForCoef(j, X, outlier.det)
    z <- y - v * sum(v*y) / sum(v*v)
    truncation <- intForZAll(method = object$method, constraint = object$constraint,
                            v = v, z = z, sigma = sigma, outlier.det = outlier.det) # save the truncation
    return(computeCI(v, y, sigma, truncation, alpha))

  }, FUN.VALUE = c(0.1, 0.2))

  res <- t(res)
  rownames(res) <- label[parm]
  colnames(res) <- c(paste(alpha*100/2, '%'), paste((1-alpha/2)*100, '%'))
  return(res)

}








#' Test for a single regression coefficient from an \code{"outference"} object
#'
#' This function does a two-sided selective test for a single regression coefficient
#'     (test whether it is zero or not) from an \code{"outference"} object.
#'
#' @export
#'
#' @param object, an object of class \code{"outference"}.
#' @param index, the index of parameter under consideration. The intercept is labelled as \code{index = 1}.
#'
#' @return If \code{sigma} is known or estimated, then this function returns a list with the following
#'     components:
#'     \item{Z, }{the value of the \eqn{Z}-statistic.}
#'     \item{sd, }{the standard error of the regression coefficient.}
#'     \item{truncation, }{the truncation set for the \eqn{Z}-statistic.}
#'     \item{pval, }{the corrected p-value.}
#'     If \code{sigma = NULL}, then this function returns a list with the following components:
#'     \item{truncation, }{the truncation set for the \eqn{F}-statistic.}
#'     \item{pval, }{the corrected p-value.}
#'
#' @seealso \code{\link{grptest}} for testing group structures from an \code{"outference"} object.
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' coeftest(fit, 2)
#'
coeftest <- function(object, index) {
  sigma <- object$sigma
  # extract stuff to be used for inference
  y <- stats::model.response(object$fit.full$model)
  X <- stats::model.matrix(object$fit.full)
  n <- length(y)
  p <- ncol(X)
  if (index > p) {
    stop("index must be an element in 1:p!")
  }
  outlier.det <- object$outlier.det

  if (is.null(sigma)) { # do selective t test

    coef.fstatistic <- coefForFStatistic(y, X, outlier.det, g = index)
    truncation <- intForFAll(method = object$method, constraint = object$constraint,
                            wDelta = coef.fstatistic$wDelta, w2 = coef.fstatistic$w2,
                            z = coef.fstatistic$z, r = coef.fstatistic$r,
                            C = coef.fstatistic$C, outlier.det = outlier.det)
    pval <- TFSurv(q = coef.fstatistic$f, df1 = 1, df2 = length(coef.fstatistic$M) - p,
                  E = truncation)
    return(list(truncation = truncation, pval = pval))
  }
  else { # do selective z-test

    v <- contrastForCoef(index, X, outlier.det)
    vTy <- sum(v*y)
    v.norm <- sqrt(sum(v*v))
    truncation <- intForZAll(method = object$method, constraint = object$constraint,
                            v = v, z = y - v * vTy/v.norm^2, sigma = sigma, outlier.det = outlier.det) # save the truncation
    se <- sigma * v.norm
    Z <- vTy/se
    pval <- giveTwoSidePval(TNSurv(q = Z, mean = 0, sd = 1, E = truncation))

    return(list(Z = Z, se = se, truncation = truncation, pval = pval))
  }
}






#' Test group structures from an \code{"outference"} object
#'
#' This function does a selective test for group structures
#'     (test whether the regression coefficients indexed by \code{g} is all zero or not)
#'     from an \code{"outference"} object.
#'
#' @export
#'
#' @param object, an object of class \code{"outference"}.
#' @param group, a vector indicating which group to test. The intercept is labelled as \code{1}.
#'
#' @return This function returns a list with the following
#'     components:
#'     \item{statistic, }{the value of the chi-squared statistic or \eqn{F}-statistic,
#'         depending on whether \code{sigma} is known}
#'     \item{df, }{the degree(s) of freedom.}
#'     \item{truncation, }{the truncation set for the statistic.}
#'     \item{pval, }{the corrected p-value.}
#'
#'
#' @seealso \code{\link{coeftest}} for testing for a single regression coefficient
#'     from an \code{"outference"} object.
#'
#' @author Shuxiao Chen <sc2667@cornell.edu>
#'
#' @examples
#' ## Brownlee’s Stack Loss Plant Data
#' data("stackloss")
#' ## fit the model
#' ## detect outlier using Cook's distance with cutoff = 4
#' fit <- outference(stack.loss ~ ., data = stackloss, method = "cook", cutoff = 4)
#' grptest(fit, 2:4)
#'
grptest <- function(object, group) {

  sigma <- object$sigma
  # extract stuff to be used for inference
  y <- stats::model.response(object$fit.full$model)
  X <- stats::model.matrix(object$fit.full)
  n <- length(y)
  p <- ncol(X)
  if (sum(group %in% 1:p) != length(group)) {
    stop("group must be a subset of 1:p!")
  }
  outlier.det <- object$outlier.det

  if (is.null(sigma)) { # do selective F test
    coef.fstatistic <- coefForFStatistic(y, X, outlier.det, g = group)
    truncation <- intForFAll(method = object$method, constraint = object$constraint,
                            wDelta = coef.fstatistic$wDelta, w2 = coef.fstatistic$w2,
                            z = coef.fstatistic$z, r = coef.fstatistic$r,
                            C = coef.fstatistic$C, outlier.det = outlier.det)
    #f.global = coef.fstatistic$f
    df <- c(length(group), length(coef.fstatistic$M) - p)
    pval <- TFSurv(q = coef.fstatistic$f, df1 = df[1], df2 = df[2], E = truncation)
    return(list(statistic = coef.fstatistic$f, df = df, truncation = truncation, pval = pval))
  }
  else { # do selective chisq test
    coef.chisqstatistic <- coefForChisqStatistic(y, X, sigma, outlier.det, group)
    truncation <- intForChisqAll(method = object$method, constraint = object$constraint,
                                w = coef.chisqstatistic$w, z = coef.chisqstatistic$z, sigma = sigma,
                                outlier.det = outlier.det)
    pval <- TChisqSurv((coef.chisqstatistic$chi)^2, df = coef.chisqstatistic$df, E = truncation)

    return(list(statistic = (coef.chisqstatistic$chi)^2, df = coef.chisqstatistic$df, truncation = truncation, pval = pval))

  }

}


# a tiny test
#set.seed(2667)
#n = 100; p = 10; sigma = 1; cutoff = 4
#trueOutlier = 1:5; shift = rep(5, 5); shift = shift * c(1,1,1,-1,-1)
#beta = c(1, 2, -2, 3, -5, rep(0, 5))
# setup the design matrix
# X = matrix(rnorm(n*(p-1)), n, (p-1))
#X = scale(X, FALSE, TRUE)*sqrt(n/(n-1)) # X is centered and has norm sqrt(n)
#X = cbind(1, X)
#u = rep(0, n)
#u[trueOutlier] = shift
#mu = X %*% beta + u
#y = mu + rnorm(n, 0, sigma)
## the data for training
#data.train = data.frame(y = y, X[, -1])
# the data for prediction
#data.pred = data.frame(matrix(rnorm(10*(p-1)), 10, (p-1)))

# fit the model
#res.unknown = lm_outlier_removed(y~., data = data.train, method = "cook", sigma = NULL)
#res.known = lm_outlier_removed(y~., data = data.train, method = "cook", sigma = sigma)
#res.est = lm_outlier_removed(y~., data = data.train, method = "cook", sigma = "estimate")
#res.unknown
#res.known
#res.est
# extract the summary
#summary(res.unknown)
#summary(res.known)
#summary(res.est)
# intervals for each coefficients
#confint(res.unknown)
# intervals for surface and prediction
#predict(res.unknown, newdata = data.pred[1:5,], interval = "confidence", level = 0.95)
#predict(res.unknown, newdata = data.pred[1:5,], interval = "prediction", level = 0.95)

