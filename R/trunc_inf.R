# ----- truncation set in the response -----

#' Compute the truncation set in the response after
#'     outlier detection using Cook's distance
#'
#' This function computes the matrices \eqn{Q_i}, such that
#'     outlier detection using Cook's distance is equivalent to
#'     \eqn{\bigcap_{i \in [n]} y^T Q_i y \ge 0}.
#'
#' Using Cook's distance as a heuristic, the \eqn{i}-th data is considered
#'     as an outlier if and only if its Cook's distance is larger than \eqn{\lambda/n},
#'     where \eqn{lambda} is the user-specified cutoff. Then we can characterize this
#'     "detection event" as an intersection of quadratic constraints in the response \eqn{y}
#'     by \eqn{\bigcap_{i \in [n]} y^T Q_i y \ge 0}.
#'
#' @keywords internal
#'
#' @param n, the number of observations.
#' @param p, the number of variables, including the intercept.
#' @param PX, the projection matrix onto the column space of the design matrix \eqn{X}.
#' @param PXperp, \code{I - PX}.
#' @param outlier.det, indexes of detected outliers, can be empty.
#' @param cutoff, the cutoff \eqn{\lambda} (see details).
#'
#' @return This function returns a list of matrices \eqn{Q_i}.
#'

constrInResponseCook <- function(n, p, PX, PXperp, outlier.det, cutoff) {
  out <- lapply(1:n, function(i) {
    hi <- PX[i, i]
    PXperp_i <- PXperp[, i] # i-th column of PXperp
    Q <- (n-p) * hi * outer(PXperp_i, PXperp_i)
    Q <- Q - (cutoff*p*(1-hi)^2/n) * PXperp
    if (i %in% outlier.det) {
      return(Q)
    }
    else {
      return(-Q)
    }
  })
  return(out)
}


#' Compute the truncation set in the response after SEQUENTIAL
#'     outlier detection using Cook's distance
#'
#' This function computes the matrices \eqn{Q_i}, such that sequential
#'     outlier detection using Cook's distance is equivalent to
#'     \eqn{\bigcap_{i \in I} y^T Q_i y \ge 0}.
#'
#' Using Cook's distance sequentially and assume there are \eqn{k} outliers,
#'     the \eqn{i}-th data is considered
#'     as an outlier if and only if its Cook's distance ranks in top-\eqn{k} among all,
#'     Then we can characterize this
#'     "detection event" as an intersection of quadratic constraints in the response \eqn{y}
#'     by \eqn{\bigcap_{i \in I} y^T Q_i y \ge 0}, where \eqn{I} is a finite index set.
#'
#'
#' @keywords internal
#'
#' @param n, the number of observations.
#' @param p, the number of variables, including the intercept.
#' @param PX, the projection matrix onto the column space of the design matrix \eqn{X}.
#' @param PXperp, \code{I - PX}.
#' @param obs.ordered, the index of observations with their cook's distance in the decreasing order.
#' @param numOfOutlier, the number of outliers assumed, must be between \eqn{0} and \eqn{n}.
#'
#' @return This function returns a list of matrices \eqn{Q_i}.
#'

constrInResponseCookSeq <- function(n, p, PX, PXperp, obs.ordered, numOfOutlier) {
  if (numOfOutlier <= 0) stop("numOfOutlier must be positive!")
  if (numOfOutlier >= n) stop("numOfOutlier cannot be greater than n!")
  out <- list()
  count <- 0
  for (i in 1:numOfOutlier) {
    index.outlier <- obs.ordered[i]
    for (j in (i+1):n) {
      count <- count + 1
      index.nonOutlier <- obs.ordered[j] # index of obs that has smaller cook's distance than index.outlier
      h.outlier <- PX[index.outlier, index.outlier]
      h.nonOutlier <- PX[index.nonOutlier, index.nonOutlier]
      PXperp.outlier <- PXperp[, index.outlier] # index.outlier-th column of PXperp
      PXperp.nonOutlier <- PXperp[, index.nonOutlier]
      out[[count]] <- outer(PXperp.outlier, PXperp.outlier) * h.outlier / (1-h.outlier)^2 - outer(PXperp.nonOutlier, PXperp.nonOutlier) * h.nonOutlier / (1-h.nonOutlier)^2
    }
  }
  return(out)
}


#' Compute the truncation set in the response after
#'     outlier detection using DFFITS
#'
#' This function computes the matrices \eqn{Q_i}, such that
#'     outlier detection using DFFITS is equivalent to
#'     \eqn{\bigcap_{i \in [n]} y^T Q_i y \ge 0}.
#'
#' Using DFFITS as a heuristic, the \eqn{i}-th data is considered
#'     as an outlier if and only if the square of its DFFITS value
#'     is larger than \eqn{\lambda p/(n-p)},
#'     where \eqn{lambda} is the user-specified cutoff. Then we can characterize this
#'     "detection event" as an intersection of quadratic constraints in the response \eqn{y}
#'     by \eqn{\bigcap_{i \in [n]} y^T Q_i y \ge 0}.
#'
#' @keywords internal
#'
#' @param n, the number of observations.
#' @param p, the number of variables, including the intercept.
#' @param PX, the projection matrix onto the column space of the design matrix \eqn{X}.
#' @param PXperp, \code{I - PX}.
#' @param outlier.det, indexes of detected outliers, can be empty.
#' @param cutoff, the cutoff \eqn{\lambda} (see details).
#'
#' @return This function returns a list of matrices \eqn{Q_i}.
#'

constrInResponseDffits <- function(n, p, PX, PXperp, outlier.det, cutoff) {
  out <- lapply(1:n, function(i) {
    hi <- PX[i, i]
    PXperp_i <- PXperp[, i] # i-th column of PXperp
    Q <- cutoff*p/(n-p) * PXperp
    Q <- Q - (hi*(n-p-1)/(1-hi)^2 + cutoff*p/((n-p)*(1-hi))) * outer(PXperp_i, PXperp_i)
    if (i %in% outlier.det) {
      return(-Q)
    }
    else {
      return(Q)
    }
  })
  return(out)
}






#' Compute the truncation set in the response after
#'     outlier detection using lasso
#'
#' This function computes a matrix \eqn{A} and a vector \eqn{b}, such that
#'     outlier detection using lasso is equivalent to
#'     \eqn{A y \ge b}.
#'
#' Consider solving the following program
#'     \deqn{minimize ||y-X\beta-u||_2^2/(2n) + \lambda ||u||_1.}
#'     The \eqn{i}-th observation is considered as an outlier
#'     if and only if \eqn{\hat u_i \neq 0}.
#'     This is equivalent to solving
#'     \deqn{minimize ||P_X^\perp (y-u)||_2^2/(2n) + \lambda ||u||_1.}
#'     Then the variable selection can be characterized by a set of affine constraints
#'     \eqn{Ay \ge b}. In essence, this function is equivalent to
#'     \code{selectiveInference:::fixedLasso.poly}, but adapted to our notations
#'     and up to some scaling factors and linear transformations.
#' @keywords internal
#'
#' @param n, the number of observations.
#' @param p, the number of variables, including the intercept.
#' @param PXperp, the projection matrix onto the orthogonal complement
#'     of the column space of the design matrix \eqn{X}.
#' @param outlier.det, indexes of detected outliers, can be empty.
#' @param outlier.det.sign, the sign of the active variable estimated by lasso.
#' @param cutoff, the cutoff \eqn{\lambda} (see details).
#'
#' @return This function returns a list (A, b).
#'
#' @references Lee, Jason D., et al. "Exact post-selection inference, with application to the lasso."
#'     The Annals of Statistics 44.3 (2016): 907-927.
#' @references Tibshirani, R., et al. "selectiveInference: Tools for Post-Selection Inference."
#'     R package version 1.3 (2016).
#'
#'
constrInResponseLasso <- function(n, p, PXperp, outlier.det, outlier.det.sign, cutoff) {
  if (length(outlier.det) == 0) {
    A <- 1/(n*cutoff) * rbind(-PXperp, PXperp)
    b <- rep(-1, nrow(A))
    # we should have return A = A %*% PXperp, but they are the same!
    return(list(A = A, b = b))
  }
  else {
    # 0 < length(outlier.det) < n (if == n, then lm_outlier_removed should have already thrown an error)
    M <- (1:n)[-outlier.det]
    # PXperp.Mc <- PXperp[, outlier.det]
    # PXperp.Mc.crosspd <- crossprod(PXperp.Mc)
    # PXperp.Mc.crosspd.inv <- chol2inv(chol(PXperp.Mc.crosspd))
    # PXperp.Mc.ginv <- PXperp.Mc.crosspd.inv %*% t(PXperp.Mc)
    # proj.Mc <- PXperp.Mc %*% PXperp.Mc.ginv
    # proj.Mc.perp <- diag(n) - proj.Mc
    # temp <- 1/(n*cutoff) * t(PXperp[, M]) %*% proj.Mc.perp
    # A <- rbind(-temp, temp)
    # temp <- t(PXperp[, M]) %*% t(PXperp.Mc.ginv) %*% outlier.det.sign
    # b <- rbind(-1 + temp, -1 - temp)
    # if (length(outlier.det) == 1) {
    #   diag.sign <- as.matrix(outlier.det.sign)
    # }
    # else {
    #   diag.sign <- diag(outlier.det.sign)
    # }
    # A <- rbind(A, diag.sign %*% PXperp.Mc.ginv)
    # # translate to polyhedron for y!
    # A <- A %*% PXperp
    # b <- as.numeric(rbind(b, cutoff*n * diag.sign %*% PXperp.Mc.crosspd.inv %*% outlier.det.sign))

    PXperp.Mc.inv <- chol2inv(chol(PXperp[-M, -M]))
    MMc.times.PXperp.Mc.inv <- PXperp[M, -M] %*% PXperp.Mc.inv
    temp <- -PXperp[M, ] + MMc.times.PXperp.Mc.inv %*% PXperp[-M, ]
    A <- 1/(n*cutoff) * rbind(temp, -temp) # A0
    temp <- as.numeric(MMc.times.PXperp.Mc.inv %*% outlier.det.sign)
    b <- c(-1+temp, -1-temp) # b0
    if (length(outlier.det) == 1) {
      diag.sign <- as.matrix(outlier.det.sign)
    }
    else {
      diag.sign <- diag(outlier.det.sign)
    }
    A <- rbind(A, diag.sign  %*% PXperp.Mc.inv %*% PXperp[-M, ]) # stack A0 and A1 together
    A <- A %*% PXperp # polyhedron is w.r.t. y (not PXperp %*% y!)
    temp <- as.numeric(n*cutoff * diag.sign %*% PXperp.Mc.inv %*% outlier.det.sign) # b1
    b <- c(b, temp)

  }

  return(list(A = A, b = b))

}


# ----- functions related to testing nu^T y = 0 when sigma is known -----


#' Compute the contrast vector for inference of regression coefficients
#'
#' This function computes the contrast vector \eqn{\nu} such that
#'     \deqn{\nu^T y = \beta^M_j.}
#'
#' @keywords internal
#'
#' @param j, the index of regression coefficients.
#' @param X, the design matrix (including the intercept).
#' @param outlier.det, indexes of detected outliers.
#'
#' @return This function returns the contrast vector \eqn{\nu}.
contrastForCoef <- function(j, X, outlier.det) {
  n <- nrow(X)
  p <- ncol(X)

  # if length(outlier.det) == 0, this function returns as.numeric(t(t(ej) %*% MASS::ginv(X)))
  # else, this function returns
  # as.numeric(t(t(ej) %*% MASS::ginv(X[-outlier.det, ]) %*% diag(n)[-outlier.det, ]))

  if (length(outlier.det) == 0) {
    temp <- t(MASS::ginv(X))
    return(temp[, j])
  }
  # there is at least one outlier detected
  temp <- matrix(0, n, p)
  temp[-outlier.det, ] <- t(MASS::ginv(X[-outlier.det, ]))
  return(temp[, j])
}


# return v such that v^T y = x0^T beta^N

#' Compute the contrast vector for inference of regression surfaces
#'
#' This function computes the contrast vector \eqn{\nu} such that
#'     \deqn{\nu^T y = x_0^T\beta^M,}
#'     where \eqn{x_0} is a \eqn{p}-dimensional new data point.
#'
#' @keywords internal
#'
#' @param x0, the new data point of length \eqn{p}.
#' @param X, the design matrix (including the intercept).
#' @param outlier.det, indexes of detected outliers.
#'
#' @return This function returns the contrast vector \eqn{\nu}.
#'
contrastForSurf = function(x0, X, outlier.det) {
  n = nrow(X)
  p = ncol(X)

  # if length(outlier.det) == 0, this function returns as.numeric(t(t(x0) %*% MASS::ginv(X))),
  # else, this function returns
  # as.numeric(t(t(x0) %*% MASS::ginv(X[-outlier.det, ]) %*% diag(n)[-outlier.det, ]))

  if (length(outlier.det) == 0) {
    temp <- t(MASS::ginv(X))
    return(temp %*% x0)
  }
  # there is at least one outlier detected
  temp <- matrix(0, n, p)
  temp[-outlier.det, ] <- t(MASS::ginv(X[-outlier.det, ]))
  return(as.numeric(temp %*% x0))
}



# compute coefficients which define *one slice* of the truncation set for Z = v^T y / (sigma * ||v||_2)
#' Compute coeffcients that defines one slice of the truncation set for the \eqn{Z} statistic
#'
#' This function computes the coefficients \eqn{(A, B, C)} that defines
#'     one slice of the truncation set for the \eqn{Z} statistic.
#'
#' Consider the \eqn{Z}-statistic \deqn{Z = \nu^T y / (\sigma ||\nu||_2).} The constraint
#'     in the response \eqn{{y^T Q y + a^T y + b \ge 0}} is equivalent to
#'     \eqn{{AZ^2 + BZ + C \ge 0}} for some scaler coefficients \eqn{(A, B, C)}.
#'
#' @keywords internal
#'
#' @param Q, a matrix.
#' @param a, a vector.
#' @param b, a scaler.
#' @param v, the contrast vector \eqn{\nu}.
#' @param z, \eqn{P_\nu^\perp y}.
#' @param sigma, the noise level \eqn{\sigma}.
#'
#' @return This function returns a vector of scaler coefficients \eqn{(A, B, C)}.
coefForZEachSlice <- function(Q, a = NULL, b = NULL, v, z, sigma) {
  Qv <- Q %*% v
  Qz <- Q %*% z
  vTv <- sum(v^2)
  A <- sigma^2 * sum(v * Qv) / vTv
  # in outlier detection setup, a, b usually equal to 0
  B <- 2 * sigma * sum(v * Qz) / sqrt(vTv)
  C <- sum(z * Qz)

  if (is.null(a) & is.null(b)) {
    return(c(A, B, C))
  }

  # the general form
  B <- B + sigma * sum(a * v) / sqrt(vTv)
  C <- C + sum(a*z) + b
  return(c(A, B, C))
}



# # this function is kept for reference,
# # see the documentation for compIntForZEachSlice
# intForZEachSlice <- function(coef, tol = 1e-15) {
#   A <- coef[1]
#   B <- coef[2]
#   C <- coef[3]
#   # return interval for x,
#   # s.t. Bx + C >= 0
#   affineIntForZEachSlice <- function(B, C, tol) {
#     # degenerate case: B = 0
#     if (abs(B) <= tol) {
#       #warning("A=B=0, ill-conditioned polynomial")
#       if (C >= -tol) { # C >=0 -> hold for all X
#         return(Intervals(c(-Inf, Inf)))
#       }
#       else {
#         return(Intervals()) # empty interval
#       }
#     }
#
#     # now we know B != 0
#     temp <- -C/B
#     if (B > tol) { # B > 0 -> X >= -C/B
#       return(Intervals(c(temp, Inf)))
#     }
#
#     # now we know B < 0 -> X <= -C/B
#     return(Intervals(c(-Inf, temp)))
#   }
#
#
#   # denegerate case: A = 0
#   if (abs(A) <= tol) {
#     return(affineIntForZEachSlice(B = B, C = C, tol = tol))
#   }
#
#   # now we know A != 0
#   disc <- B^2 - 4*A*C
#   if (disc <= tol) { # discriminant <= 0
#
#     # notice: we ignore the case when truncation is a single point
#
#     if (A > tol) { # parabola open upwards
#       # polynomial always >= 0
#       return(Intervals(c(-Inf, Inf)))
#     }
#     else { # parabola open downwars
#       # polynomial always <= 0 -> no such X exists
#       return(Intervals())
#     }
#   }
#
#   # now we know A !=0 & disc > 0
#   negB2A <- -B/(2*A)
#   rtDisc2A <- sqrt(disc)/(2*A)
#   # two real roots s.t. root1 < root2
#   roots <- sort(c(negB2A - rtDisc2A, negB2A + rtDisc2A))
#
#   if (A > tol) { # parabola open upwards
#     return(Intervals(rbind(c(-Inf, roots[1]), c(roots[2], Inf))))
#   }
#
#   # now we know A < 0 & disc > 0
#   # parabola open downwards
#   return(Intervals(roots))
# }


#' Solve the roots of quadratic polynomials related to truncated \eqn{Z} test
#'
#' This function solves the inequality \eqn{Ax^2 + Bx + C \ge 0}, then
#'     returns the complement of the solution set (w.r.t. to \eqn{\mathbb{R}}).
#'
#' The reason we compute the complement instead of the original solution set
#'     is that when taking intersection of multiple sets, the strategy of
#'     "taking union of complements, then taking comlement" is substantially faster
#'     than taking intersections directly based on the \code{intervals} package.
#'
#' @keywords internal
#'
#' @param coef, the coefficients \eqn{(A, B, C)} of the quadratic polynomial.
#' @param tol, the tolerance of roots.
#'
#' @return This function returns an "Intervals" object.
compIntForZEachSlice <- function(coef, tol = 1e-15) {
  A <- coef[1]
  B <- coef[2]
  C <- coef[3]
  # return interval for x,
  # s.t. Bx + C >= 0
  affineIntForZEachSlice <- function(B, C, tol) {
    # degenerate case: B = 0
    if (abs(B) <= tol) {
      #warning("A=B=0, ill-conditioned polynomial")
      if (C >= -tol) { # C >=0 -> hold for all X
        return(intervals::Intervals())
      }
      else {
        return(intervals::Intervals(c(-Inf, Inf)))
      }
    }

    # now we know B != 0
    temp <- -C/B
    if (B > tol) { # B > 0 -> X >= -C/B
      return(intervals::Intervals(c(-Inf, temp)))
    }

    # now we know B < 0 -> X <= -C/B
    return(intervals::Intervals(c(temp, Inf)))
  }


  # denegerate case: A = 0
  if (abs(A) <= tol) {
    return(affineIntForZEachSlice(B = B, C = C, tol = tol))
  }

  # now we know A != 0
  disc <- B^2 - 4*A*C
  if (disc <= tol) { # discriminant <= 0

    # notice: we ignore the case when truncation is a single point

    if (A > tol) { # parabola open upwards
      # polynomial always >= 0
      return(intervals::Intervals())
    }
    else { # parabola open downwars
      # polynomial always <= 0 -> no such X exists
      return(intervals::Intervals(c(-Inf, Inf)))
    }
  }

  # now we know A !=0 & disc > 0
  negB2A <- -B/(2*A)
  rtDisc2A <- sqrt(disc)/(2*A)
  # two real roots s.t. root1 < root2
  roots <- sort(c(negB2A - rtDisc2A, negB2A + rtDisc2A))

  if (A > tol) { # parabola open upwards
    return(intervals::Intervals(roots))
  }

  # now we know A < 0 & disc > 0
  # parabola open downwards
  return(intervals::Intervals(rbind(c(-Inf, roots[1]), c(roots[2], Inf))))
}





# return the truncation interval for Z statistic
# recall z := P_v^\perp y
#' Compute the truncation set for \eqn{Z} statistic
#'
#' This function computes the truncation set for the \eqn{Z} statistc.
#'
#' Consider the \eqn{Z}-statistic \deqn{Z = \nu^T y / (\sigma ||\nu||_2).}
#'     This function translates the constraints in the response into the truncation
#'     set (which is a union of intervals) for the \eqn{Z} statistic.
#'
#' @keywords internal
#'
#' @param method, the outlier detection method, must be one of "cook", "dffits", "lasso".
#' @param constraint, the constraint in the response.
#' @param v, the contrast vector \eqn{\nu}.
#' @param z, \eqn{P_\nu^\perp y}.
#' @param sigma, the noise level \eqn{\sigma}.
#' @param outlier.det, indexes of detected outliers.
#'
#' @return This function returns an "Intervals" object.
#'
intForZAll <- function(method = c("cook", "dffits", "lasso"), constraint,
                       v, z, sigma, outlier.det) {
  n <- length(z)
  method <- match.arg(method)
  Pv <- outer(v, v)/(sum(v*v))
  Pvperp <- diag(n) - Pv

  if (method %in% c("cook", "dffits")) {

    coef.all <- lapply(constraint, coefForZEachSlice, v = v, z = z, sigma = sigma)
    truncation <- lapply(coef.all, compIntForZEachSlice)
    truncation <- intervals::interval_complement(do.call(intervals::interval_union, truncation)) # notice we are taking union then taking complement!
    return(truncation)
  }
  else { # method == "lasso"
    A <- constraint$A
    b <- constraint$b
    coef.all <- lapply(1:length(b), function(i) {
      coefForZEachSlice(Q = matrix(0, n, n), a = A[i, ], b = -b[i], v = v, z = z, sigma = sigma)
    })
    truncation <- lapply(coef.all, compIntForZEachSlice)
    truncation <- intervals::interval_complement(do.call(intervals::interval_union, truncation)) # notice we are taking union then taking complement!
    return(truncation)
  }

}


# ----- functions related to testing group structures when sigma is known -----


# return list(M, Mc, g, gc, P, chi, w, z, df)
# those are used in computing the truncation for f statistic.
# require: length(outlier.det) < length(y)

#' Compute the coefficients related to truncated chi-squared test
#'
#' This function computes the coefficients related to the truncated chi-squared test.
#'
#' @keywords internal
#'
#' @param y, the response.
#' @param X, the design matrix (including intercepts).
#' @param sigma, the noise level \eqn{\sigma}.
#' @param outlier.det, indexes of outliers detected.
#' @param g, a subset of \eqn{[p]} indicating which group structure to test.
#'
#' @return This function returns a list of \code{(M, Mc, g, gc, P, chi, w, z, df)}.
#'     See the manuscript for details of those quantities.
coefForChisqStatistic <- function(y, X, sigma, outlier.det, g) {
  n <- length(y)
  p <- ncol(X)
  gc <- (1:p)[-g]
  if (length(outlier.det) == 0) {
    M <- 1:n
    Mc <- integer(0)
  }
  else {
    M <- (1:n)[-outlier.det]
    Mc <- outlier.det
  }
  gc <- (1:p)[-g]
  Xtilde <- (diag(length(M)) - proj(X[M, gc])) %*% X[M, g]
  P <- matrix(0, n, n)
  P[M, M] <- proj(Xtilde)
  chi <- sqrt(sum(y * P%*%y))/sigma
  w <- P%*%y/(sigma*chi)
  z <- (diag(n) - P)%*%y
  df <- sum(diag(P))
  return(list(M = M, Mc = Mc, g = g, gc = gc, P = P, chi = chi, w = w, z = z, df = df))
}


#' Compute coeffcients that defines one slice of the truncation set for the chi-squared statistic
#'
#' This function computes the coefficients \eqn{(A, B, C)} that defines
#'     one slice of the truncation set for the chi-squared statistic.
#'
#' The constraint in the response \eqn{{y^T Q y + a^T y + b \ge 0}} is equivalent to
#'     \eqn{{AZ^2 + BZ + C \ge 0}} for some scaler coefficients \eqn{(A, B, C)}.
#'
#' @keywords internal
#'
#' @param Q, a matrix.
#' @param a, a vector.
#' @param b, a scaler.
#' @param w,z, coefficients related to chi-squared statistic.
#' @param sigma, the noise level \eqn{\sigma}.
#'
#' @return This function returns a vector of scaler coefficients \eqn{(A, B, C)}.
coefForChisqEachSlice <- function(Q, a = NULL, b = NULL, w, z, sigma) {
  Qw <- Q %*% w
  Qz <- Q %*% z
  A <- sigma^2 * sum(w * Qw) # sigma^2 * w^T Q w
  # in outlier detection setup, a, b usually equal to 0
  B <- 2 * sigma * sum(z * Qw)
  C <- sum(z * Qz)
  if (is.null(a) & is.null(b)) {
    return(c(A, B, C))
  }

  # the general form
  B <- B + sigma * sum(a * w)
  C <- C + sum(a * z) + b

  return(c(A, B, C))
}

# # This function is kept for reference.
# # see the documentation for compIntForChisqEachSlice
# # return inverval for X^2,
# # s.t. AX^2 + BX + C >= 0 AND X >= 0
# intForChisqEachSlice <- function(coef, tol = 1e-15) {
#   A <- coef[1]
#   B <- coef[2]
#   C <- coef[3]
#   # return interval for X^2,
#   # s.t. BX + C >= 0
#   affineIntForChisqEachSlice <- function(B, C, tol) {
#     # degenerate case: B = 0
#     if (abs(B) <= tol) {
#       #warning("A=B=0, ill-conditioned polynomial")
#       if (C >= -tol) { # C >=0 -> hold for all X
#         return(Intervals(c(0, Inf)))
#       }
#       else {
#         return(Intervals())
#       }
#     }
#
#     # now we know B != 0
#     temp <- -C/B
#     if (B > tol) { # B > 0 -> X >= -C/B
#
#       # if -C/B <=0, then X^2 >= 0
#       # if -C/B > 0, then X^2 >= (-C/B)^2
#       endpt1 <- (max(0, temp))^2
#       return(Intervals(c(endpt1, Inf)))
#     }
#
#     # now we know B < 0 -> X <= -C/B
#     # if -C/B < 0, then constraint infeasible b/c we already have X >=0
#     if (temp < -tol) {return(Intervals())}
#     # we know -C/B >= 0, which gives 0<= X^2 <= (-C/B)^2
#     return(Intervals(c(0, temp^2)))
#   }
#
#   # denegerate case: A = 0
#   if (abs(A) <= tol) {
#     return(affineIntForChisqEachSlice(B = B, C = C, tol = tol))
#   }
#
#   # now we know A != 0
#   disc <- B^2 - 4*A*C
#   if (disc <= tol) { # discriminant <= 0
#
#     # notice: we ignore the case when truncation is a single point
#
#     if (A > tol) { # parabola open upwards
#       # polynomial always >= 0 -> X^2 >= 0
#       return(Intervals(c(0, Inf)))
#     }
#     else { # parabola open downwards
#       # polynomial always <= 0 -> no such X exists
#       return(Intervals())
#     }
#   }
#
#   # now we know A !=0 & disc > 0
#   negB2A <- -B/(2*A)
#   rtDisc2A <- sqrt(disc)/(2*A)
#   # two real roots s.t. root1 < root2
#   roots <- sort(c(negB2A - rtDisc2A, negB2A + rtDisc2A))
#
#   if (A > tol) { # parabola open upwards
#     if (roots[1] > tol) { # 0 < root1 < root2
#       return(Intervals(rbind(c(0, roots[1]^2), c(roots[2]^2, Inf))))
#     }
#     # we know root1 <= 0
#     if (roots[2] > tol) { # root1 <= 0 < root2
#       return(Intervals(c(roots[2]^2, Inf)))
#     }
#
#     # we know root1 <= 0 AND root2 <= 0
#     return(Intervals(c(0, Inf)))
#   }
#
#   # now we know A < 0 & disc > 0
#   # parabola open downwards
#   # root1 <= X <= root2
#   if (roots[1] > tol) { # 0 < root1 < root2
#     return(Intervals(roots^2))
#   }
#   # we know root1 <= 0
#   if (roots[2] > tol) { # root1 <= 0 < root2
#     return(Intervals(c(0, roots[2]^2)))
#   }
#   # we know root1 <= 0 AND root2 <= 0
#   return(Intervals())
#
# }



#' Solve the roots of quadratic polynomials related to truncated chi-squared test
#'
#' This function solves \eqn{{x^2: Ax^2 + Bx + C \ge 0}}, then
#'     returns the complement of the solution set (w.r.t. to \eqn{\mathbb{R}}).
#'
#' The reason we compute the complement instead of the original solution set
#'     is that when taking intersection of multiple sets, the strategy of
#'     "taking union of complements, then taking comlement" is substantially faster
#'     than taking intersections directly based on the \code{intervals} package.
#'
#' @keywords internal
#'
#' @param coef, the coefficients \eqn{(A, B, C)} of the quadratic polynomial.
#' @param tol, the tolerance of roots.
#'
#' @return This function returns an "Intervals" object.
compIntForChisqEachSlice <- function(coef, tol = 1e-15) {
  A <- coef[1]
  B <- coef[2]
  C <- coef[3]
  # return interval for X^2,
  # s.t. BX + C >= 0
  affineIntForChisqEachSlice <- function(B, C, tol) {
    # degenerate case: B = 0
    if (abs(B) <= tol) {
      #warning("A=B=0, ill-conditioned polynomial")
      if (C >= -tol) { # C >=0 -> hold for all X
        return(intervals::Intervals(c(-Inf, 0)))
      }
      else {
        return(intervals::Intervals(c(-Inf, Inf)))
      }
    }

    # now we know B != 0
    temp <- -C/B
    if (B > tol) { # B > 0 -> X >= -C/B

      # if -C/B <=0, then X^2 >= 0
      # if -C/B > 0, then X^2 >= (-C/B)^2
      endpt1 <- (max(0, temp))^2
      return(intervals::Intervals(c(-Inf, endpt1)))
    }

    # now we know B < 0 -> X <= -C/B
    # if -C/B < 0, then constraint infeasible b/c we already have X >=0
    if (temp < -tol) { return(intervals::Intervals(c(-Inf, Inf)))}
    # we know -C/B >= 0, which gives 0<= X^2 <= (-C/B)^2
    return(intervals::Intervals(rbind(c(temp^2, Inf), c(-Inf, 0))))
  }

  # denegerate case: A = 0
  if (abs(A) <= tol) {
    return(affineIntForChisqEachSlice(B = B, C = C, tol = tol))
  }

  # now we know A != 0
  disc <- B^2 - 4*A*C
  if (disc <= tol) { # discriminant <= 0

    # notice: we ignore the case when truncation is a single point

    if (A > tol) { # parabola open upwards
      # polynomial always >= 0 -> X^2 >= 0
      return(intervals::Intervals(c(-Inf, 0)))
    }
    else { # parabola open downwards
      # polynomial always <= 0 -> no such X exists
      return(intervals::Intervals(c(-Inf, Inf)))
    }
  }

  # now we know A !=0 & disc > 0
  negB2A <- -B/(2*A)
  rtDisc2A <- sqrt(disc)/(2*A)
  # two real roots s.t. root1 < root2
  roots <- sort(c(negB2A - rtDisc2A, negB2A + rtDisc2A))

  if (A > tol) { # parabola open upwards
    if (roots[1] > tol) { # 0 < root1 < root2
      return(intervals::Intervals(rbind(roots^2, c(-Inf, 0))))
    }
    # we know root1 <= 0
    if (roots[2] > tol) { # root1 <= 0 < root2
      return(intervals::Intervals(c(-Inf, roots[2]^2)))
    }

    # we know root1 <= 0 AND root2 <= 0
    return(intervals::Intervals(c(-Inf, 0)))
  }

  # now we know A < 0 & disc > 0
  # parabola open downwards
  # root1 <= X <= root2
  if (roots[1] > tol) { # 0 < root1 < root2
    return(return(intervals::Intervals(rbind(c(-Inf, roots[1]^2), c(roots[2]^2, Inf)))))
  }
  # we know root1 <= 0
  if (roots[2] > tol) { # root1 <= 0 < root2
    return(intervals::Intervals(rbind(c(roots[2]^2, Inf), c(-Inf, 0))))
  }
  # we know root1 <= 0 AND root2 <= 0
  return(intervals::Intervals(c(-Inf, Inf)))

}




#' Compute the truncation set for chi-squared statistic
#'
#' This function computes the truncation set for the chi-squared statistc.
#'
#' This function translates the constraints in the response into the truncation
#'     set (which is a union of intervals) for the chi-squared statistic.
#'
#' @keywords internal
#'
#' @param method, the outlier detection method, must be one of "cook", "dffits", "lasso".
#' @param constraint, the constraint in the response.
#' @param w,z, coefficiens related to chi-squared statistic
#' @param sigma, the noise level \eqn{\sigma}.
#' @param outlier.det, indexes of detected outliers.
#'
#' @return This function returns an "Intervals" object.
#'

intForChisqAll <- function(method = c("cook", "dffits", "lasso"), constraint, w, z, sigma, outlier.det) {
  n <- length(z)
  method <- match.arg(method)

  if (method %in% c("cook", "dffits")) {

    coef.all <- lapply(constraint, coefForChisqEachSlice, w = w, z = z, sigma = sigma)
    truncation <- lapply(coef.all, compIntForChisqEachSlice)
    truncation <- do.call(intervals::interval_union, truncation)
    truncation <- intervals::interval_union(truncation, intervals::Intervals(c(-Inf, 0))) # an ad-hoc check
    truncation <- intervals::interval_complement(truncation)
    return(truncation)
  }
  else {
    A <- constraint$A
    b <- constraint$b

    coef.all <- lapply(1:length(b), function(i) {
      coefForChisqEachSlice(Q = matrix(0, n, n), a = A[i, ], b = -b[i], w = w, z = z, sigma = sigma)
    })
    truncation <- lapply(coef.all, compIntForChisqEachSlice)
    truncation <- do.call(intervals::interval_union, truncation)
    truncation <- intervals::interval_union(truncation, intervals::Intervals(c(-Inf, 0))) # an ad-hoc check
    truncation <- intervals::interval_complement(truncation)
    return(truncation)
  }
}


# ----- functions related to testing group structures when sigma is unknown -----

#' Compute the coefficients related to truncated \eqn{F} test
#'
#' This function computes the coefficients related to the truncated \eqn{F} test.
#'
#' @keywords internal
#'
#' @param y, the response.
#' @param X, the design matrix (including intercepts).
#' @param outlier.det, indexes of outliers detected.
#' @param g, a subset of \eqn{[p]} indicating which group structure to test.
#'
#' @return This function returns a list of \code{(M, Mc, g, gc, Psub, Pfull, f, wDelta, w2, z, r, C)}.
#'     See the manuscript for details of those quantities.
coefForFStatistic <- function(y, X, outlier.det, g) {
  n <- length(y)
  p <- ncol(X)
  gc <- (1:p)[-g]
  if (length(outlier.det) == 0) {
    M <- 1:n
    Mc <- integer(0)
  }
  else {
    M <- (1:n)[-outlier.det]
    Mc <- outlier.det
  }

  Psub <- diag(n)
  Psub[M, M] <- proj(X[M, gc])
  Pfull <- diag(n)
  Pfull[M, M] <- proj(X[M, ])
  R1 <- (diag(n) - Psub) %*% y
  R2 <- (diag(n) - Pfull) %*% y
  C <- length(g) / (length(M) - p)
  R1.norm <- sqrt(sum(R1 * R1))
  R2.norm <- sqrt(sum(R2 * R2))
  f <- (R1.norm^2 - R2.norm^2) / (C * R2.norm^2)
  wDelta <- (R1 - R2) / sqrt(sum((R1 - R2)^2))
  w2 <- R2 / R2.norm
  z <- Psub %*% y
  r <- R1.norm
  return(list(M = M, Mc = Mc, g = g, gc = gc, Psub = Psub, Pfull = Pfull, f = f,
              wDelta = wDelta, w2 = w2, z = z, r = r, C = C))
}




#' Compute coeffcients that defines one slice of the truncation set for the \eqn{F} statistic
#'
#' This function computes the coefficients \eqn{(x_11, x_12, x_22, x_1, x_2, x_0)} that defines
#'     one slice of the truncation set for the \eqn{F} statistic.
#'
#' The constraint in the response \eqn{{y^T Q y + a^T y + b \ge 0}} is equivalent to
#'     an inequality constraint in the \eqn{F} statistic, where
#'     the inequality is parameterized by \eqn{(x_11, x_12, x_22, x_1, x_2, x_0)}.
#'
#' @keywords internal
#'
#' @param Q, a matrix.
#' @param a, a vector.
#' @param b, a scaler.
#' @param wDelta,w2,z,r coefficients related to \eqn{F} statistic.
#' @param tol, if the absolute value of a quantity is less than \code{tol},
#'     then set it to \eqn{0}.
#'
#' @return This function returns a list of scaler coefficients
#'     \code{(x11, x12, x22, x1, x2, x0)}.
#'
coefForFEachSlice <- function(Q, a, b, wDelta, w2, z, r, tol = 1e-10) {
  QwDelta <- Q %*% wDelta
  Qw2 <- Q %*% w2
  Qz <- Q %*% z
  r2 <- r^2

  x11 <- r2 * sum(wDelta * QwDelta)
  if (abs(x11) < tol) x11 = 0
  x12 <- 2 * r2 * sum(wDelta * Qw2)
  if (abs(x12) < tol) x12 = 0
  x22 <- r2 * sum(w2 * Qw2)
  if (abs(x22) < tol) x22 = 0
  x1 <- r * (2 * sum(wDelta * Qz) + sum(a * wDelta))
  if (abs(x1) < tol) x1 = 0
  x2 <- r * (2 * sum(w2 * Qz) + sum(a * w2))
  if (abs(x2) < tol) x2 = 0
  x0 <- sum(z * Qz) + sum(a * z) + b
  if (abs(x0) < tol) x0 = 0
  out <- list(x11 = x11, x12 = x12, x22 = x22, x1 = x1, x2 = x2, x0 = x0)
  return(out)
}





#' Solve the complement of one slice of the truncation set for the \eqn{F} statistic
#'
#' This function computes one slice of the truncation set for the \eqn{F} statistic, then
#'     returns the complement of the truncation set (w.r.t. to \eqn{\mathbb{R}}).
#'
#' This function is essentially a wrapper function for
#'     \code{selectiveInference:::TF_roots}, with minor changes so that
#'     the function does not stop when the constraint is infeasible (instead, this
#'     function returns an empty "Intervals" object).
#'     The reason we compute the complement instead of the original solution set
#'     is that when taking intersection of multiple sets, the strategy of
#'     "taking union of complements, then taking comlement" is substantially faster
#'     than taking intersections directly based on the \code{intervals} package.
#'
#' @keywords internal
#'
#' @param C,coeffs, the coefficients related to the \eqn{F} statistic.
#' @param tol,tol2 the tolerance of roots.
#'
#' @return This function returns an "Intervals" object.
#'
#' @references Tibshirani, R., et al. "selectiveInference: Tools for Post-Selection Inference."
#'     R package version 1.3 (2016).
#'
compIntForFEachSlice <- function(C, coeffs, tol = 1e-8, tol2 = 1e-6) {

  # Helper functions for TF roots
  roots_to_checkpoints <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(c(0, (checkpoints + c(checkpoints[-1], 200 + checkpoints[length(checkpoints)]))/2))
  }
  roots_to_partition <- function(roots) {
    checkpoints <- unique(sort(c(0, roots)))
    return(list(endpoints = c(checkpoints, Inf), midpoints = (checkpoints + c(checkpoints[-1], 200 + checkpoints[length(checkpoints)]))/2))
  }

  x11 <- coeffs$x11
  x22 <- coeffs$x22
  x12 <- coeffs$x12
  x1 <- coeffs$x1
  x2 <- coeffs$x2
  x0 <- coeffs$x0

  g1 <- function(t) sqrt(C*t/(1+C*t))
  g2 <- function(t) 1/sqrt(1+C*t)
  I <- function(t) x11*g1(t)^2 + x12*g1(t)*g2(t) + x22*g2(t)^2 + x1*g1(t) + x2*g2(t) + x0

  z4 <- complex(real = -x11 + x22, imaginary = -x12)/4
  z3 <- complex(real = x2, imaginary = -x1)/2
  z2 <- complex(real = x11/2+x22/2+x0)
  z1 <- Conj(z3)
  z0 <- Conj(z4)

  zcoefs <- c(z0, z1, z2, z3, z4)
  croots <- polyroot(zcoefs)
  thetas <- Arg(croots)
  # Can't specify polyroot precision :(
  modinds <- Mod(croots) <= 1 + tol2 & Mod(croots) >= 1 - tol2
  angleinds <- thetas >=0 & thetas <= pi/2
  roots <- unique(thetas[which(modinds & angleinds)])
  troots <- tan(roots)^2/C

  checkpoints <- c()
  if (length(troots) > 0) checkpoints <- roots_to_checkpoints(troots)
  checkpoints <- sort(
    c(checkpoints, 0, tol, tol2,
      seq(from = sqrt(tol2), to = 1, length.out = 50),
      seq(from = 1.2, to=50, length.out = 20),
      100, 1000, 10000))
  ## if (length(troots) == 0) {
  ##     # Polyroot didn't catch any roots
  ##     # ad-hoc check:
  ##     checkpoints <- c(0, tol, tol2,
  ##                      seq(from = sqrt(tol2), to = 1, length.out = 50),
  ##                      seq(from = 1.2, to=50, length.out = 20),
  ##                      100, 1000, 10000)
  ## } else {
  ##     checkpoints <- roots_to_checkpoints(troots)
  ## }

  signs <- sign(I(checkpoints))
  diffs <- c(0, diff(signs))
  changeinds <- which(diffs != 0)

  if (length(changeinds) > 0) {

    roots <- unlist(lapply(changeinds, function(ind) {
      stats::uniroot(I, lower = checkpoints[ind-1], upper = checkpoints[ind], tol = tol)$root
    }))

    partition <- roots_to_partition(roots)
    negative <- which(I(partition$midpoints) < 0)

    intervals <- matrix(NA, ncol=2)
    for (i in 1:length(negative)) {
      ind <- negative[i]
      if ((i > 1) && (ind == negative[i-1] + 1)) {
        # There was not a sign change at end of previous interval
        intervals[nrow(intervals), 2] <- partition$endpoints[ind+1]
      } else {
        intervals <- rbind(intervals, c(partition$endpoints[ind], partition$endpoints[ind+1]))
      }
    }

    return(intervals::Intervals(rbind(intervals[-1,], c(-Inf, 0))))
  }


  if (I(0) < 0) return(intervals::Intervals(c(-Inf, Inf))) # always negative
  return(intervals::Intervals(c(-Inf,0))) # Apparently no roots, always positive
}





#' Compute the truncation set for the \eqn{F} statistic
#'
#' This function computes the truncation set for the \eqn{F} statistc.
#'
#' This function translates the constraints in the response into the truncation
#'     set (which is a union of intervals) for the \eqn{F} statistic.
#'
#' @keywords internal
#'
#' @param method, the outlier detection method, must be one of "cook", "dffits", "lasso".
#' @param constraint, the constraint in the response.
#' @param wDelta,w2,z,r,C, coefficiens related to the \eqn{F} statistic.
#' @param outlier.det, indexes of detected outliers.
#'
#' @return This function returns an "Intervals" object.
#'
intForFAll <- function(method = c("cook", "dffits", "lasso"), constraint,
                       wDelta, w2, z, r, C, outlier.det) {
  n <- length(z)
  method <- match.arg(method)

  if (method %in% c("cook", "dffits")) {
    coef.all <- lapply(constraint, coefForFEachSlice, a = rep(0, n), b = 0, wDelta = wDelta,
                      w2 = w2, z = z, r = r)
    truncation <- lapply(coef.all, function(coef.each) {
      compIntForFEachSlice(C = C, coeffs = coef.each)
    })
    truncation <- do.call(intervals::interval_union, truncation)
    truncation <- intervals::interval_union(truncation, intervals::Intervals(c(-Inf, 0))) # consider the implementation of TF_roots, we need to ensure truncation is inside (0, Inf)
    truncation <- intervals::interval_complement(truncation)
    return(truncation)
  }
  else { # method == "lasso"
    A <- constraint$A
    b <- constraint$b
    coef.all <- lapply(1:length(b), function(i) {
      coefForFEachSlice(Q = matrix(0, n, n), a = A[i, ], b = -b[i], wDelta = wDelta,
                        w2 = w2, z = z, r = r)
    })
    truncation <- lapply(coef.all, function(coef.each) {
      compIntForFEachSlice(C = C, coeffs = coef.each)
    })
    truncation <- do.call(intervals::interval_union, truncation)
    truncation <- intervals::interval_union(truncation, intervals::Intervals(c(-Inf, 0))) # consider the implementation of TF_roots, we need to ensure truncation is inside (0, Inf)
    truncation <- intervals::interval_complement(truncation)
    return(truncation)
  }

}


# ----- computing the confidence intervals -----

# return selective ci for v^T mu (i.e. a single parameter)

#' compute selective confidence intervals
#'
#' This function computes the selective confidence intervals.
#'
#' @keywords internal
#'
#' @param v, the contrast vector.
#' @param y, the response.
#' @param sigma, the noise level \eqn{\sigma}.
#' @param truncation, the truncation set for the \eqn{Z}-statistic.
#' @param alpha, the significance level.
#'
#' @return This function returns a vector of lower and upper confidence limits.
#'
computeCI <- function(v, y, sigma, truncation, alpha) {
  #browser()
  vTv <- sum(v*v)
  scale <- sigma * sqrt(vTv)
  q <- sum(v*y) / scale

  fun <- function(x) {
    return(TNSurv(q, x/scale, 1, truncation))
  }

  # L: fun.L(L) = 0
  fun.L <- function(x) {
    return(fun(x) - alpha/2)
  }
  # U: fun.U(U) = 0
  fun.U <- function(x) {
    return(fun(x) - (1-alpha/2))
  }

  # find the starting point (x1, x2) such that
  # fun.L(x1), fun.U(x1) <= 0 AND fun.L(x2), fun.U(x2) >= 0.
  # i.e. fun(x1) <= alpha/2 AND fun(x2) >= 1-alpha/2.

  # find x1 s.t. fun(x1) <= alpha/2
  # what we know:
  # fun is monotone incresing;
  # fun(x) = NaN if x too small;
  # fun(x) > alpha/2 if x too big.
  # so we can do a modified bisection search to find x1.
  step <- 0
  x1.up <- q * scale + scale
  x1 <- q * scale - 10 * scale
  f1 <- fun(x1)
  while(step <= 20) {
    if (is.na(f1)) { # x1 is too small
      x1 <- (x1 + x1.up) / 2
      f1 <- fun(x1)
    }
    else if (f1 > alpha/2) { # x1 is too big
      x1.up <- x1
      x1 <- x1 - 10 * 1.4^step
      f1 <- fun(x1)
      step <- step + 1
    }
    else { # fun(x1) <= alpha/2, excited!
      break
    }
  }

  # find x2 s.t. fun(x2) <= 1 - alpha/2
  # what we know:
  # fun is monotone incresing;
  # fun(x) = NaN if x too big;
  # fun(x) < 1 - alpha/2 if x too small.
  # again can do a modified bisection search to find x2.
  step <- 0
  x2 = q * scale + 10 * scale
  x2.lo = q * scale - scale
  f2 = fun(x2)
  while(step <= 20) {
    if (is.na(f2)) { # x2 is too big
      x2 <- (x2 + x2.lo) / 2
      f2 <- fun(x2)
    }
    else if (f2 < 1 - alpha/2) { # x2 is too small
      x2.lo <- x2
      x2 <- x2 + 10 * 1.4^step
      f2 <- fun(x2)
      step <- step + 1
    }
    else { # fun(x2) >= 1 - alpha/2, excited!
      break
    }
  }



  # count = 0
  # while(is.na(f1)||(f1 > alpha/2)) {
  #   if(count >= 1000) {
  #     break
  #   }
  #   if (is.na(f1)) { # x1 is too small
  #     x1 = x1 + 1*scale
  #   }
  #   else { # fun(x1) > alpha/2
  #     x1 = x1 - 10 * 1.5^count * scale
  #   }
  #   f1 = fun(x1)
  #   count = count+1
  # }
  #
  # x2 = q * scale + 10 * scale
  # f2 = fun(x2)
  #
  # count = 0
  # while(is.na(f2)||(f2 < 1-alpha/2)) {
  #   if(count >= 1000) {
  #     break
  #   }
  #   if (is.na(f2)) { # x1 is too big
  #     x2 = x2 - 1*scale
  #   }
  #   else { # fun(x2) < 1-alpha/2
  #     x2 = x2 + 10*scale
  #   }
  #   f2 = fun(x2)
  #   count = count+1
  # }

  # if the above search does not work, set up a grid search
  # for starting points
  if (is.na(f1)||(f1 > alpha/2)||is.na(f2)||(f2 < 1-alpha/2)) {
    grid <- seq(from = q * scale - 1000*scale, to = q*scale + 1000*scale)
    value <- sapply(grid, fun)
    # want max x1: fun(x1) <= alpha/2
    ind1 <- rev(which(value <= alpha/2))[1]
    x1 <- grid[ind1]
    #f1 = value[ind1]
    # want min x2: fun(x2) >= 1-alpha/2
    ind2 <- which(value >= 1 - alpha/2)[1]
    x2 <- grid[ind2]
    #f2 = value[ind2]
  }

  # if the above fails, then either x1, x2 = NA, so uniroot() will throw error,
  # in which case we set (-Inf, Inf) as the CI

  # we know the functions are increasing

  L <- tryCatch({
    stats::uniroot(fun.L, c(x1, x2), extendInt = "upX", tol = 1e-5)$root
  }, error = function(e) {
    -Inf
  })


  U <- tryCatch({
    stats::uniroot(fun.U, c(x1, x2), extendInt = "upX", tol = 1e-5)$root
  }, error = function(e) {
    Inf
  })

  return(c(L, U))
}

