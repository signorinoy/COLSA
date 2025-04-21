#' Objective Function for Cox Model Optimization
#'
#' Computes the objective value, gradient, and Hessian for optimizing a Cox
#' proportional hazards model with spline-based baseline hazard estimation.
#'
#' @param par Numeric vector of parameters to be optimized, including spline
#'        coefficients for the baseline hazard (`alpha`) and regression
#'        coefficients (`beta`).
#' @param time Numeric vector of observed survival times.
#' @param status Numeric vector of event indicators.
#' @param x Numeric matrix of covariates (features) for the model.
#' @param boundary Numeric vector specifying the boundary knots for the splines.
#' @param theta Numeric vector of prior parameters for regularization.
#' @param hessian Numeric matrix of prior Hessian for regularization.
#'
#' @return A list containing:
#'   \item{value}{Scalar value of the objective function.}
#'   \item{gradient}{Numeric vector of the gradient of the objective function.}
#'   \item{hessian}{Numeric matrix of the Hessian of the objective function.}
#'
#' @details
#' The objective function combines two components:
#' - A regularization term based on prior parameters (`theta`) and prior Hessian
#'   (`hessian`).
#' - A likelihood-based term for the Cox proportional hazards model, where the
#'   baseline hazard is modeled using Bernstein polynomials.
#'
#' The function calculates the loss, gradient, and Hessian for both components
#' and combines them to return the overall objective value, gradient, and
#' Hessian. The baseline hazard is integrated using helper functions
#' (`int_basehaz`) to compute cumulative hazard and its derivatives.
objective <- function(par, time, status, x, boundary, theta, hessian) {
  n_parameters <- length(par)
  n_features <- ncol(x)
  n_basis <- n_parameters - n_features

  n_parameters_pre <- length(theta)
  n_basis_pre <- n_parameters_pre - n_features

  # Transform the parameters and Hessian matrix from n_basis_pre to n_basis
  prox <- prox_reverse(n_basis_pre, n_basis)
  prox <- rbind(
    cbind(prox, matrix(0, nrow(prox), n_features)),
    cbind(matrix(0, n_features, ncol(prox)), diag(n_features))
  )
  theta <- as.vector(prox %*% theta)
  hessian <- prox %*% hessian %*% t(prox)

  loss1 <- as.numeric((par - theta) %*% hessian %*% (par - theta) / 2)
  grad1 <- as.vector(hessian %*% (par - theta))
  hess1 <- hessian

  alpha <- par[seq_len(n_basis)]
  beta <- par[(n_basis + 1):n_parameters]
  b <- if (n_basis == 1) {
    matrix(1, nrow = length(time), ncol = 1)
  } else {
    splines2::bpoly(
      time,
      degree = n_basis - 1, intercept = TRUE,
      Boundary.knots = boundary
    )
  }
  eta <- x %*% beta

  cbh <- do.call(rbind, lapply(time, function(t) {
    int_basehaz(t, 0, alpha, n_basis, boundary)
  }))
  dcbh <- do.call(rbind, lapply(time, function(t) {
    int_basehaz(t, 1, alpha, n_basis, boundary)
  }))
  d2cbh <- t(sapply(time, function(t) {
    int_basehaz(t, 2, alpha, n_basis, boundary)
  }))
  if (n_basis == 1) {
    d2cbh <- matrix(d2cbh, ncol = 1)
  }

  dalpha <- t(dcbh) %*% exp(eta) - t(b) %*% status
  dbeta <- t(x) %*% (cbh * exp(eta) - status)
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), n_basis)
  d2beta <- t(x) %*% (as.vector(cbh * exp(eta)) * x)
  dalph_dabeta <- t(dcbh) %*% (as.vector(exp(eta)) * x)

  loss2 <- sum(cbh * exp(eta) - status * (b %*% alpha + eta))
  grad2 <- c(dalpha, dbeta)
  hess2 <- cbind(
    rbind(d2alpha, t(dalph_dabeta)),
    rbind(dalph_dabeta, d2beta)
  )

  list(
    value = loss1 + loss2,
    gradient = grad1 + grad2,
    hessian = hess1 + hess2
  )
}

#' Compute Integrated Baseline Hazard and Its Derivatives
#'
#' This function computes the integrated baseline hazard or its derivatives
#' using Gaussian quadrature and Bernstein polynomial basis functions.
#'
#' @param t Numeric. The upper limit of integration.
#' @param deriv Integer. Specifies the derivative to compute:
#'   \itemize{
#'     \item 0: Computes the cumulative baseline hazard.
#'     \item 1: Computes the gradient of the cumulative baseline hazard.
#'     \item 2: Computes the Hessian of the cumulative baseline hazard.
#'   }
#' @param alpha Numeric vector. Coefficients for the baseline hazard function.
#' @param n_basis Integer. The number of Bernstein polynomial basis functions.
#' @param boundary Numeric vector of length 2. The boundary knots for the
#' Bernstein polynomials.
#' @param n_nodes Integer. The number of quadrature nodes to use. Default is 10.
#'
#' @return Numeric. The computed integral or its derivative.
#'
#' @details
#' This function uses Gaussian quadrature to approximate the integral of the
#' baseline hazard function, which is modeled using Bernstein polynomial basis
#' functions. The `deriv` parameter determines whether the function computes
#' the integral itself or its first or second derivative.
#'
#' @importFrom statmod gauss.quad.prob
int_basehaz <- function(t, deriv, alpha, n_basis, boundary, n_nodes = 10) {
  if (!is.numeric(t) || t < 0) stop("t must be a non-negative numeric value")
  if (!is.numeric(deriv) || !(deriv %in% 0:2)) stop("deriv must be 0, 1, or 2")
  if (!is.numeric(alpha) || length(alpha) != n_basis) {
    stop("alpha must be a numeric vector of length n_basis")
  }
  if (
    !is.numeric(boundary) || length(boundary) != 2 || boundary[1] >= boundary[2]
  ) {
    stop("boundary must be a ordered numeric vector of length 2")
  }
  if (!is.numeric(n_nodes) || n_nodes <= 0 || n_nodes != as.integer(n_nodes)) {
    stop("n_nodes must be a positive integer")
  }

  quad <- gauss.quad.prob(n_nodes, "uniform", l = 0, u = t)
  b <- if (n_basis == 1) {
    matrix(1, nrow = n_nodes, ncol = 1)
  } else {
    splines2::bpoly(
      quad$nodes,
      degree = n_basis - 1, intercept = TRUE,
      Boundary.knots = boundary
    )
  }
  bh <- as.vector(exp(b %*% alpha))
  if (n_basis == 1) {
    return(t * sum(quad$weights * bh))
  }
  if (deriv == 0) {
    int <- sum(quad$weights * bh)
  } else if (deriv == 1) {
    int <- colSums(sweep(b, 1, quad$weights * bh, "*"))
  } else if (deriv == 2) {
    b_outer <- t(sapply(seq_len(n_nodes), function(i) {
      crossprod(b[i, , drop = FALSE])
    }))
    int <- colSums(sweep(b_outer, 1, quad$weights * bh, "*"))
  }
  t * int
}

#' Proximal Operator (Forward)
#'
#' Constructs a transformation matrix to map the coefficients of Bernstein
#' polynomial basis functions from one dimension to another, ensuring
#' consistency with the properties of Bernstein polynomial bases. Note that
#' the original dimension must be less than or equal to the target dimension.
#'
#'
#' @param p Integer. The number of basis functions in the original dimension.
#' @param q Integer. The number of basis functions in the target dimension.
#'
#' @return A matrix representing the transformation from the original dimension
#'         to the target dimension.
#'
#' @details This function iteratively constructs the transformation matrix by
#'          expanding the basis functions from \code{p} to \code{q}, while
#'          preserving the mathematical properties of Bernstein polynomials.
#'
#' @seealso \code{\link{prox_reverse}}
prox_forward <- function(p, q) {
  if (
    !is.numeric(p) || !is.numeric(q) || p != as.integer(p) || q != as.integer(q)
  ) {
    stop("p and q must be integers")
  }
  if (p <= 0 || q <= 0) stop("p and q must be positive integers")
  if (p > q) stop("p must be less than or equal to q")
  if (p == q) {
    return(diag(p))
  }
  prox_forward <- diag(p)
  for (i in p:(q - 1)) {
    if (i == 1) {
      prox <- matrix(1, 2, 1)
    } else {
      prox <- matrix(0, i + 1, i)
      diag(prox) <- (i - 0:(i - 1)) / i
      diag(prox[-1, ]) <- (1:i) / i
    }
    prox_forward <- prox %*% prox_forward
  }
  prox_forward
}

#' Proximal Operator (Reverse)
#'
#' Constructs a transformation matrix to map the coefficients of Bernstein
#' polynomial basis functions from one dimension to another, ensuring
#' consistency with the properties of Bernstein polynomial bases. Note that
#' the original dimension must be greater than or equal to the target dimension.
#'
#' @param q Integer. The number of basis functions in the original dimension.
#' @param p Integer. The number of basis functions in the target dimension.
#'
#' @return A matrix representing the transformation from the original dimension
#'         to the target dimension.
#'
#' @details This function iteratively constructs the transformation matrix by
#'          expanding the basis functions from \code{p} to \code{q}, while
#'          preserving the mathematical properties of Bernstein polynomials.
#'
#' @seealso \code{\link{prox_forward}}
#'
#' @importFrom MASS ginv
prox_reverse <- function(q, p) {
  if (
    !is.numeric(p) || !is.numeric(q) || p != as.integer(p) || q != as.integer(q)
  ) {
    stop("p and q must be integers")
  }
  if (p <= 0 || q <= 0) stop("p and q must be positive integers")
  if (p > q) stop("p must be less than or equal to q")
  if (p == q) {
    return(diag(p))
  }
  prox_forward <- prox_forward(p, q)
  ginv(prox_forward)
}
