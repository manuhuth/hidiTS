#' Create random matrices based on a uniform distribution.
#'
#' @param n Number of matrices to be created.
#' @param row Number of rows per matrix.
#' @param col Number of columns per matrix.
#' @param low lower bound of the uniform distribution.
#' @param up upper bound of the uniform distribution.
#' @return List including the n matrices.
#' @examples
#' random_matrices(n = 4, row = 5, col = 6)
#' @export
random_matrices <- function(n, row, col, low = -0.5, up = 0.5) {
  # to-do: add more flexible ways of matrices (more structures for coefficients)
  A <- list()
  for (j in 1:n) {
    A[[j]] <- matrix(runif(row * col, low, up), row, col)
  }
  return(A)
}


#' Check if matrices created by random_matrices generate a covariance stationary AR process..
#'
#' @param gamma list created by random_matrices
#' @return Boolean indicating if matrices yield covariance stationary process.
#' @examples
#' gamma <- random_matrices(n = 4, row = 5, col = 5)
#' cov_stationarity(gamma)
#' @export
cov_stationarity <- function(gamma, epsilon=0.0002) {
  Gamma <- do.call(cbind, gamma)
  k <- length(gamma)
  q <- ncol(gamma[[1]])
  lower_diag <- diag((k - 1) * q) # create lower block matrix
  zero_right <- matrix(rep(0, len = (k - 1) * q^2), ncol = q, nrow = (k - 1) * q)
  lower <- cbind(lower_diag, zero_right)
  M <- rbind(Gamma, lower)
  eigenvalues <- eigen(M)$values
  output <- all(abs(eigenvalues) < 1-epsilon)

  return(output)
}


#' Create factors \eqn{F_t} for \eqn{t=1, \dots, T}
#' @param dim Dimension of \eqn{F_t}.
#' @param lags Number of lags to be included for \eqn{F_t}.
#' @param T Length of the time series.
#' @param A List containing the \eqn{A_j} matrices. If NULL, \eqn{A_j} matrices are simulated.
#' @param vcv_eta Variance-covariance matrix of the errors.
#' @param only_stationary boolean indicating if only stationary TS should be returned.
#' @param epsilon Tolerance parameter that evaluates if roots of polynomial are inside unit circle. (1-epsilon)
#' @param max_it_station Integer specifying the maximal number of iterations to find a stationary process.
#' @return List of time series \eqn{F_t} for \eqn{t=1, \dots, T}.
#' @examples
#' sim_F(dim = 3, lags = 4, T = 30)
#' @export
sim_F <- function(dim, lags, T, A = NULL, low = -0.4, up = 0.4, vcv_eta = diag(dim),
                  only_stationary = TRUE, epsilon=0.0002, max_it_station = 500) {
  # dependent on MASS
  # to-do: change A creation -> see random_matrices
  #        allow for correlated errors within a time series?

  # simulate etas
  eta <- mvrnorm(n = T + lags, mu = rep(0, dim), Sigma = vcv_eta)

  # simulate A
  if (is.null(A)) {
    A <- random_matrices(n = lags, row = dim, col = dim, low = low, up = up)
    stationar <- cov_stationarity(A, epsilon = epsilon)
    if (isTRUE(only_stationary)) {
      it <- 0
      while (isFALSE(stationar)) {
        A <- random_matrices(n = lags, row = dim, col = dim, low = low, up = up)
        stationar <- cov_stationarity(A, epsilon = epsilon)
        it <- it + 1
        if (it > max_it_station) {
          stop("No stationary process was found to build the factors.")
        }
      }
    }
  }

  # simulate first and second Fs
  F <- list(eta[1, ])

  # simulate second Fs
  F[[2]] <- A[[1]] %*% F[[1]] + eta[2, ]
  # simulate until burn in is over
  if (lags > 2) {
    for (j in 3:lags) {
      F[[j]] <- do.call(cbind, A[1:j - 1]) %*% c(do.call(cbind, F)) + eta[j, ] # fast way to compute F at period t using some algebra
    }
  }
  # simulate until burn in is over
  A_stack <- do.call(cbind, A)
  for (t in 1:T) {
    F[[t + lags]] <- A_stack %*% c(do.call(cbind, F[t:(t + lags)])) + eta[t, ]
  }

  F_ts <- F[(lags + 1):(T + lags)]
  F_burn <- F[(1):(lags)]

  output <- list("F" = F_ts, "F_burn_in" = F_burn, "A" = A, "eta" = t(eta))
  return(output)
}


#' Creates matrix of lagged values of factors.
#'
#' @param F List of Factors.
#' @param lags Number of lags to be included to compute the observable variables.
#' @return Matrix containing the lagged factor values in the proper shape to multiply it with the loadings.
#' @examples
#' F <- sim_F(dim = 3, lags = 4, T = 30)
#' stack_F(F = F, lags = 3)
#' @export
stack_F <- function(F, lags) {
  par <- lags + 1
  T <- length(F)
  big_F <- do.call(cbind, F[par:T])
  for (j in 1:(par - 1)) {
    big_F <- rbind(big_F, do.call(cbind, F[(par - j):(T - j)]))
  }
  return(big_F)
}


#' Create observable variables X.
#'
#' @param F List of factors.
#' @param p Number of observable variables per period.
#' @param lags Number of lags to be included for \eqn{F_t}.
#' @param L List containing the \eqn{L_j} matrices (loadings). If NULL, \eqn{L_j} matrices are simulated.
#' @param vcv_mu Variance-covariance matrix of the errors.
#' @return List of time series \eqn{X_t} for \eqn{t=1, \dots, T}, loadings and errors.
#' @examples
#' sim_F(dim = 3, lags = 4, T = 30)
#' sim_X(F = F, p = 3, lags = 2)
#' @export
sim_X <- function(F, p, lags, L = NULL, low = -0.4, up = 0.4, vcv_mu = diag(p)) {
  # to-do: change L creation -> see random_matrices
  #        allow for correlated errors within a time series?


  num_lambda <- lags + 1
  dim_F <- length(F[[1]])
  if (is.null(L)) {
    L <- random_matrices(n = num_lambda, row = p, col = dim_F, low = low, up = up)
  }

  mu <- t(mvrnorm(n = length(F) - lags, mu = rep(0, p), Sigma = vcv_mu))

  big_F <- stack_F(F, lags = lags) # stack F such that X can be computed with one matrix multiplication
  L_stack <- do.call(cbind, L)

  X <- L_stack %*% big_F + mu
  output <- list("X" = X, "L" = L, "mu" = mu)
  return(output)
}

#' Create Y variable from factors.
#'
#' @param F List of factors.
#' @param ar_F Number of lagged factors used to compute Y.
#' @param ar_Y Number of lagged Y used to compute Y.
#' @param beta Parameter for the Factors. If NULL, \eqn{\beta_j} matrices are simulated.
#' @param gamma Parameter for the lagged Y's. If NULL, \eqn{\gamma_j} vector is simulated.
#' @param vcv_epsilon Variance-covariance matrix of the errors.
#' @return List of time series \eqn{Y_t} for \eqn{t=1, \dots, T}, parameters and errors.
#' @examples
#' F <- sim_F(dim = 3, lags = 4, T = 30)
#' sim_Y(F = F, ar_F = 3, ar_Y = 2)
#' @export
sim_Y <- function(F, ar_F, ar_Y, beta = NULL, gamma = NULL, low = -0.4, up = 0.4, vcv_epsilon = diag(length(F))) {
  # to-do: change beta/gamma creation -> see random_matrices

  dim_F <- length(F[[1]])
  T <- length(F)
  if (is.null(beta)) {
    beta <- random_matrices(n = ar_F, row = 1, col = dim_F, low = low, up = up)
  }

  if (is.null(gamma)) {
    gamma <- random_matrices(n = 1, row = 1, col = ar_Y, low = low, up = up)
  }

  lags <- max(ar_F, ar_Y)

  epsilon <- mvrnorm(n = 1, mu = rep(0, T), Sigma = vcv_epsilon)

  Y <- c(0) # necessary to initialize with 0. Is overwritten in first loop
  for (t in 1:ar_Y) {
    Y[t] <- do.call(cbind, beta) %*% c(do.call(cbind, F[(t):(t + ar_F - 1)])) + gamma[[1]][, 1:length(Y)] %*% Y + rnorm(n = 1, 0, 1)
  }

  for (t in (lags + 1):T) {
    Y[t] <- do.call(cbind, beta) %*% c(do.call(cbind, F[(t - ar_F):(t - 1)])) +
      do.call(cbind, gamma) %*% Y[(t - ar_Y):(t - 1)] + epsilon[t]
  }

  output <- list("Y" = Y, "beta" = beta, "gamma" = gamma, "epsilon" = epsilon)
  return(output)
}


#' Create whole data set.
#'
#' @param p Number of observable variables per period. Passed to \code{sim_F}
#' @param T Length of the time series. Passed to
#' @param dim_F Dimension of the factors in each period.
#' @param ar_F Number of lagged factors used to compute Y.
#' @param ar_Y Number of lagged Y used to compute Y.
#' @param A Parameters for factors in F. If NULL, \eqn{A_j} matrices are simulated.
#' @param L Factor loadings.  If NULL, \eqn{L_j} matrices are simulated.
#' @param beta Parameter for the Factors in Y. If NULL, \eqn{\beta_j} matrices are simulated.
#' @param gamma Parameter for the lagged Y's. If NULL, vector is simulated.
#' @param vcv_eta Variance-covariance matrix of the errors to create factors.
#' @param vcv_mu Variance-covariance matrix of the errors to create observable variables X.
#' @param vcv_epsilon Variance-covariance matrix of the errors to create variable to forecast Y.
#' @param burn_in Length of burn-in period that is simulated in addition to receive a time-series more stable on the data.
#' @param data_only If TRUE, only the values for factors, observable variables and the variable which is to be predicted are returned.
#' @param only_stationary boolean indicating if only stationary TS should be returned.
#' @param epsilon Tolerance parameter that evaluates if roots of polynomial are inside unit circle. (1-epsilon)
#' @param max_it_station Integer specifying the maximal number of iterations to find a stationary process.
#' @return List of factors, observable variables and the variable which is to be predicted (errors, loadings, parameter matrices)
#' @examples
#' F <- sim_F(dim = 3, lags = 4, T = 30)
#' sim_Y(F = F, ar_F = 3, ar_Y = 2)
#' @export
sim_data <- function(p, T, dim_F, lags_F, lags_X, ar_F, ar_Y, A = NULL, low_X = -0.5, up_X = 0.5,
                     L = NULL, low_F = -0.5, up_F = 0.5,
                     beta = NULL, gamma = NULL, low_Y = -0.5, up_Y = 0.5, vcv_eta = diag(dim_F), vcv_mu = diag(p),
                     vcv_epsilon = diag(length(F)), burn_in = 20, data_only = TRUE,
                     only_stationary = TRUE, epsilon=0.0002, max_it_station = 500) {
  max_ar <- lags_F + lags_X + ar_F + ar_Y
  sim_length <- T + max_ar + burn_in

  F_object <- sim_F(dim = dim_F, T = sim_length, lags = lags_F, A = A, vcv_eta = vcv_eta, low = low_F, up = up_F,
                    only_stationary = only_stationary, epsilon=epsilon, max_it_station = max_it_station)
  F <- F_object$F
  eta <- F_object$eta
  X_object <- sim_X(F = F, p = p, lags = lags_X, L = L, vcv_mu = vcv_mu, low = low_X, up = up_X)
  X <- X_object$X
  mu <- X_object$mu

  Y_object <- sim_Y(
    F = F, ar_F = ar_F, ar_Y = ar_Y, beta = beta, gamma = gamma, vcv_epsilon = vcv_epsilon,
    low = low_Y, up = up_Y
  )
  Y <- Y_object$Y
  epsilon <- Y_object$epsilon

  F <- do.call(cbind, F)
  Y_out <- Y[(length(Y) - T + 1):length(Y)]
  X_out <- X[, (ncol(X) - T + 1):ncol(X)]
  F_out <- F[, (ncol(F) - T + 1):ncol(F)]
  eta_out <- eta[(length(eta) - T + 1):length(eta)]
  mu_out <- mu[(length(mu) - T + 1):length(mu)]
  epsilon_out <- epsilon[(length(epsilon) - T + 1):length(epsilon)]

  if (isTRUE(data_only)) {
    output <- list("Y" = Y_out, "X" = X_out, "F" = F_out)
  } else {
    output <- list(
      "Y" = Y_out, "X" = X_ot, "F" = F_out, "beta" = Y_object$beta, "gamma" = Y_object$gamma,
      "A" = F_object$A, "L" = X_object$L, "eta" = eta_out, "mu" = mu_out,
      "epsilon" = epsilon_out
    )
  }
  return(output)
}
