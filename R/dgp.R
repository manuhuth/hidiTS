#' Create random matrices based on a uniform distribution.
#'
#' @param n Number of matrices to be created.
#' @param row Number of rows per matrix.
#' @param col Number of columns per matrix.
#' @param low lower bound of the uniform distribution.
#' @param up upper bound of the uniform distribution.
#' @return List including the n matrices.
#' @examples
#' random_matrices(n=4, row=5,col=6)
#' @export
random_matrices <- function(n, row, col, low = -0.2, up = 0.2) {
  #to-do: add more flexible ways of matrices (more structures for coefficients)
  A <- list()
  for (j in 1:n) {
    A[[j]] <- matrix(runif(row * col, low, up), row, col)
  }
  return(A)
}



#' Create factors \eqn{F_t} for \eqn{t=1, \dots, T}
#' @param dim Dimension of \eqn{F_t}.
#' @param lags Number of lags to be included for \eqn{F_t}.
#' @param T Length of the time series.
#' @param A List containing the \eqn{A_j} matrices. If NULL, \eqn{A_j} matrices are simulated.
#' @param vcv_eta Variance-covariance matrix of the errors.
#' @return List of time series \eqn{F_t} for \eqn{t=1, \dots, T}.
#' @examples
#' sim_F(dim = 3, lags = 4, T = 30)
sim_F <- function(dim, lags, T, A = NULL, vcv_eta = diag(dim)) {
  # dependent on MASS
  # to-do: change A creation -> see random_matrices
  #        allow for correlated errors within a time series?

  # simulate etas
  eta <- mvrnorm(n = T + lags, mu = rep(0, dim), Sigma = vcv_eta)

  # simulate A
  if (is.null(A)) {
    A <- random_matrices(n = lags, row = dim, col = dim)
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
#' stack_F(F=F,lags=3)
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
#' sim_X(F=F,p=3,lags=2)
sim_X <- function(F, p, lags, L = NULL, vcv_mu = diag(p)) {
  # to-do: change L creation -> see random_matrices
  #        allow for correlated errors within a time series?


  num_lambda <- lags + 1
  dim_F <- length(F[[1]])
  if (is.null(L)) {
    L <- random_matrices(n=num_lambda, row = p, col = dim_F)
  }

  mu <- t(mvrnorm(n = length(F)-lags, mu= rep(0,p), Sigma=vcv_mu ))

  big_F <- stack_F(F, lags = lags)
  L_stack <- do.call(cbind, L)

  X <- L_stack %*% big_F + mu
  output <- list("X"=X, "L"=L, "mu"=mu)
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
#' sim_Y(F=F, ar_F=3, ar_Y=2)
sim_Y <- function(F, ar_F, ar_Y, beta=NULL, gamma=NULL, vcv_epsilon=diag(length(F))){
  # to-do: change beta/gamma creation -> see random_matrices

  dim_F <- length(F[[1]])
  T <- length(F)
  if (is.null(beta)) {
    beta <- random_matrices(n=ar_F, row = 1, col = dim_F)
  }

  if (is.null(gamma)) {
    gamma <- random_matrices(n=1, row = 1, col = ar_Y)
  }

  lags <- max(ar_F, ar_Y)

  epsilon <- mvrnorm(n=1, mu=rep(0,T), Sigma=vcv_epsilon)
  Y <- runif(lags)

  for (t in (lags+1):T) {
    Y[t] <- do.call(cbind, beta) %*% c(do.call(cbind, F[(t-ar_F):(t-1)])) +
      do.call(cbind, gamma) %*% Y[(t-ar_Y):(t-1)] + epsilon[t]
  }

  output <- list("Y"=Y, "beta"=beta, "gamma"=gamma, "epsilon"=epsilon)
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
#' @return List of factors, observable variables and the variable which is to be predicted (errors, loadings, parameter matrices)
#' @examples
#' F <- sim_F(dim = 3, lags = 4, T = 30)
#' sim_Y(F=F, ar_F=3, ar_Y=2)
sim_data <- function(p, T, dim_F, lags_F, lags_X, ar_F, ar_Y, A = NULL, L = NULL,
                     beta = NULL, gamma=NULL, vcv_eta = diag(dim_F), vcv_mu=diag(p),
                     vcv_epsilon=diag(length(F)), burn_in =20 , data_only = TRUE) {

  max_ar <- max(lags_F, lags_X, ar_F, ar_Y)
  sim_length <- T + max_ar + burn_in

  F_object <- sim_F(dim=dim_F,T=sim_length, lags=lags_F, A=A, vcv_eta = vcv_eta)
  F <- F_object$F
  eta <- F_object$eta
  X_object <-  sim_X(F=F,p=p,lags=lags_X, L=L, vcv_mu = vcv_mu)
  X <- X_object$X
  mu <- X_object$mu

  Y_object <- sim_Y(F=F, ar_F=ar_F, ar_Y=ar_Y, beta = beta, gamma=gamma, vcv_epsilon = vcv_epsilon)
  Y <- Y_object$Y
  epsilon <- Y_object$epsilon

  F <- do.call(cbind, F)
  X <- do.call(cbind, X)
  Y_out <- Y[(length(Y)-T+1):length(Y)]
  X_out <- X[,(ncol(X)-T+1):ncol(X)]
  F_out <- F[,(ncol(F)-T+1):ncol(F)]
  eta_out <- eta[(length(eta)-T+1):length(eta)]
  mu_out <- mu[(length(mu)-T+1):length(mu)]
  epsilon_out <- epsilon[(length(epsilon)-T+1):length(epsilon)]

  if (isTRUE(data_only)) {
    output <- list("Y"=Y_out, "X"=X_out, "F"=F_out)
  } else{
    output <- list("Y"=Y_out, "X"=X_ot, "F"=F_out, "beta"=Y_object$beta, "gamma"=Y_object$gamma,
                   "A"=F_object$A, "L"=X_object$L, "eta"=eta_out, "mu"=mu_out,
                   "epsilon"=epsilon_out)
  }
  return(output)
}

