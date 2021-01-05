pca_estimator <- function(X, number) {
  sigma_X <- var(t(X))
  Lambda_hat <- eigen(sigma_X)$vectors[,1:number]
  F_hat <-  t(Lambda_hat) %*% X * 1/ncol(X)

  output <- list('F' = F_hat, 'Lambda' = matrix(Lambda_hat, nrow(X), number))
  return(output)
}



convert_pca_estimator <- function(pca_est, q){
  F <- pca_est$F
  qp1 <- nrow(F)
  lags <- qp1/q-1
  Lambda <- pca_est$Lambda

  V <- var(t(F_stack))
  theta <- eigen(V)$vectors
  D <- diag(1/eigen(V)$values^0.5, qp1)
  F_unc <- D %*% t(theta) %*% F_stack
  Lambda_unc <- Lambda %*% solve(D %*% theta)
  QR_2 <- qr(t(Lambda_unc[1:qp1, 1:qp1]))
  L_qr_2 <- qr.Q(QR_2)

  new_F_hat <- L_qr_2%*%F_unc
  Lambda_hat_new <- Lambda_unc %*% L_qr_2

  if (lags == 0) {
    F_out <- new_F_hat
  } else{
    F_out <- cbind(new_F_hat[1:q, 1:(T-lags-1)], matrix(new_F_hat[1:qp1,T-lags], q, lags+1) )
  }

  output <- list('F' = F_out, 'Lambda' = Lambda_hat_new)

  return(output)
}
