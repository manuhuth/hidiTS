pca_estimator <- function(X, number) {
  sigma_X <- var(t(X))
  Lambda_hat <- eigen(sigma_X)$vectors[,1:number]
  F_hat <-  t(Lambda_hat) %*% X * 1/ncol(X)

  output <- list('F' = F_hat, 'Lambda' = matrix(Lambda_hat, nrow(X), number))
  return(output)
}



convert_pca_estimator <- function(pca_est){
  F <- pca_est$F
  q <- nrow(F)
  Lambda <- pca_est$Lambda
  V <- var(t(F))
  theta <- eigen(V)$vectors
  D <- diag(1/eigen(V)$values^0.5, q)
  F_unc <- D %*% t(theta) %*% F
  Lambda_unc <- Lambda %*% solve(D %*% theta)
  QR_2 <- qr(t(Lambda_unc[1:q, 1:q]))
  L_qr_2 <- qr.Q(QR_2)

  new_F_hat <- L_qr_2%*%F_unc
  Lambda_hat_new <- Lambda_unc %*% t(L_qr_2)

  output <- list('F' = new_F_hat, 'Lambda' = Lambda_hat_new)

  return(output)
}
