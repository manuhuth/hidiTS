convert_factors_dgp <- function(data) {
  F_dgp <- data$F
  L <- data$L[[1]]
  q <- ncol(L)

  V <- var(t(F_dgp))
  theta <- eigen(V)$vectors
  D <- diag(1/eigen(V)$values^0.5, q)
  F_unc <- D %*% t(theta) %*% F_dgp


  L_unc <- L %*% solve(D %*% t(theta))
  QR <- qr(t(L_unc[1:q, 1:q]))
  L_qr <- qr.Q(QR)
  new_F <- L_qr%*%F_unc
  new_L <- L_unc %*% t(L_qr)

  output <- list('F' = new_F, 'Lambda' = new_L)
  return(output)
}

