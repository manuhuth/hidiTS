convert_factors_dgp <- function(data) {
  L <- data$L[[1]]
  F <- data$F
  QR <- qr(t(L[1:q, 1:q]))
  L_qr <- qr.Q(QR)
  new_F <- L_qr%*%F
  new_L <- L %*% t(L_qr)

  output <- list('F' = new_F, 'Lambda' = new_L)
  return(output)
}

