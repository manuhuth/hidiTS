stack_F_matrix <- function(F, lags, start_ar) {
  par <- lags + 1
  T <- ncol(F)
  big_F <- F[,(start_ar + 1):T]
  if (lags != 0) {
    for (j in 1:(par - 1)) {
      big_F <- rbind(big_F, F[,(start_ar + 1 - j):(T - j)])
    }

    help_trans <- c()
    for (index in 1:nrow(big_F)){
      help_trans <- rbind(help_trans, big_F[nrow(big_F) - index +1,])
    }

    big_F <- help_trans
  }

  return(big_F)
}

convert_factors_dgp <- function(data) {
  L <- do.call(cbind, data$L)
  qp1 <- ncol(L)
  q <- ncol(data$L[[1]])
  lags <- qp1/q-1
  F_dgp <- stack_F_matrix(data$F, lags, lags)
  T <- ncol(data$F)


  V <- var(t(F_dgp))
  theta <- eigen(V)$vectors
  D <- diag(1/eigen(V)$values^0.5, qp1)
  F_unc <- D %*% t(theta) %*% F_dgp


  L_unc <- L %*% solve(D %*% t(theta))
  QR <- qr(t(L_unc[1:qp1, 1:qp1]))
  L_qr <- qr.Q(QR)
  new_F <- L_qr%*%F_unc
  new_L <- L_unc %*% L_qr

  if (lags == 0) {
    F_out <- new_F
  } else{
    F_out <- cbind(new_F[1:q, 1:(T-lags-1)], matrix(new_F[1:qp1,T-lags], q, lags+1) )
  }

  output <- list('F' = F_out, 'Lambda' = new_L)
  return(output)
}


convert_factors_ML <- function(Lambda, factors, q) {
  L <- Lambda[[1]]
  
  qp1 <- ncol(L)
  lags <- qp1/q-1
  F_ML <- do.call(cbind, factors)[1:q,]
  T <- ncol(F_ML)
  
  help_stack <- c()
  for (index in 1:(lags+1)) {
    help_stack <- cbind(help_stack, factors[[1]][(qp1-q*(index)+1):(qp1-q*(index-1))])
  }
  
  F_temp <- cbind(help_stack, F_ML[,2:T])
  
  F_stacked <- F_temp[,(lags+1):(T+lags)]
  
  for(index in 1:(lags)){
    F_stacked <- rbind(F_stacked,F_temp[,(lags+1-index):(T+lags-index)])
  }
  
  V <- var(t(F_stacked))
  theta <- eigen(V)$vectors
  D <- diag(1/eigen(V)$values^0.5, qp1)
  F_unc <- D %*% t(theta) %*% F_stacked
  
  
  L_unc <- L %*% solve(D %*% t(theta))
  QR <- qr(t(L_unc[1:qp1, 1:qp1]))
  L_qr <- qr.Q(QR)
  
  
  new_F <- L_qr%*%F_unc
  new_L <- L_unc %*% L_qr
  
  if (lags == 0) {
    F_out <- new_F
  } else{
    F_out <- cbind(new_F[1:q, 1:(T-lags-1)], matrix(new_F[1:qp1,T-lags], q, lags+1) )
  }
  
  output <- list('F' = F_out, 'Lambda' = new_L)
  return(output)
}
