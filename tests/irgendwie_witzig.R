#devtools::install(dependencies = TRUE)

library(hidiTS)
rm(list = ls())



#parallel computing

for (index in 1:10) {
n=10
p=2
q=3
t=30
k=1

  data_test <- sim_data(p = n, T = t, dim_F= q, lags_F=k, lags_X=p, ar_F=1, ar_Y=1, A = NULL, low_X = -3^0.5, up_X = 3^0.5,
                        L = NULL, low_F = 0.2, up_F = 0.9,
                        beta = NULL, gamma = NULL, low_Y = -0.9, up_Y = 0.9,
                        burn_in = 20, data_only = FALSE,
                        only_stationary = TRUE, epsilon = 0.0002, max_it_station = 500, adjust_diag = FALSE,
                        geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)
  
  est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
             method = "L-BFGS-B", parallel = FALSE, max_it = 5)
  print(index)
}



its <- c(50)
est <- list()
for (index in 1:length(its)) {
  est[[index]] <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=its[index],
                     method = "L-BFGS-B", parallel = FALSE, max_it = 3)
  print(its[index])

}

matrix <- c()
#convert_ML
for (index in 1:length(its)) {
  val <- convert_factors_ML(Lambda=est[[index]]$lambda[[1]], factors = est[[index]][['f_final']], q=q)$F[1,] 
  if (val[1] > 0) {
    val <- -val
  }
  matrix <- cbind(matrix, val)
}

#convert data
dgp <- convert_factors_dgp(data_test)

matplot(matrix, type = 'l')
lines(dgp$F[1,], type='l', col='red')