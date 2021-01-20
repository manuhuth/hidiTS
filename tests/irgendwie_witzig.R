#devtools::install(dependencies = TRUE)

library(hidiTS)
#library("optimParallel")
rm(list = ls())


#parallel computing
cl <- makeCluster(detectCores()-1)

setDefaultCluster(cl = cl)


set.seed(42)
start <- Sys.time()
for (index in 1:50) {

n=10
p=0
q=2
t=20
k=1

  data_test <- sim_data(p = n, T = t, dim_F= q, lags_F=k, lags_X=p, ar_F=1, ar_Y=1, A = NULL, low_X = -3^0.5, up_X = 3^0.5,
                        L = NULL, low_F = 0.2, up_F = 0.9,
                        beta = NULL, gamma = NULL, low_Y = -0.9, up_Y = 0.9,
                        burn_in = 20, data_only = FALSE,
                        only_stationary = TRUE, epsilon = 0.0002, max_it_station = 500, adjust_diag = FALSE,
                        geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)
  clusterExport(cl, list('n', 'p', 'q', 't', 'k'))
  est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
             method = "L-BFGS-B", parallel = FALSE, max_it = 30, trace=0)
  #est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
  #                  method = "L-BFGS-B", parallel = TRUE, max_it = 30, trace=0, forward = TRUE, loginfo=FALSE)
  print(index)
}
end <- Sys.time()
end - start 


normal <- end - start
parallel_normal <- end - start
parallel_forward <- end - start

f <- convert_factors_ML(est$lambda[[1]], est$f_final, q)$F
f_d <- convert_factors_dgp(data_test)$F
plot(f[1,], type = 'l')
lines(-f_d[1,], type='l', col='blue')


