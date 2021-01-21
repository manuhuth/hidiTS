#test if mean sqaured error is strongly influenced by iterations


#devtools::install(dependencies = TRUE)
library(hidiTS)
library("optimParallel")
rm(list = ls())


#set up parallel computing
cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl = cl)


#set variables
max_it_low <- 5
max_it_moderate <- 15
max_it_high <- 40
max_it_very_high <- 100
max_it_very_very_high <- 500

max_its <- c(2, max_it_low, max_it_moderate, max_it_high, max_it_very_high, max_it_very_very_high)

n=10
p=0
q=2
t=20
k=1
clusterExport(cl, list('n', 'p', 'q', 't', 'k'))

export <- c()
for (index2 in 1:1000) {

  #set seed
  s <- index2*2
  set.seed(index2*2)

  data_test <- sim_data(p = n, T = t, dim_F= q, lags_F=k, lags_X=p, ar_F=1, ar_Y=1, A = NULL, low_X = -3^0.5, up_X = 3^0.5,
                        L = NULL, low_F = 0.2, up_F = 0.9,
                        beta = NULL, gamma = NULL, low_Y = -0.9, up_Y = 0.9,
                        burn_in = 20, data_only = FALSE,
                        only_stationary = TRUE, epsilon = 0.0002, max_it_station = 500, adjust_diag = FALSE,
                        geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

  for (index in max_its) {
    start <- Sys.time()
    est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
                      method = "L-BFGS-B", parallel = TRUE, max_it = index, trace=0, forward = TRUE, loginfo=FALSE)
    end <- Sys.time()
    f <- convert_factors_ML(est$lambda[[1]], est$f_final, q)$F[1,]
    f_d <- convert_factors_dgp(data_test)$F[1,]
    time <- end - start
    mse <- mean( (f - f_d)^2 )
    mse2 <- mean( (f + f_d)^2 )
    if (mse2 < mse ){
        mse <- mse2
      }

      to_save <- c('max_it'=index, 'mse_estimated'= mse, 'time'=time,'data_seed' = s ,'n'=n, 'p'=p, 'q'=2, 't'=20, 'k'=1, 'parallel'='yes')
      export <- rbind(export, to_save)
    }
  print(index2)

}
