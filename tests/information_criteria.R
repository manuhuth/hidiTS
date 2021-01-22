#devtools::install(dependencies = TRUE)
library(hidiTS)
library(optimParallel)


rm(list = ls())
set.seed(123)

cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl = cl)

n=50
p=1
q=2
t=40
k=1

data_test <- sim_data(p = n, T = t, dim_F= q, lags_F=k, lags_X=p, ar_F=1, ar_Y=1, A = NULL, low_X = -3^0.5, up_X = 3^0.5,
                      L = NULL, low_F = 0.2, up_F = 0.9,
                      beta = NULL, gamma = NULL, low_Y = -0.9, up_Y = 0.9,
                      burn_in = 20, data_only = FALSE,
                      only_stationary = TRUE, epsilon = 0.0002, max_it_station = 500, adjust_diag = FALSE,
                      geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)



data <- data_test
#Information.criteria(data = data,n=n,p=p,q=q,k=k,t=t, est.method = 2)
Information.criteria(data = data,n=n,p=p,q=q,k=k,t=t, est.method = 1)

