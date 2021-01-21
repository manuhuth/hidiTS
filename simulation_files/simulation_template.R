#-------------------Simulation Template--------------------------------------------
#Katrin Marques Magalhaes
#Marc Kerstan
#Manuel Huth


#-------------------Load Libraries-------------------------------------------------
#devtools::install(dependencies = TRUE)
#devtools::install_github("manuhuth/hidiTS")
library(hidiTS)
library(optimParallel)

#-------------------Specify Iterations and other Parameters------------------------
q_simulation <- c(3)            #vector of number of factors per data set
T_simulation <- seq(10,20, 5)    #vector of number of periods per data set
n_simulation <- seq(10,20, 5)    #vector of number of signals per data set

number_iterations <- 100        #number of observations per combination of (q, T, n)

#-------------------Playground Katrin---------------------------------------------





#-------------------Playground Marc-----------------------------------------------





#-------------------Playground Manu----------------------------------------------





#-------------------Define Cluster-----------------------------------------------
cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl = cl)

#------------------Simulation Study----------------------------------------------
seed_index <- 1
for (q in q_simulation) {# start for q
  for (T in T_simulation) {# start for T
    for (n in n_simulation) {# start for n
      for (iterations in 1:number_iterations){ # start for iterations
        set.seed(seed_index)
        seed_index = seed_index + 1
        clusterExport(cl, list('n', 'q', 'T'))
        data <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -3^0.5, up_X = 3^0.5,
                                low_F = 0.2, up_F = 0.9, burn_in = 20, data_only = FALSE, only_stationary = TRUE,
                                adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

        #compute Information Criteria
        ic_pca <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=1, kmax=6, ml_parallel = TRUE, ml_maxit = 5)
        pca_optimal_bai <- ic_pca['number_BaiNg']
        pca_optimal_bic <- ic_pca['number_BIC']
        ml_optimal <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=2, kmax=6, ml_parallel = TRUE, ml_maxit = 5)

        #get all optimal estimates

        #transform all estimates

        #compute MSE of f and f_hat for (true?) all estimates

        #compute explained variance for all estimates

        #compute X_hat

        #compute MSE of X and X_hat

        #create vector to append

        #append vector to huge matrix (or average over number_iterations and save rest in list)


      }# end for iterations
    }# end for n
  }# end for T
}# end for q

