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
T_simulation <- seq(10,20, 11)    #vector of number of periods per data set
n_simulation <- seq(10,20, 11)    #vector of number of signals per data set

number_iterations <- 1        #number of observations per combination of (q, T, n)

#-------------------Playground Katrin---------------------------------------------





#-------------------Playground Marc-----------------------------------------------





#-------------------Playground Manu----------------------------------------------






#-------------------Change Nothing From Here Onward----------------------------------------------------------------------------------------------------------------------


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
        ic_ml <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=2, kmax=6, ml_parallel = TRUE, ml_maxit = 5)

        #get all optimal estimates
        pca_bai_q <-ic_pca$number_BaiNg #Bai & NG optimal q pca
        pca_bic_q <- ic_pca$number_BIC #Bayesian optimal q pca
        ml_bai_q <- ic_ml$number_BaiNg #Bai & NG optimal q maximum likelihood
        ml_bic_q <- ic_ml$number_BIC #Bayesian optimal q maximum likelihood

        true_f <- data$F
        true_lambda <- data$L[[1]]

        pca_bic_f <- ic_pca$F_BIC
        pca_bic_lambda <- ic_pca$Lambda_BIC
        pca_bai_f <- ic_pca$F_BaiNg
        pca_bai_lambda <- ic_pca$Lambda_BaiNg
        pca_trueq_f <- ic_pca$F_true
        pca_trueq_lambda <- ic_pca$Lambda_true

        ml_bic_f <- ic_ml$F_BIC
        ml_bic_lambda <- ic_ml$Lambda_BIC
        ml_bai_f <- ic_ml$F_BaiNg
        ml_bai_lambda <- ic_ml$Lambda_BaiNg
        ml_trueq_f <- ic_ml$F_true
        ml_trueq_lambda <- ic_ml$Lambda_true

        #transform all estimates
        true_f_transformed <- convert_factors_dgp(data)$F
        true_lambda_transformed <- convert_factors_dgp(data)$Lambda


        #pca_bic_f_transformed <- convert_factors_ML(Lambda=pca_bic_lambda, factors=pca_bic_f, q=pca_bic_q)$F
        #pca_bic_lambda_transformed <- convert_factors_ML(Lambda=pca_bic_lambda, factors=pca_bic_f, q=pca_bic_q)$Lambda
        #pca_bai_f_transformed <- convert_factors_ML(Lambda=pca_bai_lambda, factors=pca_bai_f, q=pca_bai_q)$F
        #pca_bai_lambda_transformed <- convert_factors_ML(Lambda=pca_bai_lambda, factors=pca_bai_f, q=pca_bai_q)$Lambda
        pca_trueq_f_transformed <- convert_factors_ML(Lambda=pca_trueq_lambda, factors=pca_trueq_f, q=q)$F
        pca_trueq_lambda_transformed <- convert_factors_ML(Lambda=pca_trueq_lambda, factors=pca_trueq_f, q=q)$Lambda

        #ml_bic_f_transformed <- convert_factors_ML(Lambda=ml_bic_lambda, factors=ml_bic_f, q=ml_bic_q)$F
        #ml_bic_lambda_transformed <- convert_factors_ML(Lambda=ml_bic_lambda, factors=ml_bic_f, q=ml_bic_q)$Lambda
        #ml_bai_f_transformed <- convert_factors_ML(Lambda=ml_bai_lambda, factors=ml_bai_f, q=ml_bai_q)$F
        #ml_bai_lambda_transformed <- convert_factors_ML(Lambda=ml_bai_lambda, factors=ml_bai_f, q=ml_bai_q)$Lambda
        ml_trueq_f_transformed <- convert_factors_ML(Lambda=ml_trueq_lambda, factors=ml_trueq_f, q=q)$F
        ml_trueq_lambda_transformed <- convert_factors_ML(Lambda=ml_trueq_lambda, factors=ml_trueq_f, q=q)$Lambda

        #compute MSE of f and f_hat for (true?) all estimates
        for (index_trans in 1:q) {
          mse_pca1 <- mean((true_f_transformed[index_trans,] - pca_trueq_f_transformed[index_trans,])^2)
          mse_pca2 <- mean((true_f_transformed[index_trans,] + pca_trueq_f_transformed[index_trans,])^2)

          mse_ml1 <- mean((true_f_transformed[index_trans,] - ml_trueq_f_transformed[index_trans,])^2)
          mse_ml2 <- mean((true_f_transformed[index_trans,] + ml_trueq_f_transformed[index_trans,])^2)

          if (mse_pca1 > mse_pca2 ){pca_trueq_f_transformed[index_trans,] <- - pca_trueq_f_transformed[index_trans,]}
          if (mse_ml1 > mse_ml2 ){ml_trueq_f_transformed[index_trans,] <- - ml_trueq_f_transformed[index_trans,]}
        }
        #compute explained variance for all estimates

        #compute X_hat

        #compute MSE of X and X_hat

        #create vector to append

        #append vector to huge matrix (or average over number_iterations and save rest in list)


      }# end for iterations
    }# end for n
  }# end for T
}# end for q

