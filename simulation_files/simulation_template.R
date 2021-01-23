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
rm(list = ls())
q_simulation <- c(3)            #vector of number of factors per data set
T_simulation <- seq(10,20, 11)    #vector of number of periods per data set
n_simulation <- seq(10,15, 5)    #vector of number of signals per data set

number_iterations <- 2        #number of observations per combination of (q, T, n)
type_sigma_U <- 'diagonal'

#-------------------Playground Katrin---------------------------------------------





#-------------------Playground Marc-----------------------------------------------





#-------------------Playground Manu----------------------------------------------






#-------------------Change Nothing From Here Onward----------------------------------------------------------------------------------------------------------------------


#-------------------Define Cluster-----------------------------------------------
cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl = cl)

#------------------Simulation Study----------------------------------------------
seed_index <- 1
simulated_data <- as.data.frame(c())
start_time <- Sys.time()
for (q in q_simulation) {# start for q
  for (T in T_simulation) {# start for T
    for (n in n_simulation) {# start for n
      save_iterations <- c()
      for (iterations in 1:number_iterations){ # start for iterations
        set.seed(seed_index)
        seed_index = seed_index + 1
        clusterExport(cl, list('n', 'q', 'T'))
        data <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -3^0.5, up_X = 3^0.5,
                                low_F = 0.2, up_F = 0.9, burn_in = 20, data_only = FALSE, only_stationary = TRUE,
                                adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

        X_true <- data$X

        #compute Information Criteria
        ic_pca <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=1, kmax=5, ml_parallel = TRUE, ml_maxit = 5)
        ic_ml <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=2, kmax=5, ml_parallel = TRUE, ml_maxit = 5)

        #get all optimal estimates
        pca_bai_q <-ic_pca$number_BaiNg #Bai & NG optimal q pca
        pca_bic_q <- ic_pca$number_BIC #Bayesian optimal q pca
        ml_bai_q <- ic_ml$number_BaiNg #Bai & NG optimal q maximum likelihood
        ml_bic_q <- ic_ml$number_BIC #Bayesian optimal q maximum likelihood

        pca_bai_right = pca_bai_q == q
        pca_bic_right = pca_bic_q == q
        ml_bai_right = ml_bai_q == q
        ml_bic_right = ml_bic_q == q

        true_f <- data$F
        true_lambda <- data$L[[1]]

        pca_bic_f <- ic_pca$F_BIC*n
        pca_bic_lambda <- ic_pca$Lambda_BIC
        pca_bai_f <- ic_pca$F_BaiNg*n
        pca_bai_lambda <- ic_pca$Lambda_BaiNg
        pca_trueq_f <- ic_pca$F_true*n
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
        mses <- c()
        for (index_trans in 1:q) {
          mse_pca1 <- mean((true_f_transformed[index_trans,] - pca_trueq_f_transformed[index_trans,])^2)
          mse_pca2 <- mean((true_f_transformed[index_trans,] + pca_trueq_f_transformed[index_trans,])^2)
          mse_pca <- min(mse_pca1, mse_pca2)

          mse_ml1 <- mean((true_f_transformed[index_trans,] - ml_trueq_f_transformed[index_trans,])^2)
          mse_ml2 <- mean((true_f_transformed[index_trans,] + ml_trueq_f_transformed[index_trans,])^2)
          mse_ml <- min(mse_ml1, mse_ml2)
          mses <- c(mses, 'mse_pca_f' =mse_pca, 'mse_ml_f' = mse_ml )
        }


        #compute explained variance for all estimates
        var_X <- var(t(X_true))
        var_X_diag <- diag(var_X)

        pca_bic_f_variance_explained <- mean(diag(pca_bic_lambda %*% var(t(pca_bic_f)) %*% t(pca_bic_lambda)) / var_X_diag)
        pca_bai_f_variance_explained <- mean(diag(pca_bai_lambda %*% var(t(pca_bai_f)) %*% t(pca_bai_lambda)) / var_X_diag)
        pca_trueq_f_variance_explained <- mean(diag(pca_trueq_lambda %*% var(t(pca_trueq_f)) %*% t(pca_trueq_lambda)) / var_X_diag)

        ml_bic_f_variance_explained <- mean(diag(ml_bic_lambda %*% var(t(ml_bic_f)) %*% t(ml_bic_lambda)) / var_X_diag)
        ml_bai_f_variance_explained <- mean(diag(ml_bai_lambda %*% var(t(ml_bai_f)) %*% t(ml_bai_lambda)) / var_X_diag)
        ml_trueq_f_variance_explained <- mean(diag(ml_trueq_lambda %*% var(t(ml_trueq_f)) %*% t(ml_trueq_lambda)) / var_X_diag)

        #variance_explained_true <- true_lambda %*% var(t(true_f)) %*% t(true_lambda)

        #compute X_hat
        pca_bic_f_X_hat <- pca_bic_lambda %*% pca_bic_f
        pca_bai_f_X_hat <- pca_bai_lambda %*% pca_bai_f
        pca_trueq_f_X_hat <- pca_trueq_lambda %*% pca_trueq_f

        ml_bic_f_X_hat <- ml_bic_lambda %*% ml_bic_f
        ml_bai_f_X_hat <- ml_bai_lambda %*% ml_bai_f
        ml_trueq_f_X_hat <- ml_trueq_lambda %*% ml_trueq_f

        #compute MSE of X and X_hat
        pca_bic_mse_X <- mean((X_true - pca_bic_f_X_hat)^2)
        pca_bai_mse_X <- mean((X_true - pca_bai_f_X_hat)^2)
        pca_trueq_mse_X <- mean((X_true - pca_trueq_f_X_hat)^2)

        ml_bic_mse_X <- mean((X_true - ml_bic_f_X_hat)^2)
        ml_bai_mse_X <- mean((X_true - ml_bai_f_X_hat)^2)
        ml_trueq_mse_X <- mean((X_true - ml_trueq_f_X_hat)^2)

        #create vector to save over iterations
        save_for_iterations <- c('pca_bai_right'=pca_bai_right,'pca_bic_right'=pca_bic_right,'ml_bai_right'=ml_bai_right,'ml_bic_right'=ml_bic_right,
                                 mses, 'VE_pca_bic' =pca_bic_f_variance_explained, 'VE_pca_bic' = pca_bai_f_variance_explained,
                                 'VE_pca_trueq' = pca_trueq_f_variance_explained, 'VE_ml_bic' =pca_bic_f_variance_explained, 'VE_ml_bai' = pca_bai_f_variance_explained,
                                 'VE_ml_trueq' = pca_trueq_f_variance_explained, 'mse_pca_bic_X_Xhat' = pca_bic_mse_X, 'mse_pca_bai_X_Xhat' = pca_bai_mse_X,
                                 'mse_pca_trueq_X_Xhat' = pca_trueq_mse_X, 'mse_ml_bic_X_Xhat' = ml_bic_mse_X, 'mse_ml_bai_X_Xhat' = ml_bai_mse_X,
                                 'mse_ml_trueq_X_Xhat' = ml_trueq_mse_X )

        save_iterations <- rbind(save_iterations, save_for_iterations)

        now_time <- Sys.time()
        time <- now_time - start_time
        print(paste('q=',q,'; n=',n, '; T=',T, '; iteration=',iterations, '; time since start = ',round(time, 4),units(time), sep=''))
      }# end for iterations

      #append vector to huge matrix (or average over number_iterations and save rest in list)
      save_for_simulated_data <- c(colMeans(save_iterations), 'n' = n, 'T'=T, 'q'=q, 'p'=0, 'k'=1, 'number of iterations per setting' = number_iterations,
                                    'maxit' = 5, 'max_factors_IC' = 6, 'DGP' = 'stationary', 'sigma_u' = type_sigma_U)
      colname <- names(save_for_simulated_data)
      simulated_data <- rbind(simulated_data, save_for_simulated_data)

    }# end for n
  }# end for T
}# end for q
colnames(simulated_data) <- colname

