check_speed <- function(parallel = TRUE, seed = 123, iterations_loop = 50, max_it_optim = 20, n=10, p=0, t=20, q=2, k=1){
    library("optimParallel")
    
    
    if (isTRUE(parallel)) {
      cl <- makeCluster(detectCores()-1)
      setDefaultCluster(cl = cl)
      clusterExport(cl, list('n', 'p', 'q', 't', 'k'), envir=environment())
    }
    
    set.seed(seed)
    start <- Sys.time()
    for (index in 1:50) {
      data_test <- sim_data(p = n, T = t, dim_F= q, lags_F=k, lags_X=p, ar_F=1, ar_Y=1, A = NULL, low_X = -3^0.5, up_X = 3^0.5,
                            L = NULL, low_F = 0.2, up_F = 0.9,
                            beta = NULL, gamma = NULL, low_Y = -0.9, up_Y = 0.9,
                            burn_in = 20, data_only = FALSE,
                            only_stationary = TRUE, epsilon = 0.0002, max_it_station = 500, adjust_diag = FALSE,
                            geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)
      
      if (isTRUE(parallel)) {
        est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
                          method = "L-BFGS-B", parallel = TRUE, max_it = max_it_optim, trace=0, forward = TRUE, loginfo=FALSE)
      } else{
        est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
                          method = "L-BFGS-B", parallel = FALSE, max_it = max_it_optim, trace=0)
      }
    
      print(index)
    }
    end <- Sys.time()
    time <- end - start
    
    return(time)
}
