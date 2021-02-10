#-------------------Simulation Template--------------------------------------------
#Katrin Marques Magalhaes
#Marc Kerstan
#Manuel Huth


#-------------------Load Libraries-------------------------------------------------
#devtools::install(dependencies = TRUE)
#devtools::install_github("manuhuth/hidiTS")
library(hidiTS)
library(fbi)
library(optimParallel)

rm(list = ls())

#-------------------Load and prepare Data -----------------------------------------

prep_data_hidi <- function(data){
  
  data_hidi <- as.matrix(t(data))
  
  
  output <- list("data_fbi"= data, "X"=data_hidi )
  
  return(output)
}


Information.criteria.app <- function(data,n,p,q,k,t, est.method, kmax=8, ml_parallel = TRUE, ml_maxit = 5){
  
  data_x <- data$X
  if(est.method==1){
    p<-0
    r<-(p+1)*q
    g<-(n+t)/(n*t) *log(min(n,t))
    PC.IC <- function(num){
      
      g2<- (n+t-num) * log(n*t) / (n*t)
      
      pca_est <- pca_estimator(X=data$X,number = num)
      #pca_est_max <- pca_estimator(X=data$X,number = kmax)
      rss<-sum(diag(t(data$X-(pca_est$Lambda%*%pca_est$F))%*%(data$X-(pca_est$Lambda%*%pca_est$F)))/(n*T))
      #rss_max<-sum(diag(t(data$X-(pca_est_max$Lambda%*%pca_est_max$F))%*%(data$X-(pca_est_max$Lambda%*%pca_est_max$F)))/(n*T))
      bayesian <- log(n)/n * num  + log(rss)
      baing <- g*num + log(rss)
      
      bayesian_t <- log(t)/t * num  + log(rss)
      bayesian_nt <- log(rss) + num*(g2)
      return(list(bayesian, baing, pca_est$F, pca_est$Lambda , bayesian_t, bayesian_nt))
      
    }
    
    output1<- sapply(1:kmax, PC.IC)
    BIC.PC<-unlist(output1[1,])
    BaiNg.PC<- unlist(output1[2,])
    BIC_t.PC<- unlist(output1[5,]) 
    BIC_nt.PC <- unlist(output1[6,])
    
    num.bic<- match(min(BIC.PC),BIC.PC)
    num.BaiNg<-match(min(BaiNg.PC),BaiNg.PC)
    num.bic_t <- match(min(BIC_t.PC), BIC_t.PC) 
    num.bic_nt <- match(min(BIC_nt.PC),BIC_nt.PC)
    
    F_hat_BIC<- output1[3,num.bic][[1]]
    Lambda_hat_BIC<-output1[4,num.bic][[1]]
    
    F_hat_BaiNg<-output1[3,num.BaiNg][[1]]
    Lambda_hat_BaiNg<- output1[4,num.BaiNg][[1]]
    
    
    F_hat_BIC_t<-output1[3,num.bic_t][[1]]      
    Lambda_hat_BIC_t<- output1[4,num.bic_t][[1]]
    
    F_hat_BIC_nt<-output1[3,num.bic_nt][[1]]      
    Lambda_hat_BIC_nt<- output1[4,num.bic_nt][[1]]
    
    #true_q<- nrow(data$F)
    #F_true<- output1[3,true_q][[1]]
    #Lambda_true<- output1[4,true_q][[1]]
    
    return(list("number_BIC"=num.bic, "number_BaiNg"=num.BaiNg, "number_BIC_t"=num.bic_t, 
                "number_BIC_nt"=num.bic_nt,"F_BIC"=F_hat_BIC, "Lambda_BIC"=Lambda_hat_BIC,
                "F_BaiNg"=F_hat_BaiNg, "Lambda_BaiNg"=Lambda_hat_BaiNg, "F_BIC_t"=F_hat_BIC_t,
                "Lambda_BIC_t"= Lambda_hat_BIC_t, "F_BIC_nt"=F_hat_BIC_nt, "Lambda_BIC_nt"=Lambda_hat_BIC_nt))
    
  }else if(est.method==2){
    
    g<-(n+t)/(n*t) *log(min(n,t))
    calc_ll<- function(q,p){
      
      if(((p+1)*q) <= n ){
        if (isTRUE(ml_parallel)) {
          clusterExport(cl, list('n', 'p', 'q', 't', 'k'), envir=environment())
        }
        est <- estimate_f(data_x=data,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
                          method = "L-BFGS-B", parallel = ml_parallel, max_it = ml_maxit)
        
        return(list(est$value, est$F_final, est$lambda))
        
      }else{
        return(NaN)
      }
    }
    qvalues<- 1:kmax
    pvalues<-0
    exp_grid<- expand.grid(qvalues, pvalues)
    temp<- do.call(Vectorize(calc_ll), unname(exp_grid))
    
    save_ll<-temp[1,]
    ll<- matrix(unlist(save_ll), ncol=kmax,nrow=1, byrow = TRUE )
    ll[is.nan(ll)]<-0
    
    ll<- -ll
    bayesian<-c()
    Bai.Ng<-c()
    bic.lag<-c()
    Bai.Ng.lag<-c()
    
    ML.IC <-function(num){
      bayesian[num] <- log(T)*number_of_param(n=n,p=p,q=q,k=num,gamma_res=TRUE, lambda_res=TRUE, sigma_u_diag=TRUE) - 2*max(ll[,num])
      Bai.Ng[num]<- g*number_of_param(n=n,p=p,q=q,k=num,gamma_res=TRUE, lambda_res=TRUE, sigma_u_diag=TRUE) -2*max(ll[,num])
      return(data.frame(bayesian[num], Bai.Ng[num]))
    }
    #ML.IC.Lags<- function(num_lag){
    #  bic.lag[num_lag+1]<- log(n)*num_lag - 2*log(max(ll[num_lag+1,]))
    #  Bai.Ng.lag[num_lag+1]<- g*num_lag -2*log(max(ll[num_lag+1,]))
    #  return(data.frame(bic.lag[num_lag+1], Bai.Ng.lag[num_lag+1]))
    #}
    output<-sapply(1:kmax, ML.IC)
    BIC.ML<-unlist(output[1,])
    BaiNg.ML<- unlist(output[2,])
    num.bic<-  match(min(BIC.ML),BIC.ML)
    num.BaiNg<- match(min(BaiNg.ML),BaiNg.ML)
    
    #output.lag<- sapply(0:kmax, ML.IC.Lags)
    #BIC.ML.lag<-unlist(output.lag[1,])
    #BaiNg.ML.lag<- unlist(output.lag[2,])
    #lags.bic<- match(min(BIC.ML.lag),BIC.ML.lag) -1
    #lags.baing<- match(min(BaiNg.ML.lag),BaiNg.ML.lag)-1
    
    save_F<-temp[2,]
    save_lambda<- temp[3,]
    #F_hat_BIC<- save_F[[kmax*lags.bic + num.bic]] # if p =!0
    #Lambda_hat_BaiNg<- save_lambda[[kmax*lags.baing+num.BaiNg]]
    
    true_q<- nrow(data$F)
    F_true<- save_F[[true_q]]
    Lambda_true<- save_lambda[[true_q]][[1]]
    
    F_hat_BIC<- save_F[[num.bic]]
    F_hat_BaiNg<- save_F[[num.BaiNg]]
    Lambda_hat_BIC<- save_lambda[[num.bic]][[1]]
    Lambda_hat_BaiNg<- save_lambda[[num.BaiNg]][[1]]
    
    return(list("true_number"=true_q, "number_BIC"=num.bic,"number_BaiNg" =num.BaiNg,"F_true"=F_true,
                "Lambda_true"=Lambda_true,"F_BIC"=F_hat_BIC, "Lambda_BIC"=Lambda_hat_BIC,
                "F_BaiNg"=F_hat_BaiNg, "Lambda_BaiNg"=Lambda_hat_BaiNg))
  }
  
  
  else{
    print("Choose est.method either 1 or 2")
  }
}


df_paper=fredqd(date_end =as.Date(c("2019/12/01")), date_start=as.Date(c("1960/03/01")),transform = TRUE)

# remove dates
df_temp=df_paper[,2:(length(df_paper))]

#drop Nas
Na_names <- colnames(df_temp)[colSums(is.na(df_temp)) > 0]
df_colnames <- colnames(df_temp)
col_keeps <- setdiff(df_colnames,Na_names)
df_no_na <- subset(df_temp, select=(col_keeps))

df_final <- prep_data_hidi(df_no_na)
# df_final$X txt
# Adjust data to DGP

n=nrow(df_final$X)
t=ncol(df_final$X)

IC_rslt <- Information.criteria.app(data=df_final,n=n,p=0,q=2,k=1,t=t,est.method = 1)

#setimate hidimodel
pca_IC <- pca_estimator(df_final$X,IC_rslt$number_BaiNg) 

est_lambda <- pca_IC$Lambda
est_f <- pca_IC$F

X_est <- est_lambda %*% est_f

res <- df_final$X - X_est

#get sigma_u and correlation matrix 
est_sigma_u=var(t(res))

cor_matrix <- cor(t(res))


#-------------------Specify Iterations and other Parameters------------------------

econometrician <- 'Marc' #either 'Katrin', 'Marc' or 'Manu'

q_simulation <- seq(2,10,1)            #vector of number of factors per data set
T_simulation <- c(t)  #vector of number of periods per data set
n_simulation <- c(n)
number_iterations <- 2       #number of observations per combination of (q, T, n)
type_sigma_U <- 'Application'


#-------------------Playground Katrin---------------------------------------------
if (econometrician == 'Katrin'){
  print('Roses are red, model the shock, this code is gonna rock!')
  
  
  lamb <- 3^0.5
  
  
  
  
}
#-------------------Playground Marc-----------------------------------------------
if (econometrician == 'Marc'){
  print('May the kernels and bug fixes be with you! #Bugsafari')
  
  
  
  lamb <- 3
  
  
  
}
#-------------------Playground Manu----------------------------------------------
if (econometrician == 'Manu'){
  
  
  lamb <- 0.5
  
  
  
}
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
      
      
      Gamma_sim <- diag(-sort(-runif(q, 0.2, 0.8)))
      Lambda_sim <- matrix(runif(n*q, -lamb, lamb), n, q)
      
      for (iterations in 1:number_iterations){ # start for iterations
        tryCatch({
          set.seed(seed_index)
          seed_index = seed_index + 1
          
          clusterExport(cl, list('n', 'q', 'T'))
          data <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -lamb, up_X = lamb, A = list(Gamma_sim), L = list(Lambda_sim),
                           low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = cor_matrix,
                           adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)
          
          X_true <- data$X
          eigenvalues_x <- eigen(var(t(data$X)))$values
          
          
          #compute Information Criteria
          ic_pca <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=1, kmax=10, ml_parallel = TRUE, ml_maxit = 5)
          start_time_pca <- Sys.time()
          placeholder <- pca_estimator(X_true, q)
          end_time_pca <- Sys.time()
          time_pca <- end_time_pca - start_time_pca
          #ic_ml <- Information.criteria(data=data,n=n,p=0,q=q,k=1,t=T, est.method=2, kmax=2, ml_parallel = TRUE, ml_maxit = 5)
          
          #get all optimal estimates
          pca_bai_q <-ic_pca$number_BaiNg #Bai & NG optimal q pca
          pca_bic_q <- ic_pca$number_BIC #Bayesian optimal q pca
          ml_bai_q <- ic_pca$number_BaiNg #Bai & NG optimal q maximum likelihood
          ml_bic_q <- ic_pca$number_BIC #Bayesian optimal q maximum likelihood
          
          pca_bai_right = pca_bai_q == q
          pca_bic_right = pca_bic_q == q
          ml_bai_right = ml_bai_q == q
          ml_bic_right = ml_bic_q == q
          bai_equals_bic <- pca_bai_q == pca_bic_q
          
          true_f <- data$F
          true_lambda <- data$L[[1]]
          
          pca_bic_f <- ic_pca$F_BIC
          pca_bic_lambda <- ic_pca$Lambda_BIC
          pca_bai_f <- ic_pca$F_BaiNg
          pca_bai_lambda <- ic_pca$Lambda_BaiNg
          pca_trueq_f <- ic_pca$F_true
          pca_trueq_lambda <- ic_pca$Lambda_true
          
          
          #transform all estimates
          true_f_transformed <- convert_factors_dgp(data)$F
          true_lambda_transformed <- convert_factors_dgp(data)$Lambda
          
          
          #pca_bic_f_transformed <- convert_factors_ML(Lambda=pca_bic_lambda, factors=pca_bic_f, q=pca_bic_q)$F
          #pca_bic_lambda_transformed <- convert_factors_ML(Lambda=pca_bic_lambda, factors=pca_bic_f, q=pca_bic_q)$Lambda
          #pca_bai_f_transformed <- convert_factors_ML(Lambda=pca_bai_lambda, factors=pca_bai_f, q=pca_bai_q)$F
          #pca_bai_lambda_transformed <- convert_factors_ML(Lambda=pca_bai_lambda, factors=pca_bai_f, q=pca_bai_q)$Lambda
          pca_trueq_f_transformed <- convert_factors_ML(Lambda=pca_trueq_lambda, factors=pca_trueq_f, q=q)$F
          pca_trueq_lambda_transformed <- convert_factors_ML(Lambda=pca_trueq_lambda, factors=pca_trueq_f, q=q)$Lambda
          
          
          #compute MSE of f and f_hat for (true?) all estimates
          mses <- c()
          for (index_trans in 1:max(q_simulation)) {
            if(index_trans <= q){
              mse_pca1 <- mean((true_f_transformed[index_trans,] - pca_trueq_f_transformed[index_trans,])^2)
              mse_pca2 <- mean((true_f_transformed[index_trans,] + pca_trueq_f_transformed[index_trans,])^2)
              mse_pca <- min(mse_pca1, mse_pca2)
            }else{mse_pca=NaN}
            mses <- c(mses, 'mse_pca_f' =mse_pca)
          }
          
          
          #compute explained variance for all estimates
          var_X <- var(t(X_true))
          var_X_diag <- diag(var_X)
          
          pca_bic_f_variance_explained <- mean(diag(pca_bic_lambda %*% var(t(matrix(pca_bic_f,pca_bic_q ,T))) %*% t(pca_bic_lambda)) / var_X_diag)
          pca_bai_f_variance_explained <- mean(diag(pca_bai_lambda %*% var(t(matrix(pca_bai_f,pca_bai_q ,T))) %*% t(pca_bai_lambda)) / var_X_diag)
          pca_trueq_f_variance_explained <- mean(diag(pca_trueq_lambda %*% var(t(matrix(pca_trueq_f, q ,T))) %*% t(pca_trueq_lambda)) / var_X_diag)
          
          
          #variance_explained_true <- true_lambda %*% var(t(true_f)) %*% t(true_lambda)
          
          #compute X_hat
          pca_bic_f_X_hat <- pca_bic_lambda %*% pca_bic_f
          pca_bai_f_X_hat <- pca_bai_lambda %*% pca_bai_f
          pca_trueq_f_X_hat <- pca_trueq_lambda %*% pca_trueq_f
          
          
          #compute MSE of X and X_hat
          pca_bic_mse_X <- mean((X_true - pca_bic_f_X_hat)^2)
          pca_bai_mse_X <- mean((X_true - pca_bai_f_X_hat)^2)
          pca_trueq_mse_X <- mean((X_true - pca_trueq_f_X_hat)^2)
          
          
          #create vector to save over iterations
          ev_named <- eigenvalues_x[1:5]
          names(ev_named) <- c('ev1', 'ev2', 'ev3','ev4', 'ev5')
          
          #create vector to save over iterations
          save_for_iterations <- c(ev_named,'pca_bai_right'=pca_bai_right,'pca_bic_right'=pca_bic_right,
                                   'bic_equal_bai' = bai_equals_bic, mses, 'VE_pca_bic' =pca_bic_f_variance_explained, 'VE_pca_bai' = pca_bai_f_variance_explained,
                                   'VE_pca_trueq' = pca_trueq_f_variance_explained, 'mse_pca_bic_X_Xhat' = pca_bic_mse_X, 'mse_pca_bai_X_Xhat' = pca_bai_mse_X,
                                   'mse_pca_trueq_X_Xhat' = pca_trueq_mse_X, 'time_pca'=time_pca)
          
          save_iterations <- rbind(save_iterations, save_for_iterations)
          
          now_time <- Sys.time()
          time <- now_time - start_time
          print(paste('q=',q,'; n=',n, '; T=',T, '; iteration=',iterations, '; time since start = ',round(time, 4),units(time), sep=''))
          
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n", 'seed=', seed_index ,'q=',q,'; n=',n, '; T=',T)})
        
        
      }# end for iterations
      
      #append vector to huge matrix (or average over number_iterations and save rest in list)
      save_for_simulated_data <- data.frame(t(colMeans(save_iterations)), 'n' = n, 'T'=T, 'q'=q, 'p'=0, 'k'=1, 'number of iterations per setting' = number_iterations,
                                            'maxit' = 5, 'max_factors_IC' = 6, 'DGP' = 'stationary', 'sigma_u' = type_sigma_U)
      colname <- names(save_for_simulated_data)
      simulated_data <- rbind(simulated_data, save_for_simulated_data)
      
    }# end for n
  }# end for T
}# end for q
colnames(simulated_data) <- colname
