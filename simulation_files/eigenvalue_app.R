#-------------------Simulation Template--------------------------------------------
#Katrin Marques Magalhaes
#Marc Kerstan
#Manuel Huth


#-------------------Load Libraries-------------------------------------------------
#devtools::install(dependencies = TRUE)
#devtools::install_github("manuhuth/hidiTS")
#install.packages("xtable")
install.packages("a4Reporting")
library(a4Reporting)
library(xtable)
library(hidiTS)
library(fbi)
library(optimParallel)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
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

IC_rslt <- Information.criteria.app(data=df_final,n=n,p=0,q=2,k=1,t=t,est.method = 1,kmax = 10)

#setimate hidimodel
pca_IC <- pca_estimator(df_final$X,IC_rslt$number_BaiNg) 

est_lambda <- pca_IC$Lambda
est_f <- pca_IC$F

X_est <- est_lambda %*% est_f

res <- df_final$X - X_est

#get sigma_u and correlation matrix 
est_sigma_u=var(t(res))

cor_matrix <- cor(t(res))

#calculate eigenvalues
eigenvalues_data <- eigen(var(t(df_final$X)))$values

df_eig <- as.data.frame(eigenvalues_data)

# to do: create table 

lambda_vec <- as.vector(est_lambda)
lambda_dens <- plot(density(lambda_vec)) 


t(df_eig)


xtable(df_eig[1:10, ],auto=TRUE)

print(xtable(t(df_eig), type = "latex", digits=c(2) ), include.rownames=FALSE)


