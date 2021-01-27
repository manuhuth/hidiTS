# executed with modified estimate_f and starting_values_ML
Information.criteria <- function(data,n,p,q,k,t, est.method, kmax=8, ml_parallel = TRUE, ml_maxit = 5){

  data_x <- data$X
  if(est.method==1){
    p<-0
    r<-(p+1)*q
    g<-(n+t)/(n*t) *log(min(n,t))
    PC.IC <- function(num){
      pca_est <- pca_estimator(X=data$X,number = num)
      rss<-sum(diag(t(data$X-(pca_est$Lambda%*%pca_est$F))%*%(data$X-(pca_est$Lambda%*%pca_est$F)))/(n*T))
      bayesian <- log(n)*num + n*log(rss)
      baing <- g*num + log(rss)
      return(list(bayesian, baing, pca_est$F, pca_est$Lambda ))

    }

    output1<- sapply(1:kmax, PC.IC)
    BIC.PC<-unlist(output1[1,])
    BaiNg.PC<- unlist(output1[2,])

    num.bic<- match(min(BIC.PC),BIC.PC)
    num.BaiNg<-match(min(BaiNg.PC),BaiNg.PC)

    F_hat_BIC<- output1[3,num.bic][[1]]
    Lambda_hat_BIC<-output1[4,num.bic][[1]]

    F_hat_BaiNg<-output1[3,num.BaiNg][[1]]
    Lambda_hat_BaiNg<- output1[4,num.BaiNg][[1]]

    true_q<- nrow(data$F)
    F_true<- output1[3,true_q][[1]]
    Lambda_true<- output1[4,true_q][[1]]

    return(list("true_q"=true_q, "number_BIC"=num.bic, "number_BaiNg"=num.BaiNg,"F_true"=F_true,
                "Lambda_true"=Lambda_true, "F_BIC"=F_hat_BIC, "Lambda_BIC"=Lambda_hat_BIC,
                      "F_BaiNg"=F_hat_BaiNg, "Lambda_BaiNg"=Lambda_hat_BaiNg))

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
