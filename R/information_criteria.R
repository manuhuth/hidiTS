# executed with modified estimate_f and starting_values_ML
Information.criteria <- function(data,n,p,q,k,t, est.method, kmax=8){

  data_x <- data$X
  if(est.method==1){
    r<-(p+1)*q
    g<-(n+t)/(n*t) *log(min(n,t))
    PC.IC <- function(num){
      pca_est <- pca_estimator(X=data$X,number = num)
      rss<-sum(diag(t(data$X-(pca_est$Lambda%*%pca_est$F*n))%*%(data$X-(pca_est$Lambda%*%pca_est$F*n)))/(n*T))
      bayesian <- log(n)*num + n*log(rss)
      baing <- g*num + log(rss)
      return(data.frame(bayesian, baing))

    }

    output1<- sapply(1:kmax, PC.IC)
    BIC.PC<-unlist(output1[1,])
    BaiNg.PC<- unlist(output1[2,])
    num.bic<- match(min(BIC.PC),BIC.PC)
    num.BaiNg<-match(min(BaiNg.PC),BaiNg.PC)

    return(data.frame("true_r"=r, "number_BIC"=num.bic, "number_BaiNg"=num.BaiNg))

  }else if(est.method==2){

    g<-(n+t)/(n*t) *log(min(n,t))
     calc_ll<- function(q,p){
      
      if(((p+1)*q) <= n ){
         est <- estimate_f(data_x=data_test,n=n,p=p,q=q,k=k,t=t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=TRUE,it=1,
                          method = "L-BFGS-B", parallel = FALSE, max_it = 3)
        
        return(est$value)
        
      }else{
        return(NaN)
      }
    }
    qvalues<- 1:kmax
    pvalues<-0:kmax
    exp_grid<- expand.grid(qvalues, pvalues)
    temp<- do.call(Vectorize(calc_ll), unname(exp_grid))
    ll<- matrix(temp, ncol = kmax, byrow = TRUE )
    ll[is.nan(ll)]<-0

    ll<- -ll
    bayesian<-c()
    Bai.Ng<-c()
    bic.lag<-c()
    Bai.Ng.lag<-c()

    ML.IC <-function(num){
      bayesian[num] <- log(n)*num - 2*log(max(ll[,num]))
      Bai.Ng[num]<- g*num -2*log(max(ll[,num]))
      return(data.frame(bayesian[num], Bai.Ng[num]))
    }
    ML.IC.Lags<- function(num_lag){
      bic.lag[num_lag]<- log(n)*num_lag - 2*log(max(ll[num_lag,]))
      Bai.Ng.lag[num_lag]<- g*num_lag -2*log(max(ll[,num_lag]))
      return(data.frame(bic.lag[num_lag], Bai.Ng.lag[num_lag]))
    }
    output<-sapply(1:kmax, ML.IC)
    BIC.ML<-unlist(output[1,])
    BaiNg.ML<- unlist(output[2,])
    num.bic<-  match(min(BIC.ML),BIC.ML)
    num.baing<- match(min(BaiNg.ML),BaiNg.ML)

    output.lag<- sapply(1:kmax, ML.IC.Lags)
    BIC.ML.lag<-unlist(output.lag[1,])
    BaiNg.ML.lag<- unlist(output.lag[2,])
    lags.bic<- match(min(BIC.ML.lag),BIC.ML.lag)
    lags.baing<- match(min(BaiNg.ML.lag),BaiNg.ML.lag)

    return(data.frame("true_number"=q, "true_lag"=p, "number_BIC"=num.bic,"number_BaiNg" =num.baing, "lags_BIC"=lags.bic, "lags_BaiNg"=lags.baing))
  }


  else{
    print("Choose est.method either 1 or 2")
  }
}
