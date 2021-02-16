library(hidiTS)
library(ggplot2)

###write ics

set.seed(1)
n <- 211
T <- 104
q <- 3


BIC_n_df <- c()
BIC_T_df <- c()
BIC_nT_df <- c()
BNIC_df <- c()

for (index in 1:1000) {
  
  data <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -3^0.5, up_X = 3^0.5,
                   low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = diag(n),
                   adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)
  
  X <- data$X
  
  
  ic_BIC_n  <- c()
  ic_BIC_T  <- c()
  ic_BIC_nT  <- c()
  ic_BNIC  <- c()
  
  for (r in 1:n){
    pca_est <- pca_estimator(X, r)
    f_hat <- pca_est$F
    l_hat <- pca_est$Lambda
    x_hat <- l_hat %*% f_hat
    u_hat <- X - x_hat
    RSS <- sum(u_hat^2)
    
    RSS_part <- log(RSS/n/T)
    
    BIC_n <- RSS_part + r * log(n)/n
    BIC_T <- RSS_part + r * log(T)/T
    BIC_nT <- RSS_part + r * log(n*T)/n/T * (n+T-r)
    BNIC <- RSS_part + r *  r * log(min(n,T))/n/T * (n+T)
    
    ic_BIC_n[r]  <- BIC_n
    ic_BIC_T[r]  <-  BIC_T
    ic_BIC_nT[r]  <-  BIC_nT
    ic_BNIC[r]  <- BNIC
    print(paste('sample=', index, ' r=',r, sep=''))
  }
  
  
  BIC_n_df <- rbind(BIC_n_df, ic_BIC_n)
  BIC_T_df <- rbind(BIC_T_df, ic_BIC_T)
  BIC_nT_df <- rbind(BIC_nT_df, ic_BIC_nT)
  BNIC_df <- rbind(BNIC_df, ic_BNIC)
}


alpha = 0.1
up <- 1 - alpha/2
low <-  alpha/2




mean_n <- colMeans(BIC_n_df)
high_n <- apply(BIC_n_df, 2, quantile, up)
low_n <- apply(BIC_n_df, 2, quantile, low)
BIC_n_df_plot <- cbind(low_n, mean_n, high_n, 'r'= 1:n)[1:(n-20),]

mean_T <- colMeans(BIC_T_df)
high_T <- apply(BIC_T_df, 2, quantile, up)
low_T <- apply(BIC_T_df, 2, quantile, low)
BIC_T_df_plot <- cbind(low_T, mean_T, high_T,  'r'= 1:n)[1:(n-20),]


mean_nT <- colMeans(BIC_nT_df)
high_nT <- apply(BIC_nT_df, 2, quantile, up)
low_nT <- apply(BIC_nT_df, 2, quantile, low)
BIC_nT_df_plot <- cbind(low_nT, mean_nT, high_nT,  'r'= 1:n)[1:(n-20),]

mean_bai <- colMeans(BNIC_df)
high_bai <- apply(BNIC_df, 2, quantile, up)
low_bai <- apply(BNIC_df, 2, quantile, low)
BNIC_df_plot <- cbind(low_bai, mean_bai, high_bai,  'r'= 1:n)[1:15,]



ggplot(as.data.frame(BIC_n_df_plot) , aes(x = r , y = mean_n)) +
  geom_ribbon(aes(ymin = low_n, ymax = high_n), fill = "steelblue2") +
  ylab("BIC_n") + geom_point(aes(x=r, y = mean_n), size=1)

ggplot(as.data.frame(BNIC_df_plot) , aes(x = r , y = mean_bai)) +
  geom_ribbon(
    aes(ymin = low_bai, ymax = high_bai), fill = "steelblue2")+
  ylab("BNIC") + geom_point(aes(x=r, y = mean_bai), size=2, shape=1)

ggplot(as.data.frame(BIC_T_df_plot) , aes(x = r , y = mean_T)) +
  geom_ribbon(
    aes(ymin = low_T, ymax = high_T), fill = "steelblue2")+
  ylab("BIC_T") + geom_point(aes(x=r, y = mean_T), size=1)

ggplot(as.data.frame(BIC_nT_df_plot) , aes(x = r , y = mean_nT)) +
  geom_ribbon(aes(ymin = low_nT, ymax = high_nT), fill = "steelblue2")+
  ylab("BIC_nT") + geom_point(aes(x=r, y = mean_nT), size=2, shape=1)
