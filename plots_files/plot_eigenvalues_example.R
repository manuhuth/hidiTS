library(ggplot2)
library("reshape2")
library(hidiTS)

n <- 20
T <- 20
q <-3 
lamb <- 3
set.seed(123)
data <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -lamb, up_X = lamb, 
                 low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = diag(n),
                 adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)


vec <- c()
vec_bic <- c()
for (k in 1:20) {
  pc <- pca_estimator(data$X, k)
  x_hat <- pc$Lambda %*% pc$F
  rss <- mean((data$X - x_hat)^2)
  #print(rss/n*t)
  bic <- n * log(rss/(n*T)) + k * log(n)
  bai <- log(rss/(n*T)) + k *(n+T)/(n*T) * log(min(n,T))
  vec <- cbind(vec, bai)
  vec_bic <- c(vec_bic, bic)
}
#plot(t(vec))

#plot(vec_bic)

#plot(eigen(var(t(data$X)))$values)

data_smooth <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -0.001, up_X = 0.001, 
                 low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = diag(n),
                 adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

#eigen(var(t(data_smooth$X)))$values

df <- data.frame(cbind(eigen(var(t(data$X)))$values,eigen(var(t(data_smooth$X)))$values, 1:length(eigen(var(t(data$X)))$values)))
colnames(df) <- c('Cut','Smooth', 'Eigenvalue')

df_long <- melt(df, id='Eigenvalue')
colnames(df_long) <- c('Eigenvalue', 'Type', 'Magnitude')

ggplot(df_long, aes(x=Eigenvalue, y=Magnitude, colour=Type)) + geom_point() +
  ggtitle(' ', subtitle = 'q = 3, t = 20, n = 20') 


ggplot(df, aes(x=Eigenvalue, y=Magnitude)) + geom_point() +
  ggtitle('Number of parameter to estimate', subtitle = 'q=3, t = 20') 