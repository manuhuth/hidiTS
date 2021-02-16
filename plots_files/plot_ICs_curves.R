#devtools::install(dependencies = TRUE)
#devtools::install_github("manuhuth/hidiTS")
#install.packages("xtable")
install.packages("a4Reporting")
library(a4Reporting)
library(hidiTS)
library(fbi)
library(ggplot2)
library("reshape2")
library(ggpubr)


rm(list = ls())

#-------------------Load and prepare Data -----------------------------------------

prep_data_hidi <- function(data){

  data_hidi <- as.matrix(t(data))


  output <- list("data_fbi"= data, "X"=data_hidi )

  return(output)
}


df_paper=fredqd(date_end =as.Date(c("1985/12/01")), date_start=as.Date(c("1960/03/01")),transform = TRUE)

# remove dates
df_temp=df_paper[,2:(length(df_paper))]

#drop Nas
Na_names <- colnames(df_temp)[colSums(is.na(df_temp)) > 0]
df_colnames <- colnames(df_temp)
col_keeps <- setdiff(df_colnames,Na_names)
df_no_na <- subset(df_temp, select=(col_keeps))

df_final <- prep_data_hidi(df_no_na)


X <- df_final$X

T <- ncol(X)
n <- nrow(X)


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
  print(paste(' r=',r, sep=''))
}

df <- cbind(ic_BIC_n, ic_BIC_T, ic_BIC_nT)[1:200,]

colnames(df) <- c('BIC_n', 'BIC_T', 'BIC_nT')

df_long <- melt(df)
colnames(df_long) <- c('r', 'IC', 'value')
(figure1 <- ggplot(df_long, aes(x=r, y=value, color=IC)) + geom_line() +
  theme(text = element_text(size=28), legend.position="bottom") + labs(color='') )
ggsave("static/application_BICs_full.png")

df_long2 <- melt(df[1:20,])
colnames(df_long2) <- c('r', 'IC', 'value')
(figure2 <- ggplot(df_long2, aes(x=r, y=value, color=IC)) + geom_line() + geom_point() +
    theme(text = element_text(size=28), legend.position="bottom") + labs(color='') )
ggsave("static/application_BICs_n20.png")

(figure3 <- ggplot(as.data.frame(cbind('BNIC'=ic_BNIC, 'r'=1:n)), aes(x=r, y=BNIC, lty = 'BNIC')) +
  geom_line()  + theme(text = element_text(size=28), legend.position="bottom") + ylab('value') + scale_linetype(''))
ggsave("static/application_Bai_full.png")

(figure4 <- ggplot(as.data.frame(cbind('BNIC'=ic_BNIC, 'r'=1:n)[1:20,]), aes(x=r, y=BNIC, lty = 'BNIC')) +
  geom_line() + geom_point() + theme(text = element_text(size=28), legend.position="bottom") + ylab('value')+ scale_linetype(''))
ggsave("static/application_Bai_n20.png")











