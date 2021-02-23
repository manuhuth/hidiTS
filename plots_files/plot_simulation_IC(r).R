library(ggplot2)
library("reshape2")
library(ggpubr)


load("simulated_data/ics.RData")
alpha = 0.1
up <- 1 - alpha/2
low <-  alpha/2

mean_n <- colMeans(BIC_n_df)
high_n <- apply(BIC_n_df, 2, quantile, up)
low_n <- apply(BIC_n_df, 2, quantile, low)
BIC_n_df_plot <- cbind('low'=low_n, 'mean'=mean_n, 'high'=high_n, 'r'= 1:n, 'IC' = rep(1, n) )

mean_T <- colMeans(BIC_T_df)
high_T <- apply(BIC_T_df, 2, quantile, up)
low_T <- apply(BIC_T_df, 2, quantile, low)
BIC_T_df_plot <-cbind('low'=low_T, 'mean'=mean_T, 'high'=high_T, 'r'= 1:n, 'IC' = rep(2, n) )


mean_nT <- colMeans(BIC_nT_df)
high_nT <- apply(BIC_nT_df, 2, quantile, up)
low_nT <- apply(BIC_nT_df, 2, quantile, low)
BIC_nT_df_plot <- cbind('low'=low_nT, 'mean'=mean_nT, 'high'=high_nT, 'r'= 1:n, 'IC' = rep(3, n) )

mean_bai <- colMeans(BNIC_df)
high_bai <- apply(BNIC_df, 2, quantile, up)
low_bai <- apply(BNIC_df, 2, quantile, low)
BNIC_df_plot <- as.data.frame(cbind('low'=low_bai, 'mean'=mean_bai, 'high'=high_bai, 'r'= 1:n))



df <- as.data.frame(rbind(BIC_n_df_plot[1:200,], BIC_T_df_plot[1:200,], BIC_nT_df_plot[1:200,]))
df['IC'][df['IC'] == 1] <- 'BIC_n'
df['IC'][df['IC'] == 2] <- 'BIC_T'
df['IC'][df['IC'] == 3] <- 'BIC_nT'
#df['IC'][df['IC'] == 4] <- 'BNIC'

(fig1 <- ggplot(df, aes(x=r, y=mean,  color=IC)) + geom_line() +
  geom_ribbon(aes(x=r, ymin=low, ymax=high, color=IC),linetype=2, alpha = 0.2) +
  scale_color_manual(
    values = c(BIC_n="#F8766D", BIC_nT="#619CFF", BIC_T="#00BA38"))  +
  theme(text = element_text(size=28), legend.position="bottom") + labs(color='') ) + ylab('Value')
ggsave("static/simulation_BICs_full.png")



df2 <- as.data.frame(rbind(BIC_n_df_plot[1:20,], BIC_T_df_plot[1:20,], BIC_nT_df_plot[1:20,]))
df2['IC'][df2['IC'] == 1] <- 'BIC_n'
df2['IC'][df2['IC'] == 2] <- 'BIC_T'
df2['IC'][df2['IC'] == 3] <- 'BIC_nT'
#df2['IC'][df2['IC'] == 4] <- 'BNIC'

(fig2 <- ggplot(df2, aes(x=r, y=mean,  color=IC)) + geom_line() +
    geom_ribbon(aes(x=r, ymin=low, ymax=high, color=IC),linetype=2, alpha = 0.2) +
    scale_color_manual(
      values = c(BIC_n="#F8766D", BIC_nT="#619CFF", BIC_T="#00BA38"))  +
    theme(text = element_text(size=28), legend.position="bottom") + labs(color='') ) + ylab('Value')
ggsave("static/simulation_BICs_n20.png")


(fig3 <- ggplot(BNIC_df_plot[1:200,], aes(x=r, y=mean,  lty = 'BNIC')) + geom_line() +
    geom_ribbon(aes(x=r, ymin=low, ymax=high),linetype=2, alpha = 0.2) +
    theme(text = element_text(size=28), legend.position="bottom") + scale_linetype('')  + ylab('Value'))
ggsave("static/simulation_BNIC_full.png")

df_BNIC_n20 <- BNIC_df_plot[1:20,]
(fig4 <- ggplot(df_BNIC_n20, aes(x=r, y=mean,  lty = 'BNIC')) + geom_line() +
    geom_ribbon(aes(x=r, ymin=low, ymax=high),linetype=2, alpha = 0.2) +
    theme(text = element_text(size=28), legend.position="bottom") + scale_linetype('')  + ylab('Value'))
ggsave("static/simulation_BNIC_n20.png")






