#install.packages("ggpubr")
#install.packages("reshape2")
#install.packages("foreign")
library(foreign)
library(ggplot2)
library(reshape2)
library(ggpubr)
#load("C:/Users/Mhuth/Desktop/hidiTS/simulated_data/Lambda1point2_with_ML.RData")
simulated_data <- readRDS(file = "blueberry_medium_lambda_complete.rds")

simulated_data['mse_diff_f'] <- simulated_data$mse_pca_f-simulated_data$mse_ml_f
simulated_data['mse_diff_f.1'] <- simulated_data$mse_pca_f.1-simulated_data$mse_ml_f.1
simulated_data['mse_diff_f.2'] <- simulated_data$mse_pca_f.2-simulated_data$mse_ml_f.2

Ts <- unique(simulated_data$T)

df <- simulated_data[c('mse_diff_f','n', 'T' )]
colnames(df) <- c('MSE', 'n', 'T')
df <- df[ which( df$T %in% Ts) , ]

df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
names <- c()
for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

colnames(df) <- c('n', names)

df_long <- melt(df, id = 'n')
colnames(df_long) <- c('n', 'Periods', 'MSE')
(figure1 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line() + labs(y=paste("Delta MSE",sep = "")) + ggtitle("Average of 500 Samples")
  +expand_limits(y=c(-0.05,0, 0.05))) + scale_y_continuous(breaks = c(-0.05,0,0.05))
ggsave("C:/Users/McMüll/Desktop/plot_mse_vs_pca_factor1.png")

#-----------fig2-----------------------------
df <- simulated_data[c('mse_diff_f.1','n', 'T' )]
colnames(df) <- c('MSE', 'n', 'T')
df <- df[ which( df$T %in% Ts) , ]

df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
names <- c()
for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

colnames(df) <- c('n', names)

df_long <- melt(df, id = 'n')
colnames(df_long) <- c('n', 'Periods', 'MSE')

(figure2 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line() + labs(y=paste("Delta MSE",sep = ""))) + ggtitle("Average of 500 Samples")
ggsave("C:/Users/McMüll/Desktop/plot_mse_vs_pca_factor2.png")

#-----------fig3-----------------------------
df <- simulated_data[c('mse_diff_f.2','n', 'T' )]
colnames(df) <- c('MSE', 'n', 'T')
df <- df[ which( df$T %in% Ts) , ]

df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
names <- c()
for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

colnames(df) <- c('n', names)

df_long <- melt(df, id = 'n')
colnames(df_long) <- c('n', 'Periods', 'MSE')

(figure3 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line() + labs(y=paste("Delta MSE",sep = ""))) + ggtitle("Average of 500 Samples")
ggsave("C:/Users/McMüll/Desktop/plot_mse_vs_pca_factor3.png")
© 2021 GitHub, Inc.