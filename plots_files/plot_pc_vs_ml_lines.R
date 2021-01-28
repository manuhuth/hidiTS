library(ggplot2)
library("reshape2")
library(ggpubr)
load("C:/Users/Mhuth/Desktop/hidiTS/simulated_data/Lambda1point2_with_ML.RData")


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

(figure1 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line())
ggsave("C:/Users/Mhuth/Desktop/hidiTS/static/plot_mse_vs_pca_factor1.png")

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

(figure2 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line())
ggsave("C:/Users/Mhuth/Desktop/hidiTS/static/plot_mse_vs_pca_factor2.png")

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

(figure3 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line())
ggsave("C:/Users/Mhuth/Desktop/hidiTS/static/plot_mse_vs_pca_factor3.png")
