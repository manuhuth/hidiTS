library(ggplot2)
library("reshape2")
library(ggpubr)
#load("C:/Users/Mhuth/Desktop/hidiTS/simulated_data/Lambda1point2_only_pca.RData")

Ts <- unique(simulated_data$T)
Ts <- c(7,15,30, 50, 100, 300)

df <- simulated_data[c('mse_pca_f','n', 'T' )]
colnames(df) <- c('MSE', 'n', 'T')
df <- df[ which( df$T %in% Ts) , ]

df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
names <- c()
for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

colnames(df) <- c('n', names)

df_long <- melt(df, id = 'n')
colnames(df_long) <- c('n', 'Periods', 'MSE')

figure1 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line() + ggtitle("Average of 1000 Samples")
ggsave("C:/Users/McMüll/Desktop/plot_mse_factor1.png")


#----------fig3-------------------------------------------------
df <- simulated_data[c('mse_pca_f.1','n', 'T' )]
colnames(df) <- c('MSE', 'n', 'T')
df <- df[ which( df$T %in% Ts) , ]

df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
names <- c()
for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

colnames(df) <- c('n', names)

df_long <- melt(df, id = 'n')
colnames(df_long) <- c('n', 'Periods', 'MSE')

figure2 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line() +  ggtitle("Average of 1000 Samples")
ggsave("C:/Users/McMüll/Desktop/plot_mse_factor2.png")

#----------fig3-------------------------------------------------
df <- simulated_data[c('mse_pca_f.2','n', 'T' )]
colnames(df) <- c('MSE', 'n', 'T')
df <- df[ which( df$T %in% Ts) , ]

df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
names <- c()
for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

colnames(df) <- c('n', names)

df_long <- melt(df, id = 'n')
colnames(df_long) <- c('n', 'Periods', 'MSE')

figure3 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) + geom_line() +  ggtitle("Average of 1000 Samples")
ggsave("C:/Users/McMüll/Desktop/plot_mse_factor3.png")