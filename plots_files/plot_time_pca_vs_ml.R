library(ggplot2)
library("reshape2")

T <- 20

df_help <- simulated_data[simulated_data$T==T,]
df_plot <- df_help[c('time_pca', 'time_ml', 'n')]
colnames(df_plot) <- c('PCA', 'ML', 'n')
df_plot_long <- melt(df_plot, id='n')
colnames(df_plot_long) <- c('n', 'Method', 'Time')


ggplot(df_plot_long, aes(x=n, y=Time, colour=Method)) + geom_line() +
  ggtitle(' ',
          subtitle = 'q = 3, t = 20, maximal optimizer iterations = 5, parallelized optimization') +
  labs(fill = "Methods")


n_val <- df_plot['n']
par <- 3*(n_val + T) + n_val

df_par <- data.frame(cbind(n_val, par))
colnames(df_par) <- c('n', 'Parameter')

ggplot(df_par, aes(x=n, y=Parameter)) + geom_line() +
  ggtitle('Number of parameter to estimate', subtitle = 'q=3, t = 20') 