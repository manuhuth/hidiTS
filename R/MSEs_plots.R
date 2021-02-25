MSE_plots <- function(simulated_data, Ts = c(7,15,30, 50, 100,300)) {

  df <- simulated_data[c('mse_pca_f','n', 'T' )]
  colnames(df) <- c('MSE', 'n', 'T')
  df <- df[ which( df$T %in% Ts) , ]

  df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
  names <- c()
  for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

  colnames(df) <- c('n', names)

  df_long <- melt(df, id = 'n')
  colnames(df_long) <- c('n', 'Periods', 'MSE')

  figure1 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) +theme(axis.text=element_text(size=20),
                                                                    axis.title=element_text(size=20),
                                                                    legend.title = element_text(size = 20),
                                                                    legend.text = element_text(size = 20),
                                                                    legend.position="bottom") + geom_line()

  df <- simulated_data[c('mse_pca_f.1','n', 'T' )]
  colnames(df) <- c('MSE', 'n', 'T')
  df <- df[ which( df$T %in% Ts) , ]

  df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
  names <- c()
  for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

  colnames(df) <- c('n', names)

  df_long <- melt(df, id = 'n')
  colnames(df_long) <- c('n', 'Periods', 'MSE')

  figure2 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) +theme(axis.text=element_text(size=20),
                                                                    axis.title=element_text(size=20),
                                                                    legend.title = element_text(size = 20),
                                                                    legend.text = element_text(size = 20),
                                                                    legend.position="bottom") + geom_line()


  df <- simulated_data[c('mse_pca_f.2','n', 'T' )]
  colnames(df) <- c('MSE', 'n', 'T')
  df <- df[ which( df$T %in% Ts) , ]

  df <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
  names <- c()
  for (index in 1:length(Ts)) {names[index] = paste('T=',Ts[index],sep='')}

  colnames(df) <- c('n', names)

  df_long <- melt(df, id = 'n')
  colnames(df_long) <- c('n', 'Periods', 'MSE')

  figure3 <- ggplot(df_long, aes(x=n, y=MSE, color=Periods)) +theme(axis.text=element_text(size=20),
                                                                    axis.title=element_text(size=20),
                                                                    legend.title = element_text(size = 20),
                                                                    legend.text = element_text(size = 20),
                                                                    legend.position="bottom") + geom_line()



  return(list('factor1' = figure1, 'factor2' = figure2, 'factor3' = figure3))
}
