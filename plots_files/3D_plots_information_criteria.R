#### plot difference between MSE of Maximum Likelihood and PCA
library(tidyr)
library(plotly)
load("C:/Users/katri/Desktop/hidiTS/raspberry_paper_medium_lambda_ext.RData")


#mse_diff_f <- as.matrix(abs(simulated_data$mse_pca_f - simulated_data$mse_ml_f))
#mse_diff_f.1<- abs(simulated_data$mse_pca_f.1 - simulated_data$mse_ml_f.1)
#mse_diff_f.2<- abs(simulated_data$mse_pca_f.2 - simulated_data$mse_ml_f.2)

#simulated_data$mse_diff_f <- mse_diff_f
#simulated_data$mse_diff_f.1 <- mse_diff_f.1
#simulated_data$mse_diff_f.2 <- mse_diff_f.2

#n <- 150, y = 100
#-------------------------- BNIC-------------
maxn <- 300
maxt <- 300

df <- simulated_data[c('pca_bai_right','n', 'T')]
df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

ns <- unique(df$n )
ns <- ns[ns <= maxn]
Ts <- unique(df$T )
Ts <- Ts[Ts <= maxt]

df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
df_wide <- subset(df_wide, select=-c(n))

fig <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
fig <- fig %>% add_surface(showscale = FALSE)

fig<- fig %>% layout(title = "", scene = list(xaxis = list(title =list(text=" n ", font=list(size=35))),
                                              yaxis = list(title = list(text= "T", font=list(size=35))),
                                              zaxis = list(title =list(text="Accuracy", font=list(size=32)))
), font= list(family = "sans serif", size = 20) )
fig


#-------------BIC_n------------------------------------------------
maxn <- 300
maxt <- 300
df <- simulated_data[c('pca_bic_right','n', 'T')]
df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

ns <- unique(df$n )
ns <- ns[ns <= maxn]
Ts <- unique(df$T )
Ts <- Ts[Ts <= maxt]

df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
df_wide <- subset(df_wide, select=-c(n))
#####

fig <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
fig <- fig %>% add_surface(showscale = FALSE)
fig<- fig %>% layout(title = "", scene = list(xaxis = list(title =list(text="n", font=list(size=35))),
                                              yaxis = list(title = list(text= "T", font=list(size=35))),
                                              zaxis = list(title = list(text="Accuracy", font=list(size=32)))
), font= list(family = "sans serif", size = 20) )

fig

#------------------BIC_t -------------------------------------------
maxn <- 300
maxt <- 300
df <- simulated_data[c('pca_bic_T_right','n', 'T')]
df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

ns <- unique(df$n )
ns <- ns[ns <= maxn]
Ts <- unique(df$T )
Ts <- Ts[Ts <= maxt]

df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
df_wide <- subset(df_wide, select=-c(n))
#####

fig <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
fig <- fig %>% add_surface(showscale = FALSE)
fig<- fig %>% layout(title = "", scene= list(xaxis = list(title =list(text="n", font=list(size=35))),
                                      yaxis = list(title = list(text= "T", font=list(size=35))),
                                      zaxis = list(title = list(text="Accuracy", font=list(size=32)))
), font= list(family = "sans serif", size = 20) )
fig



#---------------------------BIC_nt-------------------------------
maxn <- 300
maxt <- 300
df <- simulated_data[c('pca_bic_nT_right','n', 'T')]
df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

ns <- unique(df$n )
ns <- ns[ns <= maxn]
Ts <- unique(df$T )
Ts <- Ts[Ts <= maxt]

df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
df_wide <- subset(df_wide, select=-c(n))
#####

fig <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
fig <- fig %>% add_surface(showscale = FALSE)
fig<- fig %>% layout(title = "", scene = list(xaxis = list(title =list(text="n", font=list(size=35))),
                                              yaxis = list(title = list(text= "T", font=list(size=35))),
                                              zaxis = list(title = list(text="Accuracy", font=list(size=32)))
), font= list(family = "sans serif", size = 20) )
fig
