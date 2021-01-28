#### plot difference between MSE of Maximum Likelihood and PCA
library(tidyr)
library(plotly)
load("hidiTS/simulated_data/Lambda1point2_only_pca.RData")

maxn <- 80
maxt <- 80

df <- simulated_data[c('pca_bai_right','n', 'T')]
df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

ns <- unique(df$n )
ns <- ns[ns <= maxn]
Ts <- unique(df$T )
Ts <- Ts[Ts <= maxt]

df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
df_wide <- subset(df_wide, select=-c(n))

fig <- plot_ly(df_wide, x = ns, y = Ts, z = as.matrix(df_wide))
fig <- fig %>% add_surface()
fig<- fig %>% layout(title = "", scene = list(xaxis = list(title = "n"),
  yaxis = list(title = "T"),
  zaxis = list(title = "Bai & Ng IC")
))
fig


#-------------BIC------------------------------------------------
maxn <- 150
maxt <- 150
df <- simulated_data[c('pca_bic_right','n', 'T')]
df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

ns <- unique(df$n )
ns <- ns[ns <= maxn]
Ts <- unique(df$T )
Ts <- Ts[Ts <= maxt]

df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
df_wide <- subset(df_wide, select=-c(n))
#####

fig <- plot_ly(df_wide, x = ns, y = Ts, z = as.matrix(df_wide))
fig <- fig %>% add_surface(showscale = FALSE)
fig<- fig %>% layout(title = "", scene = list(xaxis = list(title = "n"),
                                              yaxis = list(title = "T"),
                                              zaxis = list(title = "BIC")))
fig
