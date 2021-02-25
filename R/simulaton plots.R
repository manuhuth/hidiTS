accuracy_plots <- function(simulated_data, maxn=300, maxt=300) {


    df <- simulated_data[c('pca_bai_right','n', 'T')]
    df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

    ns <- unique(df$n )
    ns <- ns[ns <= maxn]
    Ts <- unique(df$T )
    Ts <- Ts[Ts <= maxt]

    df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
    df_wide <- subset(df_wide, select=-c(n))

    fig_bai <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
    fig_bai <- fig_bai %>% add_surface(showscale = TRUE)

    fig_bai <- fig_bai %>% layout(title = "", scene = list(xaxis = list(title =list(text=" n ", font=list(size=35))),
                                                  yaxis = list(title = list(text= "T", font=list(size=35))),
                                                  zaxis = list(title =list(text="Accuracy", font=list(size=32)))
    ), font= list(family = "sans serif", size = 20) )


    #-------------BIC_n------------------------------------------------
    df <- simulated_data[c('pca_bic_right','n', 'T')]
    df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

    ns <- unique(df$n )
    ns <- ns[ns <= maxn]
    Ts <- unique(df$T )
    Ts <- Ts[Ts <= maxt]

    df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
    df_wide <- subset(df_wide, select=-c(n))
    #####

    fig_n <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
    fig_n <- fig_n %>% add_surface(showscale = TRUE)
    fig_n <- fig_n %>% layout(title = "", scene = list(xaxis = list(title =list(text="n", font=list(size=35))),
                                                   yaxis = list(title = list(text= "T", font=list(size=35))),
                                                   zaxis = list(title = list(text="Accuracy", font=list(size=32)))
    ), font= list(family = "sans serif", size = 20) )


    #------------------BIC_t -------------------------------------------
    df <- simulated_data[c('pca_bic_T_right','n', 'T')]
    df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

    ns <- unique(df$n )
    ns <- ns[ns <= maxn]
    Ts <- unique(df$T )
    Ts <- Ts[Ts <= maxt]

    df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
    df_wide <- subset(df_wide, select=-c(n))
    #####

    fig_t <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
    fig_t <- fig_t %>% add_surface(showscale = TRUE)
    fig_t <- fig_t %>% layout(title = "", scene= list(xaxis = list(title =list(text="n", font=list(size=35))),
                                                 yaxis = list(title = list(text= "T", font=list(size=35))),
                                                 zaxis = list(title = list(text="Accuracy", font=list(size=32)))
    ), font= list(family = "sans serif", size = 20) )



    #---------------------------BIC_nt-------------------------------
    df <- simulated_data[c('pca_bic_nT_right','n', 'T')]
    df <- df[which( (df$T <= maxt) & (df$n <=maxn) ),]

    ns <- unique(df$n )
    ns <- ns[ns <= maxn]
    Ts <- unique(df$T )
    Ts <- Ts[Ts <= maxt]

    df_wide <- reshape(df, idvar = "n", timevar = "T", direction = "wide")
    df_wide <- subset(df_wide, select=-c(n))
    #####

    fig_nT <- plot_ly(df_wide, x = ns, y = Ts, z = t(as.matrix(df_wide)))
    fig_nT <- fig_nT %>% add_surface(showscale = TRUE)
    fig_nT <- fig_nT %>% layout(title = "", scene = list(xaxis = list(title =list(text="n", font=list(size=35))),
                                                  yaxis = list(title = list(text= "T", font=list(size=35))),
                                                  zaxis = list(title = list(text="Accuracy", font=list(size=32)))
    ), font= list(family = "sans serif", size = 20) )

    return(list('BNIC' = fig_bai, 'BIC_nT' = fig_nT, 'BIC_n' = fig_n, 'BIC_T' = fig_t))

}
