#-------------------------------------------------------------------------------------------------#
#                                                                                                 #
#                                           Simulation Study                                      #
#                                                                                                 #
#-------------------------------------------------------------------------------------------------#


# Load Packages
#-------------------------------------------------------------------------------------------------
install.packages("devtools") #load/install
devtools::install_github("manuhuth/hidiTS")

library(hidiTS)

library(optimParallel) #load/install
library(ggplot2) #load/install
library(xtable) #load/install
library(reshape2) #load/install
library(ggpubr) #load/install
library(plotly) #load/install
library(tidyr) #load/install


# Prepare Multi Core Work
#-------------------------------------------------------------------------------------------------
cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl = cl)



# 1. -  Data Simulation Difference MSE PCA vs. ML - 500 Samples
#-------------------------------------------------------------------------------------------------

    # 1.1 - b = 0.5, Set-up as described in Table 1
    #-------------------------------------------------------------------------------------------------
    #link to function
    df_table_1_1 <- simulation_mse_vs_pca(b=0.5, n_simulation=c(5, 15, 20, 30), T_simulation=c(7, 15, 20, 30, 40), number_samples = 500, start_seed = 1)



    # 1.2 - b = 3^0.5, Set-up as described in Table 1
    #-------------------------------------------------------------------------------------------------
    df_table_1_2 <-simulation_mse_vs_pca(b=3^0.5, n_simulation=c(5, 15, 20, 30), T_simulation=c(7, 15, 20, 30, 40), number_samples = 500, start_seed = 1)




# 2. -  Data Simulation PCA Simulation Study - 1000 Samples
#-------------------------------------------------------------------------------------------------


    # 2.1 b = 0.5, Set-up as described in Figure 3
    #-------------------------------------------------------------------------------------------------
    df_pca_b_low <- simulation_pca(b=0.5, T_simulation = c(7,15,30,40,seq(50,300,50)), n_simulation = c(5,10,seq(10,40,10),seq(50,300,50)), number_samples = 2, start_seed = 1)



    # 2.2 b = 3^0.5, Set-up as described in Figure 3
    #-------------------------------------------------------------------------------------------------
    df_pca_b_medium <- simulation_pca(b=3^0.5, T_simulation = c(7,15,30,40,seq(50,300,50)), n_simulation = c(5,10,seq(10,40,10),seq(50,300,50)), number_samples = 2, start_seed = 1)



# 3. - Tables and Figures
#-------------------------------------------------------------------------------------------------


    # 3.1 - Table 1 - MSEs PCA vs. ML
    #-------------------------------------------------------------------------------------------------
    df_table_1_2['mse_diff_f'] <- df_table_1_2$mse_pca_f-df_table_1_2$mse_ml_f
    df_table_1_2['mse_diff_f.1'] <- df_table_1_2$mse_pca_f.1-df_table_1_2$mse_ml_f.1
    df_table_1_2['mse_diff_f.2'] <- df_table_1_2$mse_pca_f.2-df_table_1_2$mse_ml_f.2

    sim_data_medium <- df_table_1_2[!duplicated(df_table_1_1[c('n','T')]),]

    ns <- 1:400
    Ts <- 1:400

    data_medium <- sim_data_medium[ which(sim_data_medium$T %in% Ts & sim_data_medium$n %in% ns), ]
    ics_medium <- data_medium[,c('mse_diff_f', 'mse_diff_f.1', 'mse_diff_f.2')]
    names(ics_medium) <- c('factor1', 'factor2', 'factor3')

    df_table_1_1['mse_diff_f'] <- df_table_1_1$mse_pca_f-df_table_1_1$mse_ml_f
    df_table_1_1['mse_diff_f.1'] <- df_table_1_1$mse_pca_f.1-df_table_1_1$mse_ml_f.1
    df_table_1_1['mse_diff_f.2'] <- df_table_1_1$mse_pca_f.2-df_table_1_1$mse_ml_f.2

    sim_data_small <- df_table_1_1[!duplicated(df_table_1_2[c('n','T')]),]

    data_small <- sim_data_small[ which(sim_data_small$T %in% Ts & sim_data_small$n %in% ns), ]
    ics_small<- data_small[,c('n', 'T','mse_diff_f', 'mse_diff_f.1', 'mse_diff_f.2')]
    names(ics_small) <- c('n', 'T','factor1', 'factor2', 'factor3')

    data_final <- cbind(ics_small, ics_medium)

    data_final['n'] <- as.integer(data_final[,'n'])
    data_final['T'] <- as.integer(data_final[,'T'])

    print(xtable(data_final, type = "latex", digits=c(0,0,0,4,4,4, 4,4,4) ), include.rownames=FALSE)


    # 3.2 - Figure 2 - Eigenvalue Structure of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------

    df <- df_pca_b_low[c("ev1","ev2" ,"ev3","ev4", "ev5", 'n', 'T' )]
    for (index in 1:5) {
      df[c(paste('ev',index, sep=''))] <- df[c(paste('ev',index, sep=''))] / df['n']
    }

    ns <- c(5, 20, 50, 100)
    Ts <- c(7, 15, 30, 100)
    width <- 0.6

    list_df <- list()

    for (index in 1:length(Ts)) {
      list_df[[index]] <- melt(df[ which( (df$T == Ts[index]) & (df$n %in% ns) ),
                                   c("ev1","ev2" ,"ev3","ev4", "ev5", 'n')], id = 'n')
      list_df[[index]]$n <- as.factor(list_df[[index]]$n)
      colnames(list_df[[index]]) <- c('n', 'Eigenvalue','Magnitude')
    }

    t7 <- ggplot(list_df[[1]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
      geom_bar(stat="identity", width = width, position=position_dodge())+
      scale_fill_brewer(palette="Paired") + ggtitle('T=7') + theme(text = element_text(size=22))

    #t15 <- ggplot(list_df[[2]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
    #  geom_bar(stat="identity",width = width, position=position_dodge())+
    #  scale_fill_brewer(palette="Paired") + ggtitle('T=15') + theme(text = element_text(size=20)) + ylab('')

    t30 <- ggplot(list_df[[3]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
      geom_bar(stat="identity",width = width, position=position_dodge())+
      scale_fill_brewer(palette="Paired") + ggtitle('T=30') + theme(text = element_text(size=22))
    #t100 <- ggplot(list_df[[4]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
    #  geom_bar(stat="identity",width = width, position=position_dodge())+
    #  scale_fill_brewer(palette="Paired") + ggtitle('T=100') + theme(text = element_text(size=20)) + ylab('')

    (figure <- ggarrange(t7,t30, ncol = , nrow=2, common.legend = TRUE, legend="bottom"))


    df <- df_pca_b_medium[c("ev1","ev2" ,"ev3","ev4", "ev5", 'n', 'T' )]
    for (index in 1:5) {
      df[c(paste('ev',index, sep=''))] <- df[c(paste('ev',index, sep=''))] / df['n']
    }

    ns <- c(5, 20, 50, 100)
    Ts <- c(7, 15, 30, 100)
    width <- 0.6

    list_df <- list()

    for (index in 1:length(Ts)) {
      list_df[[index]] <- melt(df[ which( (df$T == Ts[index]) & (df$n %in% ns) ),
                                   c("ev1","ev2" ,"ev3","ev4", "ev5", 'n')], id = 'n')
      list_df[[index]]$n <- as.factor(list_df[[index]]$n)
      colnames(list_df[[index]]) <- c('n', 'Eigenvalue','Magnitude')
    }

    t7_med <- ggplot(list_df[[1]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
      geom_bar(stat="identity", width = width, position=position_dodge())+
      scale_fill_brewer(palette="Paired") + ggtitle('T=7') + theme(text = element_text(size=22))

    #t15 <- ggplot(list_df[[2]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
    #  geom_bar(stat="identity",width = width, position=position_dodge())+
    #  scale_fill_brewer(palette="Paired") + ggtitle('T=15') + theme(text = element_text(size=20)) + ylab('')

    t30_med <- ggplot(list_df[[3]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
      geom_bar(stat="identity",width = width, position=position_dodge())+
      scale_fill_brewer(palette="Paired") + ggtitle('T=30') + theme(text = element_text(size=22))
    #t100 <- ggplot(list_df[[4]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
    #  geom_bar(stat="identity",width = width, position=position_dodge())+
    #  scale_fill_brewer(palette="Paired") + ggtitle('T=100') + theme(text = element_text(size=20)) + ylab('')

    (figure_med <- ggarrange(t7_med,t30_med, ncol = , nrow=2, common.legend = TRUE, legend="bottom"))
    #ggsave("static/eigenvalues_medium_lambda.png")


    # 3.3 - Figure 7 & Figure 8 - Average Accuracy of the BIC_nT and the BNIC of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------
    list_accuracy_plots_low_b <- accuracy_plots(df_pca_b_low, maxn=300, maxt=300)
    list_accuracy_plots_medium_b <- accuracy_plots(df_pca_b_medium, maxn=300, maxt=300)
    list_accuracy_plots_low_b[[1]] #3D plots can be rotated manually
    list_accuracy_plots_low_b[[2]]
    list_accuracy_plots_low_b[[3]]
    list_accuracy_plots_low_b[[4]]

    list_accuracy_plots_medium_b[[1]]
    list_accuracy_plots_medium_b[[2]]
    list_accuracy_plots_medium_b[[3]]
    list_accuracy_plots_medium_b[[4]]


    # 3.4 - Table 3 - Average Accuracy of the BIC_nT and the BNIC of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------
    sim_data_medium <- df_pca_b_medium[!duplicated(df_pca_b_medium[c('n','T')]),]

    ns <- c(5, 10, 30, 50, 100,300)
    Ts <- c(7, 15, 30, 50, 100, 300)

    data_medium <- sim_data_medium[ which(sim_data_medium$T %in% Ts
                                          & sim_data_medium$n %in% ns), ]

    ics_medium <- data_medium[,c('pca_bai_right', 'pca_bic_right', 'pca_bic_T_right', 'pca_bic_nT_right')]
    names(ics_medium) <- c('BNIC_med', 'BIC_n_med', 'BIC_T_med', 'BIC_nT_med')

    sim_data_small <- simulated_small[!duplicated(df_pca_b_small[c('n','T')]),]

    data_small <- sim_data_small[ which(sim_data_small$T %in% Ts
                                        & sim_data_small$n %in% ns), ]
    ics_small<- data_small[,c('n', 'T','pca_bai_right', 'pca_bic_right', 'pca_bic_T_right', 'pca_bic_nT_right')]
    names(ics_small) <- c('n', 'T','BNIC_small', 'BIC_n_small', 'BIC_T_small', 'BIC_nT_small')

    data_final <- cbind(ics_small, ics_medium)

    data_final['n'] <- as.integer(data_final[,'n'])
    data_final['T'] <- as.integer(data_final[,'T'])

    print(xtable(data_final, type = "latex", digits=c(0,0,0,2,2,2,2,2, 2,2,2)), include.rownames=FALSE)


    # 3.5 - Figure 3 & Figure 9 - MSEs of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------
    list_MSE_plots_low_b <- MSE_plots(df_pca_b_low, Ts = c(7,15,30, 50, 100,300))
    list_MSE_plots_medium_b <- MSE_plots(df_pca_b_medium, Ts = c(7,15,30, 50, 100,300))
    list_MSE_plots_low_b[[1]]
    list_MSE_plots_low_b[[2]]
    list_MSE_plots_low_b[[3]]

    list_MSE_plots_medium_b[[1]]
    list_MSE_plots_medium_b[[2]]
    list_MSE_plots_medium_b[[3]]


    # 3.6 - Table 4 - MSEs of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------
    sim_data_medium <- df_pca_b_medium[!duplicated(df_pca_b_medium[c('n','T')]),]

    ns <- c(5, 10, 30, 50, 100, 300)
    Ts <- c(7, 15, 30, 50, 100, 300)

    data_medium <- sim_data_medium[ which(sim_data_medium$T %in% Ts & sim_data_medium$n %in% ns), ]

    ics_medium <- data_medium[,c('mse_pca_f', 'mse_pca_f.1', 'mse_pca_f.2')]
    names(ics_medium) <- c('factor1', 'factor2', 'factor3')

    sim_data_small <- df_pca_b_low[!duplicated(df_pca_b_low[c('n','T')]),]

    data_small <- sim_data_small[ which(sim_data_small$T %in% Ts & sim_data_small$n %in% ns), ]
    ics_small<- data_small[,c('n', 'T','mse_pca_f', 'mse_pca_f.1', 'mse_pca_f.2')]
    names(ics_small) <- c('n', 'T','factor1', 'factor2', 'factor3')

    data_final <- cbind(ics_small, ics_medium)

    data_final['n'] <- as.integer(data_final[,'n'])
    data_final['T'] <- as.integer(data_final[,'T'])

    print(xtable(data_final, type = "latex", digits=c(0,0,0,4,4,4, 4,4,4) ), include.rownames=FALSE)


    # 3.7 - Figure 5 - Time ML vs PCA
    #-------------------------------------------------------------------------------------------------
    df_Figure_5 <-simulation_mse_vs_pca(b=3^0.5, n_simulation=c(seq(5, 125,10)),
                                        T_simulation=c(20), number_samples = 20, start_seed = 1)
    df_help <- df_Figure_5[df_Figure_5$T==20,]
    df_plot <- df_help[c('time_pca', 'time_ml', 'n')]
    colnames(df_plot) <- c('PCA', 'ML', 'n')
    df_plot_long <- melt(df_plot, id='n')
    colnames(df_plot_long) <- c('n', 'Method', 'Time')


    ggplot(df_plot_long, aes(x=n, y=Time, colour=Method)) + geom_line() +
      ggtitle(' ',
              subtitle = 'q = 3, T = 20, maximal optimizer iterations = 5, parallelized optimization') +
      labs(fill = "Methods")


    # 3.8 - Figure 6 - Example Eigenvalue Structure
    #-------------------------------------------------------------------------------------------------
    n <- 20
    T <- 20

    set.seed(123)
    data <- sim_data(p = n, T = T, dim_F= 3, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -3, up_X = 3,
                     low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = diag(n),
                     adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)


    vec <- c()
    vec_bic <- c()
    for (k in 1:20) {
      pc <- pca_estimator(data$X, k)
      x_hat <- pc$Lambda %*% pc$F
      rss <- mean((data$X - x_hat)^2)

      bic <- n * log(rss/(n*T)) + k * log(n)
      bai <- log(rss/(n*T)) + k *(n+T)/(n*T) * log(min(n,T))
      vec <- cbind(vec, bai)
      vec_bic <- c(vec_bic, bic)
    }

    data_smooth <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -1, up_X = 1,
                            low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = 10*diag(n),
                            adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

    df <- data.frame(cbind(eigen(var(t(data$X)))$values,eigen(var(t(data_smooth$X)))$values, 1:length(eigen(var(t(data$X)))$values)))
    colnames(df) <- c('Sharp Cut-Off','Smooth', 'Eigenvalue')

    df_long <- melt(df, id='Eigenvalue')
    colnames(df_long) <- c('Eigenvalue', 'Shape', 'Magnitude')

    ggplot(df_long, aes(x=Eigenvalue, y=Magnitude, colour=Shape)) + geom_point(aes(size = 2)) +
      theme(text = element_text(size=28))













