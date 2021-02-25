#-------------------------------------------------------------------------------------------------#
#                                                                                                 #
#                                     Replicate the Paper                                         #
#                                                                                                 #
#-------------------------------------------------------------------------------------------------#


# Install Packages if Necessary
#-------------------------------------------------------------------------------------------------
list_packages <- c("optimParallel", "ggplot2", "xtable", "reshape2",
                   "ggpubr", "plotly", "tidyr", "fbi", "devtools")
new_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)


# Load Packages
#-------------------------------------------------------------------------------------------------
library(devtools)
install_github("manuhuth/hidiTS") #devtools needed to execute this line
library(hidiTS) #our written package containing all functions

library(optimParallel)
library(ggplot2)
library(xtable)
library(reshape2)
library(ggpubr)
library(plotly)
library(tidyr)
library(fbi)


# Prepare Multi Core Work
#-------------------------------------------------------------------------------------------------
cl <- makeCluster(detectCores()-1)
setDefaultCluster(cl = cl)



# 1. -  Data Simulation Difference MSE PCA vs. ML - 500 Samples
#-------------------------------------------------------------------------------------------------

#--------------set_up---------------------------
n_simulation_MSE_vs_PCA <- c(5, 15, 20, 30)
T_simulation_MSE_vs_PCA <- c(7, 15, 20, 30, 40)
number_samples_MSE_vs_PCA <- 500
start_seed_MSE_vs_PCA <- 1
#-----------------------------------------------


    # 1.1 - b = 0.5, Set-up as described in Table 1
    #-------------------------------------------------------------------------------------------------
    #link to function
    df_table_1_1 <- simulation_mse_vs_pca(b=0.5, n_simulation=n_simulation_MSE_vs_PCA, T_simulation=T_simulation_MSE_vs_PCA,
                                          number_samples = number_samples_MSE_vs_PCA, start_seed = start_seed_MSE_vs_PCA)



    # 1.2 - b = 3^0.5, Set-up as described in Table 1
    #-------------------------------------------------------------------------------------------------
    df_table_1_2 <-simulation_mse_vs_pca(b=3^0.5, n_simulation=n_simulation_MSE_vs_PCA, T_simulation=T_simulation_MSE_vs_PCA,
                                         number_samples = number_samples_MSE_vs_PCA, start_seed = start_seed_MSE_vs_PCA)




# 2. -  Data Simulation PCA Simulation Study - 1000 Samples
#-------------------------------------------------------------------------------------------------

#--------------set_up---------------------------
n_simulation_PCA <- c(5,10,seq(10,40,10),seq(50,300,50)) #n smaller 5 is not recommended, since this can lead to errors with the eigenvalue plots
T_simulation_PCA <- c(7,15,30,40,seq(50,300,50))
number_samples_PCA <- 1000
start_seed_PCA <- 1
#-----------------------------------------------

    # 2.1 b = 0.5, Set-up as described in Figure 3
    #-------------------------------------------------------------------------------------------------
    df_pca_b_low <- simulation_pca(b=0.5, T_simulation = T_simulation_PCA, n_simulation = n_simulation_PCA ,
                                   number_samples = number_samples_PCA, start_seed = start_seed_PCA)



    # 2.2 b = 3^0.5, Set-up as described in Figure 3
    #-------------------------------------------------------------------------------------------------
    df_pca_b_medium <- simulation_pca(b=3^0.5, T_simulation = T_simulation_PCA, n_simulation = n_simulation_PCA ,
                                      number_samples = number_samples_PCA, start_seed = start_seed_PCA)



# 3. - Tables and Figures
#-------------------------------------------------------------------------------------------------


    # 3.1 - Table 1 - MSEs PCA vs. ML
    #-------------------------------------------------------------------------------------------------
        df_table_1_2['mse_diff_f'] <- df_table_1_2$mse_pca_f-df_table_1_2$mse_ml_f
        df_table_1_2['mse_diff_f.1'] <- df_table_1_2$mse_pca_f.1-df_table_1_2$mse_ml_f.1
        df_table_1_2['mse_diff_f.2'] <- df_table_1_2$mse_pca_f.2-df_table_1_2$mse_ml_f.2

        sim_data_medium <- df_table_1_2[!duplicated(df_table_1_1[c('n','T')]),]

        ns <- 1:max(n_simulation_PCA)
        Ts <- 1:max(T_simulation_PCA)

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

    #--------------set_up---------------------------
    eigenvalues_to_be_plotted <- 5 #must be smaller than smallest value of n_simulation_PCA in line 63
    ns <- c(5, 20, 50, 100) # must be subset of n_simulation_PCA in line 63
    Ts <- c(7, 15, 30, 100) # must be subset of T_simulation_PCA in line 64
    #-----------------------------------------------

        df <- df_pca_b_low[c("ev1","ev2" ,"ev3","ev4", "ev5", 'n', 'T' )]
        for (index in 1:eigenvalues_to_be_plotted) {
          df[c(paste('ev',index, sep=''))] <- df[c(paste('ev',index, sep=''))] / df['n']
        }

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
        for (index in 1:eigenvalues_to_be_plotted) {
          df[c(paste('ev',index, sep=''))] <- df[c(paste('ev',index, sep=''))] / df['n']
        }

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


    # 3.3 - Figure 7 & Figure 8 - Average Accuracy of the BIC_nT and the BNIC of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------

    #--------------set_up---------------------------
    maxn <- 300
    maxt <- 300
    #-----------------------------------------------

        list_accuracy_plots_low_b <- accuracy_plots(df_pca_b_low, maxn=maxn, maxt=maxt)
        list_accuracy_plots_medium_b <- accuracy_plots(df_pca_b_medium, maxn=maxn, maxt=maxt)
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

    #--------------set_up---------------------------
    ns <-  c(5, 10, 30, 50, 100,300) # must be subset of n_simulation_PCA in line 63
    Ts <-  c(7, 15, 30, 50, 100, 300) # must be subset of T_simulation_PCA in line 64
    #-----------------------------------------------

        sim_data_medium <- df_pca_b_medium[!duplicated(df_pca_b_medium[c('n','T')]),]

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

    #--------------set_up---------------------------
    Ts <-  c(7,15,30, 50, 100,300) # must be subset of T_simulation_PCA in line 64
    #-----------------------------------------------

        list_MSE_plots_low_b <- MSE_plots(df_pca_b_low, Ts = Ts)
        list_MSE_plots_medium_b <- MSE_plots(df_pca_b_medium, Ts = Ts)
        list_MSE_plots_low_b[[1]]
        list_MSE_plots_low_b[[2]]
        list_MSE_plots_low_b[[3]]

        list_MSE_plots_medium_b[[1]]
        list_MSE_plots_medium_b[[2]]
        list_MSE_plots_medium_b[[3]]


    # 3.6 - Table 4 - MSEs of PCA Simulation Data
    #-------------------------------------------------------------------------------------------------

    #--------------set_up---------------------------
    ns <-  c(5, 10, 30, 50, 100,300) # must be subset of n_simulation_PCA in line 63
    Ts <-  c(7, 15, 30, 50, 100, 300) # must be subset of T_simulation_PCA in line 64
    #-----------------------------------------------

        sim_data_medium <- df_pca_b_medium[!duplicated(df_pca_b_medium[c('n','T')]),]
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

    #--------------set_up---------------------------
    Ts <-  c(20) # must be subset of T_simulation_PCA in line 64
    #-----------------------------------------------

        df_Figure_5 <-simulation_mse_vs_pca(b=3^0.5, n_simulation=c(seq(5, 125,10)),
                                            T_simulation=Ts, number_samples = 20, start_seed = 1)
        df_help <- df_Figure_5[df_Figure_5$T==Ts[1],]
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


# 4. - Application
#-------------------------------------------------------------------------------------------------
#rm(list = ls()) #if desired remove all elements from environment


    # 4.1 - Figure 4 a)/b) & Figure 10
    #-------------------------------------------------------------------------------------------------
        prep_data_hidi <- function(data){
          data_hidi <- as.matrix(t(data))
          output <- list("data_fbi"= data, "X"=data_hidi )
          return(output)
        }


        df_paper=fredqd(date_end =as.Date(c("1985/12/01")), date_start=as.Date(c("1960/03/01")),transform = TRUE)

        # remove dates
        df_temp=df_paper[,2:(length(df_paper))]

        #drop Nas
        Na_names <- colnames(df_temp)[colSums(is.na(df_temp)) > 0]
        df_colnames <- colnames(df_temp)
        col_keeps <- setdiff(df_colnames,Na_names)
        df_no_na <- subset(df_temp, select=(col_keeps))

        df_final <- prep_data_hidi(df_no_na)

        X <- df_final$X

        T <- ncol(X)
        n <- nrow(X)

        ic_BIC_n  <- c()
        ic_BIC_T  <- c()
        ic_BIC_nT  <- c()
        ic_BNIC  <- c()

        for (r in 1:n){
          pca_est <- pca_estimator(X, r)
          f_hat <- pca_est$F
          l_hat <- pca_est$Lambda
          x_hat <- l_hat %*% f_hat
          u_hat <- X - x_hat
          RSS <- sum(u_hat^2)

          RSS_part <- log(RSS/n/T)

          BIC_n <- RSS_part + r * log(n)/n
          BIC_T <- RSS_part + r * log(T)/T
          BIC_nT <- RSS_part + r * log(n*T)/n/T * (n+T-r)
          BNIC <- RSS_part + r * log(min(n,T))/n/T * (n+T)

          ic_BIC_n[r]  <- BIC_n
          ic_BIC_T[r]  <-  BIC_T
          ic_BIC_nT[r]  <-  BIC_nT
          ic_BNIC[r]  <- BNIC
          print(paste(' r=',r, sep=''))
        }

        df <- cbind(ic_BIC_n, ic_BIC_T, ic_BIC_nT)[1:200,]

        colnames(df) <- c('BIC_n', 'BIC_T', 'BIC_nT')

        df_long <- melt(df)
        colnames(df_long) <- c('r', 'IC', 'Value')
        (figure1 <- ggplot(df_long, aes(x=r, y=Value, color=IC)) + geom_line() +
            theme(text = element_text(size=28), legend.position="bottom") + labs(color='') )
        ggsave("static/application_BICs_full.png")

        df_long2 <- melt(df[1:20,])
        colnames(df_long2) <- c('r', 'IC', 'Value')
        (figure2 <- ggplot(df_long2, aes(x=r, y=Value, color=IC)) + geom_line() + geom_point() +
            theme(text = element_text(size=28), legend.position="bottom") + labs(color='') )
        ggsave("static/application_BICs_n20.png")

        (figure3 <- ggplot(as.data.frame(cbind('BNIC'=ic_BNIC, 'r'=1:n)[1:200,]), aes(x=r, y=BNIC, lty = 'BNIC')) +
            geom_line()  + theme(text = element_text(size=28), legend.position="bottom") + ylab('Value') + scale_linetype(''))
        ggsave("static/application_Bai_full.png")

        (figure4 <- ggplot(as.data.frame(cbind('BNIC'=ic_BNIC, 'r'=1:n)[1:20,]), aes(x=r, y=BNIC, lty = 'BNIC')) +
            geom_line() + geom_point() + theme(text = element_text(size=28), legend.position="bottom") + ylab('Value')+ scale_linetype(''))
        ggsave("static/application_Bai_n20.png")




    # 4.2 - Figure 4 c)/D) & Figure 11
    #-------------------------------------------------------------------------------------------------
    #--------------set_up---------------------------
    set.seed(1)
    n <- 211
    T <- 104
    q <- 10
    #-----------------------------------------------

        BIC_n_df <- c()
        BIC_T_df <- c()
        BIC_nT_df <- c()
        BNIC_df <- c()

        for (index in 1:1000) {
          data <- sim_data(p = n, T = T, dim_F= q, lags_F=1, lags_X=0, ar_F=1, ar_Y=1, low_X = -3^0.5, up_X = 3^0.5,
                           low_F = 0.3, up_F = 0.6, burn_in = 20, data_only = FALSE, only_stationary = TRUE, vcv_mu = diag(n),
                           adjust_diag = FALSE, geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

          X <- data$X
          ic_BIC_n  <- c()
          ic_BIC_T  <- c()
          ic_BIC_nT  <- c()
          ic_BNIC  <- c()

          for (r in 1:n){
            pca_est <- pca_estimator(X, r)
            f_hat <- pca_est$F
            l_hat <- pca_est$Lambda
            x_hat <- l_hat %*% f_hat
            u_hat <- X - x_hat
            RSS <- sum(u_hat^2)
            RSS_part <- log(RSS/n/T)

            BIC_n <- RSS_part + r * log(n)/n
            BIC_T <- RSS_part + r * log(T)/T
            BIC_nT <- RSS_part + r * log(n*T)/n/T * (n+T-r)
            BNIC <- RSS_part + r * log(min(n,T))/n/T * (n+T)

            ic_BIC_n[r]  <- BIC_n
            ic_BIC_T[r]  <-  BIC_T
            ic_BIC_nT[r]  <-  BIC_nT
            ic_BNIC[r]  <- BNIC
            print(paste('sample=', index, ' r=',r, sep=''))
          }

          BIC_n_df <- rbind(BIC_n_df, ic_BIC_n)
          BIC_T_df <- rbind(BIC_T_df, ic_BIC_T)
          BIC_nT_df <- rbind(BIC_nT_df, ic_BIC_nT)
          BNIC_df <- rbind(BNIC_df, ic_BNIC)
        }

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


        (fig3 <- ggplot(BNIC_df_plot[1:200,], aes(x=r, y=mean,  lty = 'BNIC')) + geom_line() +
            geom_ribbon(aes(x=r, ymin=low, ymax=high),linetype=2, alpha = 0.2) +
            theme(text = element_text(size=28), legend.position="bottom") + scale_linetype('')  + ylab('Value'))


        df_BNIC_n20 <- BNIC_df_plot[1:20,]
        (fig4 <- ggplot(df_BNIC_n20, aes(x=r, y=mean,  lty = 'BNIC')) + geom_line() +
            geom_ribbon(aes(x=r, ymin=low, ymax=high),linetype=2, alpha = 0.2) +
            theme(text = element_text(size=28), legend.position="bottom") + scale_linetype('')  + ylab('Value'))
