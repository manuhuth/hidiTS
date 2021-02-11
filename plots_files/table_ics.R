#Latex table MSES
library(xtable)

load("simulated_data/raspberry_paper_medium_lambda.RData")
sim_data_medium <- simulated_data[!duplicated(simulated_data[c('n','T')]),]

ns <- c(5, 10, 30, 50, 100,300)
Ts <- c(7, 15, 30, 50, 100, 300)

data_medium <- sim_data_medium[ which(sim_data_medium$T %in% Ts
                         & sim_data_medium$n %in% ns), ]

ics_medium <- data_medium[,c('pca_bai_right', 'pca_bic_right', 'pca_bic_T_right', 'pca_bic_nT_right')]
names(ics_medium) <- c('BNIC_med', 'BIC_n_med', 'BIC_T_med', 'BIC_nT_med')

load("simulated_data/raspberry_paper_small_lambda.RData")
sim_data_small <- simulated_data[!duplicated(simulated_data[c('n','T')]),]

data_small <- sim_data_small[ which(sim_data_small$T %in% Ts
                        & sim_data_small$n %in% ns), ]
ics_small<- data_small[,c('n', 'T','pca_bai_right', 'pca_bic_right', 'pca_bic_T_right', 'pca_bic_nT_right')]
names(ics_small) <- c('n', 'T','BNIC_small', 'BIC_n_small', 'BIC_T_small', 'BIC_nT_small')

data_final <- cbind(ics_small, ics_medium)

data_final['n'] <- as.integer(data_final[,'n'])
data_final['T'] <- as.integer(data_final[,'T'])

print(xtable(data_final, type = "latex", digits=c(0,0,0,4,4,4,4,4, 4,4,4)), include.rownames=FALSE)
