#Latex table MSES
library(xtable)

load("simulated_data/raspberry_paper_medium_lambda.RData")

sim_data_medium <- simulated_data[!duplicated(simulated_data[c('n','T')]),]

ns <- c(5, 10, 30, 50, 100, 300)
Ts <- c(7, 15, 30, 50, 100, 300)

data_medium <- sim_data_medium[ which(sim_data_medium$T %in% Ts
                                      & sim_data_medium$n %in% ns), ]

ics_medium <- data_medium[,c('mse_pca_f', 'mse_pca_f.1', 'mse_pca_f.2')]
names(ics_medium) <- c('factor1', 'factor2', 'factor3')

load("simulated_data/raspberry_paper_small_lambda.RData")
sim_data_small <- simulated_data[!duplicated(simulated_data[c('n','T')]),]

data_small <- sim_data_small[ which(sim_data_small$T %in% Ts
                                    & sim_data_small$n %in% ns), ]
ics_small<- data_small[,c('n', 'T','mse_pca_f', 'mse_pca_f.1', 'mse_pca_f.2')]
names(ics_small) <- c('n', 'T','factor1', 'factor2', 'factor3')

data_final <- cbind(ics_small, ics_medium)

data_final['n'] <- as.integer(data_final[,'n'])
data_final['T'] <- as.integer(data_final[,'T'])

print(xtable(data_final, type = "latex", digits=c(0,0,0,4,4,4, 4,4,4) ), include.rownames=FALSE)
