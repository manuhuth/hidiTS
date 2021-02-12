library(xtable)
load("simulated_data/blueberry_paper_small_lambda.RData")
simulated_data['mse_diff_f'] <- simulated_data$mse_pca_f-simulated_data$mse_ml_f
simulated_data['mse_diff_f.1'] <- simulated_data$mse_pca_f.1-simulated_data$mse_ml_f.1
simulated_data['mse_diff_f.2'] <- simulated_data$mse_pca_f.2-simulated_data$mse_ml_f.2

sim_data_medium <- simulated_data[!duplicated(simulated_data[c('n','T')]),]

ns <- 1:400
Ts <- 1:400

data_medium <- sim_data_medium[ which(sim_data_medium$T %in% Ts
                                      & sim_data_medium$n %in% ns), ]
ics_medium <- data_medium[,c('mse_diff_f', 'mse_diff_f.1', 'mse_diff_f.2')]
names(ics_medium) <- c('factor1', 'factor2', 'factor3')



load("simulated_data/blueberry_paper_medium_lambda.RData")
simulated_data['mse_diff_f'] <- simulated_data$mse_pca_f-simulated_data$mse_ml_f
simulated_data['mse_diff_f.1'] <- simulated_data$mse_pca_f.1-simulated_data$mse_ml_f.1
simulated_data['mse_diff_f.2'] <- simulated_data$mse_pca_f.2-simulated_data$mse_ml_f.2

sim_data_small <- simulated_data[!duplicated(simulated_data[c('n','T')]),]

data_small <- sim_data_small[ which(sim_data_small$T %in% Ts
                                    & sim_data_small$n %in% ns), ]
ics_small<- data_small[,c('n', 'T','mse_diff_f', 'mse_diff_f.1', 'mse_diff_f.2')]
names(ics_small) <- c('n', 'T','factor1', 'factor2', 'factor3')

data_final <- cbind(ics_small, ics_medium)

data_final['n'] <- as.integer(data_final[,'n'])
data_final['T'] <- as.integer(data_final[,'T'])

print(xtable(data_final, type = "latex", digits=c(0,0,0,4,4,4, 4,4,4) ), include.rownames=FALSE)
