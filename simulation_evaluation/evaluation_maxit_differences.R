#load data
load('simulated_data/data_maxit_difference.RData')
library(ggplot2)
library(plotly)


low <- export[export[, 'max_it'] == max_it_low,]
moderate <- export[export[, 'max_it'] == max_it_moderate,]
high <- export[export[, 'max_it'] == max_it_high,]
very_high <- export[export[, 'max_it'] == max_it_very_high,]
very_very_high <- export[export[, 'max_it'] == max_it_very_very_high,]


begin <- 1
last <- 100
plot(low[,'mse_estimated'][begin:last])
points(moderate[,'mse_estimated'][begin:last], col = 'blue', pch = 20)
points(high[,'mse_estimated'][begin:last], col = 'purple', pch = 17)
points(very_high[,'mse_estimated'][begin:last], col = 'orange', pch = 4)
points(very_very_high[,'mse_estimated'][begin:last], col = 'red', pch = 11)



