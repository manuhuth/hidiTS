library(ggplot2)
library("reshape2")
library(ggpubr)

load("simulated_data/raspberry_paper_medium_lambda.RData")

df <- simulated_data[c("ev1","ev2" ,"ev3","ev4", "ev5", 'n', 'T' )]
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
              scale_fill_brewer(palette="Paired") + ggtitle('T=7') + theme(text = element_text(size=20))

t15 <- ggplot(list_df[[2]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
  geom_bar(stat="identity",width = width, position=position_dodge())+
  scale_fill_brewer(palette="Paired") + ggtitle('T=15') + theme(text = element_text(size=20)) + ylab('')

t30 <- ggplot(list_df[[3]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
  geom_bar(stat="identity",width = width, position=position_dodge())+
  scale_fill_brewer(palette="Paired") + ggtitle('T=30') + theme(text = element_text(size=20))
t100 <- ggplot(list_df[[4]], aes(x=n, y=Magnitude, fill=Eigenvalue)) +
  geom_bar(stat="identity",width = width, position=position_dodge())+
  scale_fill_brewer(palette="Paired") + ggtitle('T=100') + theme(text = element_text(size=20)) + ylab('')

figure <- ggarrange(t7,t15,t30,t100, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
