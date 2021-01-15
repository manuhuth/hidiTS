#rm(list = ls())

#install.packages("devtools") #only once


#library(hidiTS)
#library(matrixcalc)
#set.seed(1234)
set.seed(4511)
#set.seed(98531)
#set.seed(445)


n=20
p=3
q=4
t=10
k=3

#gamma_matrix(p=p,q=q,k=k,data=runif(q,1,2),restricted=TRUE)

elements <- number_of_param(n=n,p=p,q=q,k=k,lambda_res = TRUE,gamma_res = TRUE,sigma_u_diag=FALSE)

data_param_init <- runif(elements,0.05,0.2)

data_test <- sim_data(p=n,T=t)

matrices <- matrix_form(data=data_param_init,n=n,p=p,q=q,k=k,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=FALSE)

lambda <- list(matrices$lambda)

gamma <- list(matrices$gamma)

sigma_e <- list(matrices$sigma_e)

sigma_u <- list(matrices$sigma_u)  

Kalman_first <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda, gamma=gamma, Sigma_e=sigma_e, Sigma_u=sigma_u ,start_ar=0,X=data_test)


fsmooth1 <- Kalman_first$Fsmooth
psmooth1 <- Kalman_first$Psmooth

list_param= vector(mode = "list", length = 0)
list_f= vector(mode = "list", length = 0)
list_sigma_f= vector(mode = "list", length = 0)

list_param[[1]] <- data_param_init
list_f[[1]] <- fsmooth1
list_sigma_f[[1]] <- psmooth1


rslt1 <- optim(par=data_param_init,likelihood_wrapper,data_x=data_test,n=n,p=p,q=q,k=k,t=t,post_F=fsmooth1,post_P=psmooth1,control=list(fnscale=-1),method = "Nelder-Mead")


matrices_2 <- matrix_form(rslt1$par,n=n,p=p,q=q,k=k)

lambda2 <- list(matrices_2$lambda)

gamma2 <- list(matrices_2$gamma)

sigma_e2 <- list(matrices_2$sigma_e)

sigma_u2 <- list(matrices_2$sigma_u)


Kalman_second <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda2, gamma=gamma2, Sigma_e=sigma_e2, Sigma_u=sigma_u2 ,start_ar=0,X=data_test)

fsmooth2 <- Kalman_second$Fsmooth
psmooth2 <- Kalman_second$Psmooth


do.call(cbind,fsmooth2)
do.call(cbind,fsmooth2)[1:q,]


#rslt2=optim(par=rslt1$par,likelihood_wrapper,data_x=data_test,n=n,p=p,q=q,k=k,t=t,post_F=fsmooth2,post_P=psmooth2,control=list(fnscale=-1),method = "Nelder-Mead")
#matrices_3 <- matrix_form(rslt2$par,n=n,p=p,q=q,k=k)

#lambda3 <- list(matrices_3$lambda)

#gamma3 <- list(matrices_3$gamma)

#sigma_e3 <- list(matrices_3$sigma_e)

#sigma_u3 <- list(matrices_3$sigma_u)



#Kalman_third <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda3, gamma=gamma3, Sigma_e=sigma_e3, Sigma_u=sigma_u3 ,start_ar=0,X=data_test)

#fsmooth3 <- Kalman_third$Fsmooth
#psmooth3 <- Kalman_third$Psmooth
