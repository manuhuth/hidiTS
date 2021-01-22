#devtools::install(dependencies = TRUE)
library(hidiTS)
rm(list = ls())
set.seed(123)

n=10
p=0
q=2
t=15
k=1

data_test <- sim_data(p = n, T = t, dim_F= q, lags_F=k, lags_X=p, ar_F=1, ar_Y=1, A = NULL, low_X = -3^0.5, up_X = 3^0.5,
                      L = NULL, low_F = 0.2, up_F = 0.9,
                      beta = NULL, gamma = NULL, low_Y = -0.9, up_Y = 0.9,
                      burn_in = 20, data_only = FALSE,
                      only_stationary = TRUE, epsilon = 0.0002, max_it_station = 500, adjust_diag = FALSE,
                      geometric_F =TRUE, diag_F = TRUE, geometric_X =FALSE, geometric_Y =FALSE)

data_param_init <- starting_values_ML(data_test,sigma_u_diag=TRUE, sigma_u_ID=TRUE, sigma_eta_ID=TRUE)

matrices <- matrix_form(data=data_param_init,n=n,p=p,q=q,k=k,gamma_res=FALSE,lambda_res=TRUE,sigma_u_diag=TRUE)

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


rslt1 <- optim(par=data_param_init,likelihood_wrapper,data_x=data_test,n=n,
               p=p,q=q,k=k,t=t,gamma_res=FALSE,lambda_res=TRUE,sigma_u_diag=TRUE,post_F=fsmooth1,
               post_P=psmooth1,control=list(fnscale=-1, trace=list(1)),method = "Nelder-Mead")


matrices_2 <- matrix_form(rslt1$par,n=n,p=p,q=q,k=k,gamma_res=FALSE,lambda_res=TRUE,sigma_u_diag=TRUE)

lambda2 <- list(matrices_2$lambda)
gamma2 <- list(matrices_2$gamma)
sigma_e2 <- list(matrices_2$sigma_e)
sigma_u2 <- list(matrices_2$sigma_u)

Kalman_second <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda2, gamma=gamma2,
                        Sigma_e=sigma_e2, Sigma_u=sigma_u2 ,start_ar=0,X=data_test)

fsmooth2 <- Kalman_second$Fsmooth
psmooth2 <- Kalman_second$Psmooth

conv <- convert_factors_ML(Lambda=lambda2, factors=fsmooth2,q=q)

con_dgp <- convert_factors_dgp(data_test)
plot(conv$F[2,], type = 'l')
lines(con_dgp$F[2,], type='l', col='red')


#rslt2=optim(par=rslt1$par,likelihood_wrapper,data_x=data_test,n=n,p=p,q=q,k=k,t=t,post_F=fsmooth2,post_P=psmooth2,control=list(fnscale=-1),method = "Nelder-Mead")
#matrices_3 <- matrix_form(rslt2$par,n=n,p=p,q=q,k=k)

#lambda3 <- list(matrices_3$lambda)

#gamma3 <- list(matrices_3$gamma)

#sigma_e3 <- list(matrices_3$sigma_e)

#sigma_u3 <- list(matrices_3$sigma_u)



#Kalman_third <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda3, gamma=gamma3, Sigma_e=sigma_e3, Sigma_u=sigma_u3 ,start_ar=0,X=data_test)

#fsmooth3 <- Kalman_third$Fsmooth
#psmooth3 <- Kalman_third$Psmooth
