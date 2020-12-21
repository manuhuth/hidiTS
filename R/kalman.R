############ Kalman ####
rm(list = ls())
library(hidiTS)
library(matrixcalc)
set.seed(42)


######## dgp von manu ######
q=3 #number of dynamic factors
p=3 # number of lags of the factors
T=34 # time when the AR part starts =34-4=30
n=6 # number of observed variables
lags_obs_var=2 # number of static factors
start_ar=4
gamma <- random_matrices(n = 1, row = 9, col = 9) #n= 1 einfachshalber
cov_stationarity(gamma,epsilon = 0.05)
F <- sim_F(dim = 3, lags = p, T = T)$F
Stack<- stack_F(F = F, lags = 3, start_ar=start_ar)
X= sim_X(F = F, p = n, lags = lags_obs_var, start_ar = start_ar)

lambda<- random_matrices(n=1, row = n, col = q*p)
#Sigma_eta=matrix(rep(diag(q),p), nrow = q, ncol = p*q)
#Sigma_e=list(matrix(data = rbind(Sigma_eta,matrix(data=0,nrow = p*q-q, ncol = p*q)), nrow = p*q, ncol=p*q))
Sigma_e=list(diag(p*q))
Sigma_u= list(diag(n))



Kalman= function(q,p,T,n, lambda, gamma, Sigma_e, Sigma_u ){ # plug in all  input variables

  ####### set the lists #########
  prior_F= vector(mode = "list", length = T-start_ar)
  prior_P= vector(mode = "list", length = T-start_ar)
  prior_x= vector(mode = "list", length = T-start_ar)
  prior_V= vector(mode = "list", length = T-start_ar)
  post_F= vector(mode = "list", length = T-start_ar)
  post_P= vector(mode = "list", length = T-start_ar)

  ########## Set the first values ###########
  prior_F[[1]]= matrix(data = 0, nrow = q*p, ncol = 1)
  prior_x[[1]]= lambda[[1]]%*%prior_F[[1]]
  vec_prior_P= solve(diag((p*q)^2) -kronecker(gamma[[1]], gamma[[1]], FUN = "*")) %*% vec(Sigma_e[[1]])
  prior_P[[1]]= matrix(vec_prior_P, ncol = p*q, nrow=p*q)
  prior_V[[1]]=(lambda[[1]]%*%prior_P[[1]]%*%t(lambda[[1]])) + Sigma_u[[1]]

  ######## Updating first values before starting the loop ##########
  post_F[[1]]= prior_F[[1]] +prior_P[[1]]%*%t(lambda[[1]])%*% solve(prior_V[[1]])%*% (X[[1]][,1] - prior_x[[1]]) #V anschauen
  post_P[[1]]= prior_P[[1]] - prior_P[[1]]%*%t(lambda[[1]])%*% solve(prior_V[[1]]) %*% lambda[[1]]%*%prior_P[[1]] #nochmal checken

  Kalman.predF <- function(i){
    a <- i-1
    prior_F[[i]] <<-gamma[[1]]%*% post_F[[i-1]]
    return(prior_F)

  }

  Kalman.predP<-function(i){
    prior_P[[i]] <<- gamma[[1]]%*% prior_P[[i-1]]%*%t(gamma[[1]]) + Sigma_e[[1]]
    return(prior_P)
  }

  Kalman.obsx <-function(i){
    prior_x[[i]] <<- lambda[[1]]%*%prior_F[[i]]
    return(prior_x)

  }
  Kalman.obsV <- function(i){
    prior_V[[i]] <<- (lambda[[1]]%*%prior_P[[i]]%*%t(lambda[[1]])) + Sigma_u[[1]]
    return(prior_V[[i]])
  }
  Kalman.updF <- function(i){
    post_F[[i]] <<- prior_F[[i]] +prior_P[[i]]%*%t(lambda[[1]])%*% solve(prior_V[[i]])%*% (X[[1]][,i] - prior_x[[i]])
    return(post_F[[i]])
  }
  Kalman.updP <- function(i){
    post_P[[i]] <<- prior_P[[i]] - prior_P[[i]]%*%t(lambda[[1]])%*% solve(prior_V[[i]]) %*% lambda[[1]]%*%prior_P[[i]]
    return(post_P[[i]])
  }

  Looping <- function(i){
    Pred1 <-Kalman.predF(i)
    Pred2 <-Kalman.predP(i)

    Obs1 <- Kalman.obsx(i)
    Obs2 <- Kalman.obsV(i)

    Upd1 <- Kalman.updF(i)
    Upd2 <- Kalman.updP(i)
    return(post_F)
  }
  res = lapply(2:(T-start_ar ), Looping)
  prior_F <<- prior_F
  prior_P <<- prior_P
  prior_x <<- prior_x
  prior_V <<- prior_V
  post_F <<- post_F
  post_P <<- post_P
  #return(list(prior_F,prior_P, prior_x, prior_V,post_F, post_P))

}

Result=Kalman(q=q, p=p, T=T, n=n, lambda = lambda, gamma = gamma, Sigma_e = Sigma_e, Sigma_u = Sigma_u)
