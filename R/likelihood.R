#rm(list = ls())
################################functions################################


################################theta_wrapper
# maps initial values from column vector into matrix form for optim

theta_wrapper <- function(theta,n,p,q){
  
  
  last_gamma <- p*q*p*q
  gamma <- matrix(data=theta[1:last_gamma],nrow=p*q,ncol = p*q)
  
  
  
  last_lambda <- last_gamma + (n*p*q)
  lambda <- matrix(data=theta[((last_gamma)+1):(last_lambda)],nrow=n,ncol = p*q)
  
  
  #sigma_u
  last_sigma_u <- last_lambda +(n*n)
  sigma_u <- matrix(data=theta[((last_lambda)+1):(last_sigma_u) ],nrow=n,ncol = n)
  
  #sigma_e
  last_sigma_e <- last_sigma_u + (p*q*p*q)
  sigma_e <- matrix(data=theta[((last_sigma_u)+1):(last_sigma_e) ],nrow=p*q,ncol = p*q)
  
  output <- list("gamma" = gamma, "lambda" = lambda, "sigma_u" = sigma_u, "sigma_e" = sigma_e)
  
  return(output) 
}


################################theta_rows
# gives number of rows for theta vector to create initial values for optimization


theta_rows <- function(n,p,q){
  rows=2*(p*q)^2+n*p*q+n^2
  
  return(rows)
  
}


################################initial guess for kalman


sigma_f_init  <- function(p,q,gamma,sigma_e){
  
  temp <- solve(diag((p*q)^2) -kronecker(gamma, gamma, FUN = "*")) %*% vec(sigma_e)
  sigma_f <- matrix(temp, ncol = p*q, nrow=p*q)
  
  return(sigma_f)    
}

f_init <- function(p,q){
  
  f_matrix <- matrix(data = 0, nrow = q*p, ncol = 1)
  
  return (f_matrix)
}


################################Kalman filter functions

post_prediction_sigma_f <- function(sigma_f_prior,lambda,sigma_x_prior){
  
  sigma_f_post <- sigma_f_prior - sigma_f_prior %*% t(lambda) %*% solve(sigma_x_prior) %*% lambda %*% sigma_f_prior 
  
  return(sigma_f_post)
  
}


post_prediction_f <-  function(lambda,f_prior,sigma_f_prior,sigma_x_prior,data,x_prior){
  
  f_post  <- f_prior + sigma_f_prior %*% t(lambda) %*%solve(sigma_x_prior) %*%  (data - x_prior)
  
  return(f_post)
}


x_prediction <- function(lambda,f_prior){
  x_prior <- lambda %*% f_prior
  
  return(x_prior)
}


sigma_x_prediction <- function(lambda,sigma_u,sigma_f_prior){
  sigma_x_prior <- lambda %*% sigma_f_prior %*% t(lambda) + sigma_u 
  return(sigma_x_prior)
}


f_prediciton <- function(gamma,f_post){
  f_prior <- gamma %*% f_post
  
  return(f_prior)
}


sigma_f_prediction <- function(gamma,sigma_f_post,sigma_e){
  sigma_f_prior <- gamma %*% sigma_f_post %*% t(gamma)  + sigma_e
  
  return(sigma_f_prior)
}


################################Kalman wrapper

kalman_wrapper_inital_period <- function(lambda,sigma_u,f_inital,sig_f_inital,data){
  
  #prior observation prediction 
  x_pred <-  x_prediction(lambda = lambda, f_prior = f_inital)
  sig_x_pred <- sigma_x_prediction(lambda=lambda, sigma_u = sigma_u, sigma_f_prior = sig_f_inital)
  # post facor predictiom
  f_post <- post_prediction_f(lambda=lambda,f_prior=f_inital,sigma_f_prior=sig_f_inital,x_prior=x_pred,sigma_x_prior=sig_x_pred,data=data)
  sigma_f_post <- post_prediction_sigma_f(sigma_f_prior=sig_f_inital,sigma_x_prior=sig_x_pred,lambda=lambda)
  
  output <- list("f_post"=f_post,"sigma_f_post"=sigma_f_post)
  
  return(output)
  
}


kalman_wrapper <- function(gamma,lambda,sigma_u,sigma_e,f_post_last,sigma_f_post_last,data){
  
  #prior factor prediction
  f_pred <- f_prediciton(gamma=gamma,f_post = f_post_last)
  sig_f_pred <- sigma_f_prediction(gamma=gamma,sigma_f_post = sigma_f_post_last,sigma_e = sigma_e)
  #prior observation prediction 
  x_pred <-  x_prediction(lambda = lambda, f_prior = f_pred)
  sig_x_pred <- sigma_x_prediction(lambda=lambda, sigma_u = sigma_u, sigma_f_prior = sig_f_pred)
  # post facor predictiom
  f_post <- post_prediction_f(lambda=lambda,f_prior=f_pred,sigma_f_prior=sig_f_pred,x_prior=x_pred,sigma_x_prior=sig_x_pred,data = data)
  sigma_f_post <- post_prediction_sigma_f(sigma_f_prior=sig_f_pred,sigma_x_prior=sig_x_pred,lambda=lambda)
  
  output <- list("f_post"=f_post,"sigma_f_post"=sigma_f_post)
  
  
  return(output)
}


################################likelihood_init
# calcualtes the log density given paramters using the initial guess from kalman filter


likelihood_init <- function(gamma,lambda,sigma_u,sigma_e,data,p,q){
  
  (-n/2 * log(2 * pi) -0.5 * log(det(lambda %*%(sigma_f_init(p=p,q=q,gamma=gamma,sigma_e=sigma_e )) %*% t(lambda) + sigma_u)) 
   - 0.5 * t(((data) -lambda %*% (f_init(q=q,p=p)))) %*% solve(lambda %*%sigma_f_init(p=p,q=q,gamma=gamma,sigma_e=sigma_e)%*% t(lambda) + sigma_u) %*% ((data) - lambda %*%(f_init(q=q,p=p))))
}


################################likelihood
# calcualtes the log density given paramters and kalman posterior f and posterior sigma_f


likelihood <- function(gamma,lambda,sigma_u,sigma_e,f_post,sigma_f_post,data){
  
  (-n/2 * log(2 * pi) -0.5 * log(det(lambda %*%(gamma %*% sigma_f_post %*% t(gamma) +sigma_e) %*% t(lambda) + sigma_u)) 
   - 0.5 * t(((data)-lambda %*% gamma %*% f_post )) %*% solve(lambda %*%( gamma %*% sigma_f_post %*% t(gamma) + sigma_e)%*% t(lambda) + sigma_u) %*% ((data) - lambda %*% gamma %*% f_post))
}


################################wrapper
# function used in optim
# to do: replace loop with apply, check that it produces maximum

wrapper <- function(theta,data,n,p,q,t){
  
  matrices <- theta_wrapper(theta=theta,n=n,p=p,q=q)
  
  
  gamma <- matrices$gamma
  lambda <- matrices$lambda
  sigma_u <- matrices$sigma_u
  sigma_e <- matrices$sigma_e
  
  # initliaze list 
  post_f= vector(mode = "list", length = 0)
  post_sigma_f= vector(mode = "list", length = 0)
  likelihood_list <- vector(mode = "numeric", length = 0)
  
  
  #first period
  likelihood_list[[1]] <- likelihood_init(gamma=gamma,lambda=lambda,sigma_u=sigma_u,sigma_e=sigma_e,data = data[,1],p=p,q=q) 
  
  # save first values
  post_f[[1]] <- kalman_wrapper_inital_period(lambda=lambda,sigma_u=sigma_u,f_inital=f_init(q=q,p=p),sig_f_inital=sigma_f_init(p=p,q=q,gamma=gamma,sigma_e=sigma_e),data = data[,1])$f_post
  post_sigma_f[[1]] <- kalman_wrapper_inital_period(lambda=lambda,sigma_u=sigma_u,f_inital=f_init(q=q,p=p),sig_f_inital=sigma_f_init(p=p,q=q,gamma=gamma,sigma_e=sigma_e),data = data[,1])$sigma_f_post
  
  
  #loop for second period and so on and so on
  
  for(i in 2:t){
    
    likelihood_list[[i]] <- likelihood(gamma=gamma,lambda=lambda,sigma_u=sigma_u,sigma_e=sigma_e,f_post=post_f[[i-1]],sigma_f_post=post_sigma_f[[i-1]],data = data[,i])
    #save second values
    post_f[[i]] <- kalman_wrapper(gamma=gamma,lambda=lambda,sigma_u=sigma_u,sigma_e=sigma_e,f_post_last=post_f[[i-1]],sigma_f_post_last=post_sigma_f[[i-1]],  data = data[,i])$f_post
    post_sigma_f[[i]] <- kalman_wrapper(gamma=gamma,lambda=lambda,sigma_u=sigma_u,sigma_e=sigma_e,f_post_last=post_f[[i-1]],sigma_f_post_last=post_sigma_f[[i-1]],  data = data[,i])$sigma_f_post
  }
  
  output <- sum(likelihood_list)
  
  return(-output)
  
}

###############################short script   ################################

###############################execution   

################################parameter   

#n=3 #number of time series
#q=2 #number of factors
#p=2#factor dimension
#t=5 #lags of xs
#start_ar=4
#rows=theta_rows(n=n,p=p,q=q)

################################data       
#data <-sim_data(p=n,T=t,dim_F=3,lags_F=2,lags_X=1,ar_F = 2, ar_Y = 3)
#data_x <- data$X

####################### define initial inputs for optim 

#init <- matrix(runif(rows, min = 1, max = 10), nrow = rows, ncol = 1)


########## optim function 

#rslt=optim(par=init,wrapper,data=data_x,n=n,p=p,q=q,t=t,method = "SANN")

#rslt$par

#theta_wrapper(rslt$par,n=n,p=p,q=q)
