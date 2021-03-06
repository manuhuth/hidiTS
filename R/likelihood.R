#'auxiliary function in gamma_matrix that creates diag small gamma matrices
#'@param chuncks list of equially sized data clusters
#'@return list of gamma matrices
create_gammas <- function(chunks,q){

  gammas <- vector(mode = "list", length = length(chunks))
  for (i in 1:length(chunks)){

    gammas[[i]] <- diag(x=chunks[[i]],q,q)
    #print(gammas[[i]])

  }

  return(gammas)
}

#' auxiliary function in gamma_matrix that creates zero matrices
#' @param p number of lags for factors
#' @param q number of factors per period
#' @param k order of gamma lag polynomial/number of gamma matrices
#' @return zeros matrices
create_zeros <- function(p,q,k){

  number_zeros <- (p+1)-k
  if(number_zeros <= 0){
    zeros <- vector(mode = "list", length = 0)
  } else {
    zeros <- vector(mode = "list", length = number_zeros)
    for (i in 1:number_zeros){

      zeros[[i]] <- diag(x=0,q,q)
      #print(gammas[[i]])
    }
  }

  return(zeros)
}

#' creates Gamma matrix as  described in paper
#' @param p number of lags for factors
#' @param q number of factors per period
#' @param k order of gamma lag polynomial
#' @param data vector of data points to supply
#' @param restricted if TRUE restricts small gammas to lag structure
#' @return Gamma matrix
gamma_matrix <- function(p,q,k,data,restricted=TRUE){

  if(isTRUE(restricted)){

    gammas <- vector(mode = "list", length = k)

    for(i in 1:k){
      gammas[[i]] <- diag(x=data^i,q)
    }


  } else {
    chunks <- split(data, ceiling(seq_along(data)/q))

    gammas <- create_gammas(chunks = chunks)


  }

  zeros <- create_zeros(p=p,q=q,k=k)

  matrices <- c(gammas[1:(k)],zeros)

  upper <- do.call(cbind, matrices)

  m <- length(matrices)
  n <- ncol(matrices[[1]])
  lower_diag <- diag((m - 1) * q) # create lower block matrix

  zero_right <- matrix(rep(0, len = (m - 1) * q^2), ncol = n, nrow = (m - 1) * n)
  lower <- cbind(lower_diag, zero_right)

  output <- rbind(upper, lower)

  return(output)

}

#' creates Lambda matrix as  described in paper
#' @param n number of observable x/number of time series
#' @param p number of lags for factors
#' @param q number of factors per period
#' @param data vector of data points to supply
#' @param restricted if TRUE restricts Lambda to structure by Bai,Ng (2013)
#' @return Gamma matrix
lambda_matrix <-function(data,n,p,q,restricted=TRUE){

  if( (n < (p+1)*q) & (restricted==TRUE) ){

    stop('Restriction of lambda not possible given n,p,q')

  }

  if(isTRUE(restricted)){


    lambda <- matrix(data = 0,nrow=n,ncol=(p+1)*q)

    vectors <- vector(mode ="list",length = 0)
    data_temp <- data

    for (i in 0:((p+1)*q-1)){

      vectors[[(i+1)]] <- data_temp[1:(n-i)]
      data_temp <- tail(data_temp, -(length(vectors[[(i+1)]])))

    }

    for (i in 1:((p+1)*q)){
      lambda[(i:n),i]=vectors[[i]]
    }



    return(lambda)


  }else{

    lambda <- matrix(data=data, nrow = n, ncol = (p+1)*q)

    return(lambda)

  }
}


#' creates symmetric matrix of dimension n x n
#' @param data vector of data points to supply
#' @param n dimension
#' @return symmetric matrix
sigma_u_matrix <- function(data,n,diag_res=FALSE){

  if(isTRUE(diag_res)){
    sigma_u <-diag(x=data,n,n)
    return(sigma_u)
  } else  {

    sigma_u <- matrix(rep(0,n*n),n,n)
    sigma_u[lower.tri(sigma_u,diag=TRUE)] <- data
    sigma_u <- sigma_u+t(sigma_u)

    diag(sigma_u)<-diag(sigma_u)/2

    return(sigma_u)
  }

}

#' creates sigma_e matrix as described in paper
#' @param p number of lags for factors
#' @param q number of factors per period
#' @param data vector of data points to supply
#' @return symmetric matrix
sigma_e_matrix <- function(p,q,data){
  #create 0 matrix
  sigma_e <- matrix(data=0, nrow = (p+1)*q, ncol = (p+1)*q)
  #create sigma_eta
  diag_m <- diag(data,q,q)
  #place sigma_eta
  sigma_e[1:q,1:q] <- diag_m

  return(sigma_e)

}




# maps initial values from column vector into matrix form for optim
#' wraps data points in matrices in theta
#' @param data vector of data points to supply
#' @param n number of observable x/number of time series
#' @param p number of lags for factors
#' @param q number of factors per period
#' @param k order of gamma lag polynomial
#' @param gamma_res if TRUE restricts small gammas to lag structure
#' @param lambda_res if TRUE restricts Lambda to structure by Bai,Ng (2013)
#' @return number of paramaters over all matrices in theta
matrix_form <- function(data,n,p,q,k,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=FALSE){

  op=options(warn=2)


  if(isTRUE(gamma_res)){
    last_gamma <- q
    gamma <- gamma_matrix(p=p,q=q,k=k,data = data[1:q],restricted=TRUE)

  }  else {
    last_gamma <- k*q
    gamma <- gamma_matrix(p=p,q=q,k=k,data = data[1:(k*q)],restricted=FALSE)

  }



  if(isTRUE(lambda_res)){

    if( (n < (p+1)*q) ){stop('Restriction of lambda not possible given n,p,q')}

    last_lambda <- last_gamma + 0.5*(1+(p+1)*q)*((p+1)*q) + ((n-(p+1)*q)*(p+1)*q)
    lambda <- lambda_matrix(data=data[(last_gamma+1):last_lambda],n=n,p=p,q=q,restricted = TRUE)

  }  else {

    last_lambda <- last_gamma + n*(p+1)*q
    lambda <- lambda_matrix(data=data[(last_gamma+1):last_lambda],n=n,p=p,q=q,restricted = FALSE)

  }



  if(sigma_u_diag==FALSE){
    #sigma_u
    last_sigma_u <- last_lambda + 0.5*(1+n)*n
    sigma_u <- sigma_u_matrix(data=data[((last_lambda)+1):(last_sigma_u) ],n=n,diag_res=FALSE)

  } else {
    last_sigma_u <- last_lambda + n
    sigma_u <- sigma_u_matrix(data=data[((last_lambda)+1):(last_sigma_u) ],n=n,diag_res=TRUE)

  }
  #sigma_e


  sigma_e <- sigma_e_matrix(data=data[(last_sigma_u+1):(length(data))],p=p,q=q)

  output <- list("gamma" = gamma, "lambda" = lambda, "sigma_u" = sigma_u, "sigma_e" = sigma_e)

  options(op)
  return(output)

}



# gives number of rows for theta vector to create initial values for optimization
number_of_param <- function(n,p,q,k,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=FALSE){

  F_dim <- (p+1)*q

  #lambda paramaters

  #restricted lambda
  if(isTRUE(lambda_res)) {

    if( (n < (p+1)*q) ){stop('Restriction of lambda not possible given n,p,q')}

    p_lambda <- 0.5*(1+F_dim)*F_dim + (n-F_dim)*F_dim

  } else {
    #unrestricted lambda
    p_lambda <- n*F_dim
  }
  #gamma paramaters

  #restricted gamma
  if(isTRUE(gamma_res)) {

    p_gamma <- q

  } else {

    #unrestricted gamma
    p_gamma <- q*k
  }

  #sigma_u paramaters
  if(sigma_u_diag==TRUE) {

    p_sigma_u <- n

  } else {

    p_sigma_u <- 0.5*(1+n)*n

  }


  #sigma_e paramters
  p_sigma_e <- q


  num_param =  p_lambda + p_gamma + p_sigma_u + p_sigma_e



  return(num_param)(n,p,q,k,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=FALSE)

}



#initial guess for kalman sigma
sigma_f_init  <- function(p,q,gamma,sigma_e){

  temp <- solve(diag(((p+1)*q)^2) - kronecker(gamma, gamma, FUN = "*")) %*% vec(sigma_e)
  sigma_f <- matrix(temp, ncol = (p+1)*q, nrow=(p+1)*q)

  return(sigma_f)
}

f_init <- function(p,q){

  f_matrix <- matrix(data = 0, nrow = (p+1)*q, ncol = 1)

  return (f_matrix)
}



# calcualtes the log density given paramters using the initial guess from kalman filter
likelihood_init <- function(gamma,lambda,sigma_u,sigma_e,data,p,q){
  d <- det(lambda %*%(sigma_f_init(p=p,q=q,gamma=gamma,sigma_e=sigma_e )) %*% t(lambda) + sigma_u)

  return(-n/2 * log(2 * pi) -0.5 * log(d)
   - 0.5 * t(((data) -lambda %*% (f_init(q=q,p=p)))) %*% solve(lambda %*%sigma_f_init(p=p,q=q,gamma=gamma,sigma_e=sigma_e)%*% t(lambda) + sigma_u) %*% ((data) - lambda %*%(f_init(q=q,p=p))))
}



# calcualtes the log density given paramters and kalman posterior f and posterior sigma_f
likelihood <- function(gamma,lambda,sigma_u,sigma_e,f_post,sigma_f_post,data){

  (-n/2 * log(2 * pi) -0.5 * log(det(lambda %*%(gamma %*% sigma_f_post %*% t(gamma) +sigma_e) %*% t(lambda) + sigma_u))
   - 0.5 * t(((data)-lambda %*% gamma %*% f_post )) %*% solve(lambda %*%( gamma %*% sigma_f_post %*% t(gamma) + sigma_e)%*% t(lambda) + sigma_u) %*% ((data) - lambda %*% gamma %*% f_post))
}


#' calculates the likelihood as  specified in paper
#' @param data_param vector of data points to supply
#' @param data_x matrix data of the observable time series
#' @param n number of observable x/number of time series
#' @param p number of lags for factors
#' @param q number of factors per period
#' @param k order of gamma lag polynomial
#' @param t time dimension of observable x/time series
#' @param gamma_res if TRUE restricts small gammas to lag structure
#' @param lambda_res if TRUE restricts Lambda to structure by Bai,Ng (2013)
#' @param post_F list of post_F estimations by Kalman Filter
#' @param post_P list of post_P/sigma_F estimations by Kalman Filter
#' @return sum of likelihood values across time
likelihood_wrapper <- function(data_param,data_x,n,p,q,k,t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=FALSE,post_F,post_P){

  matrices <- matrix_form(data=data_param,n=n,p=p,q=q,k=k,gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag)


  gamma <- matrices$gamma
  lambda <- matrices$lambda
  sigma_u <- matrices$sigma_u
  sigma_e <- matrices$sigma_e

  likelihood_list <- vector(mode = "numeric", length = 0)
  #first period
  likelihood_list[[1]] <- likelihood_init(gamma=gamma,lambda=lambda,sigma_u=sigma_u,sigma_e=sigma_e,data = data_x[[2]][,1],p=p,q=q)

  #loop for second period and so on and so on

  for(i in 2:t){

    likelihood_list[[i]] <- likelihood(gamma=gamma,lambda=lambda,sigma_u=sigma_u,sigma_e=sigma_e,f_post=post_F[[i-1]],sigma_f_post=post_P[[i-1]],data = data_x[[2]][,i])

  }

  output <- sum(likelihood_list)

  return(output)

}

optim_wrapper <- function(data_param,optim_func,data_x,n,p,q,k,t,gamma_res=FALSE,lambda_res=TRUE,sigma_u_diag=TRUE,post_F,post_P,
                          lower,upper,method = "L-BFGS-B", max_it, parallel=FALSE, trace=0, forward = FALSE, loginfo=FALSE){

  if (isTRUE(parallel)) {
    rslt=optimParallel(par=data_param, fn=optim_func,data_x=data_x,n=n,p=p,q=q,k=k,t=t,gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag,post_F=post_F,post_P=post_P,
               control=list(fnscale=-1, maxit=max_it, trace=trace),method = method, lower=lower, upper=upper, parallel = list(forward = forward, loginfo=loginfo))
  } else{
    rslt=optim(par=data_param, fn=optim_func,data_x=data_x,n=n,p=p,q=q,k=k,t=t,gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag,post_F=post_F,post_P=post_P,
               control=list(fnscale=-1,trace=trace, maxit=max_it),method = method, lower=lower, upper=upper)
  }
  #rslt$par

  params <- rslt$par

  value <- rslt$value

  matrices <- matrix_form(rslt$par,n=n,p=p,q=q,k=k, gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag)

  gamma <- list(matrices$gamma)

  lambda <- list(matrices$lambda)

  sigma_e <- list(matrices$sigma_e)

  sigma_u <- list(matrices$sigma_u)

  Kalman_rslt <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda, gamma=gamma, Sigma_e=sigma_e, Sigma_u=sigma_u ,start_ar=0,X=data_x)

  Fsmooth_upd <- Kalman_rslt$Fsmooth
  Psmooth_upd <- Kalman_rslt$Psmooth


  output <- list("params" = params, "post_F" = Fsmooth_upd, "post_P" =Psmooth_upd, "value"=value)
  return(output)



}

#set cluster -> must be outside
#cl <- makeCluster(detectCores()-1)
#clusterExport(cl, 'n')
#setDefaultCluster(cl = cl)

estimate_f <- function(data_x,n,p,q,k,t,gamma_res=TRUE,lambda_res=TRUE,sigma_u_diag=FALSE,it=4,method = "L-BFGS-B", parallel = FALSE,
                       max_it = 50, trace=0, forward = FALSE, loginfo = FALSE){

  start_object <- starting_values_ML(t=t,q=q,n=n,k=k,p=p,data_x,sigma_u_diag=sigma_u_diag, sigma_u_ID=TRUE, sigma_eta_ID=TRUE)
  data_param_init <- start_object$data
  data_param_names <- start_object$names
  matrices <- matrix_form(data=data_param_init,n=n,p=p,q=q,k=k,gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag)

  lambda <- list(matrices$lambda)

  gamma <- list(matrices$gamma)

  sigma_e <- list(matrices$sigma_e)

  sigma_u <- list(matrices$sigma_u)

  Kalman_first <- Kalman(q=q,p=p,T=t,n=n, lambda=lambda, gamma=gamma, Sigma_e=sigma_e, Sigma_u=sigma_u ,start_ar=0,X=data_x)

  fsmooth1 <- Kalman_first$Fsmooth
  psmooth1 <- Kalman_first$Psmooth

  list_param= vector(mode = "list", length = 0)
  list_f= vector(mode = "list", length = 0)
  list_sigma_f= vector(mode = "list", length = 0)
  value=vector(mode = "list", length = 0)

  list_param[[1]] <- data_param_init
  list_f[[1]] <- fsmooth1
  list_sigma_f[[1]] <- psmooth1

  lower_gamma <- -0.8
  upper_gamma <- 0.8
  lower_lambda <- -Inf
  upper_lambda <- Inf
  lower_sigma_u <- 0.01
  upper_sigma_u <- Inf
  lower_sigma_e <- 0.01
  upper_sigma_e <- Inf

  names <- c('gamma', 'lambda', 'sigma_u', 'sigma_e')
  lowers <- c('gamma'=lower_gamma, 'lambda'=lower_lambda, 'sigma_u'=lower_sigma_u, 'sigma_e'=lower_sigma_e )
  uppers <- c('gamma'=upper_gamma, 'lambda'=upper_lambda, 'sigma_u'=upper_sigma_u, 'sigma_e'=upper_sigma_e )

  lower <- rep(NaN, length(data_param_init))
  upper <- rep(NaN, length(data_param_init))

  for (index in names){
    lower <- replace(lower, data_param_names == index, lowers[index])
    upper <- replace(upper, data_param_names == index, uppers[index])
  }

  # perform loop with likelihood estimation  for parameters and kalman filter

  counter <- 1

  while (counter <= it){

    rslt_while_counter <- optim_wrapper(data_param=list_param[[(length(list_param))]],optim_func=likelihood_wrapper,data_x=data_x,n=n,
                                        p=p,q=q,k=k,t=t,gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag,
                                        post_F=list_f[[(length(list_f))]],post_P=list_sigma_f[[(length(list_sigma_f))]],
                                        method = method, lower=lower, upper=upper, max_it=max_it, parallel = parallel,
                                        trace=trace, forward = forward, loginfo=loginfo)

    list_param <- append(list_param,list(rslt_while_counter$params))
    list_f[[counter+1]] <- rslt_while_counter$post_F
    list_sigma_f <- append(list_sigma_f, list(rslt_while_counter$post_P))
    value <- rslt_while_counter$value

    counter <- counter + 1

  }


  last_param <- list_param[[length(list_param)]]

  last_matrices<- matrix_form(data=last_param,n=n,p=p,q=q,k=k,gamma_res=gamma_res,lambda_res=lambda_res,sigma_u_diag=sigma_u_diag)

  last_lambda <- list(last_matrices$lambda)

  last_gamma <- list(last_matrices$gamma)

  last_sigma_e <- list(last_matrices$sigma_e)

  last_sigma_u <- list(last_matrices$sigma_u)


  # map results from kalman f into f vector

  last_fs <- list_f[[(length(list_f))]]


  f_final <- do.call(cbind,last_fs)[1:q,]
  F_final <- do.call(cbind,last_fs)


  output =list("f_final"=f_final,"F_final"=F_final,"list_f"=list_f,"param" =list_param, "lambda"=last_lambda,
               "gamma" =last_gamma, "sigma_e"=last_sigma_e,"sigma_u" =last_sigma_u, "value"=value)

  return(output)


}


starting_values_ML <- function(t,q,n,k,p,data_test, sigma_u_diag=TRUE, sigma_u_ID=TRUE, sigma_eta_ID=TRUE) {
  #input must be data list created with data_only = FALSE

 # t <- ncol(data_test$F)
 # q <- nrow(data_test$F)
 # n <- nrow(data_test$X)
 # k <- length(data_test$A)
 # p <- length(data_test$L) - 1

  pca_est <- pca_estimator(data_test$X, (p+1)*q)
  converted_pca_est <- convert_pca_estimator(pca_est, q, big_F = TRUE)

  Lambda_start <- converted_pca_est$Lambda
  F_start <- converted_pca_est$F


  X_hat <- Lambda_start %*% F_start
  X <- data_test$X
  u_hat <- X - X_hat
  u_VCV <- var(t(u_hat))
  if (isTRUE(sigma_u_ID)) {
    u_VCV <- diag(1,nrow=nrow(u_VCV), ncol=ncol(u_VCV))
  }

  #make regression with small fs -> build small fs and do gamma regression separate
  F_estimated <- c()
  gammas <- c()
  for (fact in 1:q) {

    F_lags <- F_start[fact, 1:(t)]
    for (index in 1:k) {
      F_lags <- cbind(F_lags, c(rep(0,index),F_start[fact, 1:(t-index)]) )
    }
    F_Y <- F_lags[(k+1):t,1]
    F_X <- F_lags[(k+1):t,2:(k+1)]

    Gamma <- solve(t(F_X) %*% F_X) %*% t(F_X) %*% F_Y

    #if (Gamma > 1)

    gammas <- rbind(gammas, t(Gamma))
    F_hat <- F_X %*% Gamma

    F_estimated <- rbind(F_estimated, t(F_hat))
  }

  gammas[gammas > 0.9] <- 0.9

  eta_hat <- F_start[1:q,1:(t-k)] - F_estimated
  sigma_e_hat <- var(t(eta_hat))

  if (isTRUE(sigma_eta_ID)) {
    sigma_e_hat <- diag(1,nrow=nrow(sigma_e_hat), ncol=ncol(sigma_e_hat))
  }

  data_param_init <- vec(gammas)
  params_names <- rep('gamma', length(vec(gammas)))

  for (index in 1:((p+1)*q)) {
    data_param_init <- c(data_param_init, Lambda_start[(1+index-1):n,index])
    params_names <- c(params_names, rep('lambda', length(Lambda_start[(1+index-1):n,index])))
  }


  if (isTRUE(sigma_u_diag)) {
    data_param_init <- c(data_param_init, diag(u_VCV))
    params_names <- c(params_names, rep('sigma_u', length(diag(u_VCV))) )
  } else{
    for (index in 1:(ncol(u_VCV)) ) {
      data_param_init <- c(data_param_init, u_VCV[(1+index-1):n,index])
      params_names <- c(params_names, rep('sigma_u', length(u_VCV[(1+index-1):n,index])) )
    }
  }

  data_param_init <- c(data_param_init, diag(sigma_e_hat))
  params_names <- c(params_names, rep('sigma_e', length(diag(sigma_e_hat))) )

  return(list("data"=data_param_init, "names"=params_names))
}
