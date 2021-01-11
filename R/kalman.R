Kalman= function(q,p,T,n, lambda, gamma, Sigma_e, Sigma_u ,start_ar,X){

  prior_F= vector(mode = "list", length = T-start_ar)
  prior_P= vector(mode = "list", length = T-start_ar)
  prior_x= vector(mode = "list", length = T-start_ar)
  prior_V= vector(mode = "list", length = T-start_ar)
  post_F= vector(mode = "list", length = T-start_ar)
  post_P= vector(mode = "list", length = T-start_ar)

  prior_F[[1]]= matrix(data = 0, nrow = q*(p+1), ncol = 1)
  prior_x[[1]]= lambda[[1]]%*%prior_F[[1]]
  vec_prior_P= solve(diag(((p+1)*q)^2) -kronecker(gamma[[1]], gamma[[1]], FUN = "*")) %*% vec(Sigma_e[[1]])
  prior_P[[1]]= matrix(vec_prior_P, ncol = (p+1)*q, nrow=(p+1)*q)
  prior_V[[1]]=(lambda[[1]]%*%prior_P[[1]]%*%t(lambda[[1]])) + Sigma_u[[1]]

  post_F[[1]]= prior_F[[1]] +prior_P[[1]]%*%t(lambda[[1]])%*% solve(prior_V[[1]])%*% (X[[1]][,1] - prior_x[[1]]) #V anschauen
  post_P[[1]]= prior_P[[1]] - prior_P[[1]]%*%t(lambda[[1]])%*% solve(prior_V[[1]]) %*% lambda[[1]]%*%prior_P[[1]] #nochmal checken

  Kalman.predF <- function(i){
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

  Fsmooth <- vector(mode = "list", length = T-start_ar)
  Psmooth <- vector(mode = "list", length = T-start_ar)

  Fsmooth[[T-start_ar]]<- post_F[[T-start_ar]]
  Psmooth[[T-start_ar]]<- post_P[[T-start_ar]]

  Smoother <- function(i){

    Fsmooth[[i-1]]<<- post_F[[i-1]] + post_P[[i-1]]%*%t(gamma[[1]])%*%solve(prior_P[[i]])%*%(Fsmooth[[i]]- prior_F[[i]])
    Psmooth[[i-1]]<<- post_P[[i-1]] - post_P[[i-1]]%*%t(gamma[[1]])%*%solve(prior_P[[i]])%*%(Psmooth[[i]]-prior_P[[i]])%*%solve(prior_P[[i]])%*%gamma[[1]]%*%post_P[[i-1]]
  }
  lapply((T-start_ar):2,Smoother)
  Fsmooth <<- Fsmooth
  Psmooth <<- Psmooth


}


