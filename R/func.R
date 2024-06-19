library(statmod)
library(MASS)
library(simsurv)

#' simulate data from a mixture of Weibull distributions
#'
#' @param N number of observations
#' @param lambdas vector of shape parameters
#' @param gammas vector of scale parameters
#' @param pmix mixture probabilities
#' @param beta vector of regression coefficients
#' @param rateC rate of the censoring distribution
#'
#' @return simulated dataset
#' @export
#' @import survival
#' @import MASS
#' @import simsurv
#' @examples
simulmixweibull <- function(N, lambdas, gammas, pmix, beta, rateC){

  # Generating the covariates
  Sigma =matrix(c(10,3,3,2),2,2)
  x1x2 = as.data.frame(mvrnorm(n=N,mu=c(5,5),Sigma=Sigma) )
  x1 = x1x2[,1]
  x2 = x1x2[,2]
  x3 = (rbinom(N,1,0.8))
  x4 = rep(0,N)
  x4[x3==0] = sample(c(1:4),table(x3)[1], replace = TRUE, prob = c(0.2,0.2,0.3,0.3))
  x4[x3==1] = sample(c(1:4),table(x3)[2], replace = TRUE, prob = c(0.1,0.2,0.4,0.5))

  x_data = as.data.frame(cbind(x1,x2,x3,x4))
  x_data$x3 = factor(x_data$x3)
  x_data$x4 = as.factor(x_data$x4)
  covs = model.matrix(~x1+x2+x3+x4,x_data)
  covs <- data.frame(id = 1:N, covs[,2:dim(covs)[2]])
  # Generate event times using simsurv
  if(length(lambdas)==2){
    s3 <- simsurv(lambdas = lambdas, gammas = gammas,interval = c(1e-08, 10000),betas = c(x1 =beta[1],x2=beta[2],x31 = beta[3],x42= beta[4],x43 = beta[5],x44 =beta[6]),
                  mixture = TRUE, pmix = 0.5, x = covs)
  }else{
    s3 <- simsurv(lambda = lambdas, gamma = gammas,interval = c(1e-08, 10000),betas = c(x1 =beta[1],x2=beta[2],x31 = beta[3],x42= beta[4],x43 = beta[5],x44 =beta[6]),
                  mixture = FALSE, pmix = 0.5, x = covs)

  }


  Tlat <- s3$eventtime
  # Generate censoring time from exponetial distribution
  C <- rexp(n=N, rate=rateC)

  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)

  # Create the simulated data set
  data.frame(id=1:N,
             time=time,
             status=status,
             covs [,2:dim(covs)[2]])
}


#' Caculate the negative log-likelihood function
#'
#' @param par
#' @param t
#' @param d
#' @param X
#' @param boundaryknots
#'
#' @return
#' @export
#' @import splines2
#' @examples
density_surv_spline_eval.hazard<-function(par,t,d,X,boundaryknots){

  n_cov = dim(X)[2]
  n = length(t)
  n_spline_coeff =length(par) - n_cov # number of spline coefficients
  gamma = par[1:n_spline_coeff]
  beta = par[(n_spline_coeff+1):length(par)]
  betaX =as.vector(X%*%beta)
  Basist0=bernsteinPoly(t, degree = n_spline_coeff-1,Boundary.knots
                                 =boundaryknots, derivs = 0, intercept = TRUE)

  intgamma=rep(NA,n)
  for(i in 1:n){

    intgamma[i]=intBerns_eval(0,t[i],gamma,0,0,0,t,d,X, boundaryknots)

  }

  # for gamma's score
  ll=d*(betaX+Basist0%*%gamma)-intgamma*exp(betaX)


  return(-sum(ll))
}
#' Get the first order derivative of the log likelihood function
#'
#' @param par
#' @param t
#' @param d
#' @param X
#' @param boundaryknots
#'
#' @return
#' @export
#'
#' @examples
dloglik_eval.hazard<-function(par,t,d,X,boundaryknots){#par,n_gamma,degree,t,d,X,knots
  n_cov = dim(X)[2]
  n = length(t)
  n_spline_coeff =length(par) - n_cov # number of spline coefficients
  gamma = par[1:n_spline_coeff]
  beta = par[(n_spline_coeff+1):length(par)]
  betaX =as.vector(X%*%beta)

  pbeta = dim(X)[2]
  Basist0=bernsteinPoly(t, degree = n_spline_coeff-1,Boundary.knots
                                  =boundaryknots, derivs = 0, intercept = TRUE)

  b_n = dim(Basist0)[2]
  intgamma=rep(NA,n)
  dintgamma=matrix(NA,ncol=b_n,nrow=n)
  for(i in 1:n){

    intgamma[i]=intBerns_eval(0,t[i],gamma,0,0,0,t,d,X, boundaryknots)
    for(I in 1:b_n){
      dintgamma[i,I]=intBerns_eval(0,t[i],gamma,1,I,0,t,d,X, boundaryknots)

    }
  }


  # for beta's score
  d_beta=d-intgamma*exp(betaX)
  d_beta_matrix=matrix(rep(d_beta,pbeta),ncol=pbeta,byrow=F)
  g_beta=d_beta_matrix*X

  # for gamma's score
  g_gamma=matrix(rep(d,b_n),ncol=b_n,byrow=F)*Basist0-
    matrix(rep(exp(betaX),b_n),ncol=b_n,byrow=F)*dintgamma

  res=cbind(g_gamma,g_beta)



  return(-colSums(res))
}


#' Get the Hessian matrix of the likelihood function.
#'
#' @param par
#' @param t
#' @param d
#' @param X
#' @param boundaryknots
#'
#' @return
#' @export
#'
#' @examples
ddloglik_eval.hazard<-function(par,t,d,X,boundaryknots){
  n_cov = dim(X)[2]
  n = length(t)
  n_spline_coeff =length(par) - n_cov # number of spline coefficients
  gamma = par[1:n_spline_coeff]
  beta = par[(n_spline_coeff+1):length(par)]
  betaX =as.vector(X%*%beta)

  p = dim(X)[2]
  Basist0=bernsteinPoly(t, degree = n_spline_coeff-1,Boundary.knots
                                =boundaryknots, derivs = 0, intercept = TRUE)


  b_n = dim(Basist0)[2]
  inttgamma=rep(NA,n)
  dinttgamma=matrix(NA,ncol=b_n,nrow=n)
  ddinttgamma=array(0,c(b_n,b_n,n))
  start = Sys.time()

  for(i in 1:n){
    inttgamma[i]=intBerns_eval(0,t[i],gamma,0,0,0,t,d,X, boundaryknots)

    for(I in 1:b_n){
      dinttgamma[i,I]=intBerns_eval(0,t[i],gamma,1,I,0,t,d,X, boundaryknots)
      for(J in 1:b_n){
        ddinttgamma[J,I,i]=intBerns_eval(0,t[i],gamma,2,J,I,t,d,X, boundaryknots)
      }
    }
  }
  end = Sys.time()


  # for \gamma\gamma
  # coeff for XX^T
  expbetaX = exp(betaX)
  g_beta=-inttgamma*expbetaX

  # for \gamma\alpha
  g_gamma=-expbetaX

  ddbb=matrix(0,nrow=p,ncol=p)
  ddbg=matrix(0,nrow=p,ncol=b_n)
  dgg=matrix(0,nrow=b_n,ncol=b_n)
  for(i in 1:n){
    ddbb=ddbb+g_beta[i]*as.numeric(X[i,])%*%t(as.numeric(X[i,]))
    ddbg=ddbg+g_gamma[i]*(matrix(as.numeric(X[i,]),ncol=1)%*%matrix(dinttgamma[i,],nrow=1))
    dgg=dgg+g_gamma[i]*ddinttgamma[,,i]
  }


  res_1=cbind(dgg,t(ddbg))
  res_2=cbind((ddbg),ddbb)
  res=rbind(res_1,res_2)

  return(-res)

}

#' Numerical integration for the polynomial
#'
#' @param lb
#' @param ub
#' @param gamma
#' @param k
#' @param I
#' @param J
#' @param t
#' @param d
#' @param X
#' @param breaks
#'
#' @return
#' @export
#' @import statmod
#' @examples
intBerns_eval<-function(lb,ub,gamma,k,I,J,t,d,X,breaks){

  out=gauss.quad.prob(10,"uniform",l=lb, u=ub)
  fx=Berns_eval(out$nodes,gamma,k,I,J,t,d,X, breaks)
  ans=(ub-lb)*sum(out$weights * fx)
  return(ans)
}

#' Evaluate Bernstein Polyomials
#'
#' @param x0
#' @param gamma
#' @param k
#' @param I
#' @param J
#' @param t
#' @param d
#' @param X
#' @param boundaryknots
#'
#' @return
#' @export
#'
#' @examples

Berns_eval<-function(x0,gamma,k,I,J,t,d,X,boundaryknots){
  dg = length(gamma)-1
  Basis0=splines2::bernsteinPoly(x0, degree = dg,Boundary.knots
                                 =boundaryknots, derivs = 0, intercept = TRUE)
  eg=as.numeric(exp(Basis0%*%gamma))
  ans=NULL
  if(k==0 && I==0 && J==0){ans=eg}
  if(k==1 && I!=0 && J==0){ans=eg*Basis0[,I]}
  if(k==2 && I!=0 && J!=0){ans=eg*Basis0[,I]*Basis0[,J]}
  return(ans)
}

#' Finding initial value for COLSA
#'
#' @param dg
#' @param data_first
#' @param form
#' @param boundaryknots
#'
#' @return
#' @export
#' @import flexsurv
#' @examples
find_inits<-function(dg,data_first,form,boundaryknots){

  splineDesign_mtrix = splines2::bernsteinPoly(data_first$time, degree = dg,
                                               Boundary.knots=boundaryknots, derivs = 0, intercept = TRUE)

  fit = flexsurv::flexsurvspline(form,data =data_first,k = dg)
  n_gamma = 2+dg
  hazard = hsurvspline(data_first$time,gamma= fit$coefficients[1:(n_gamma)],knots = fit$knots,scale="hazard")
  initial_val = as.matrix(lm(log(hazard) ~ splineDesign_mtrix-1)$coefficients)[,1]
  return(c(initial_val, fit$coefficients[(n_gamma+1):length(fit$coefficients)]))
}


#' Get the estimated survival probabilities
#'
#' @param par
#' @param t
#' @param X
#' @param boundaryknots
#'
#' @return
#' @export
#'
#' @examples
get_est_surv<-function(par,t,X,boundaryknots){

  n_cov = dim(X)[2]
  n = length(t)
  n_spline_coeff =length(par) - n_cov # number of spline coefficients
  gamma = par[1:n_spline_coeff]
  beta = par[(n_spline_coeff+1):length(par)]
  betaX =as.vector(X%*%beta)

  intgamma=rep(NA,n)

  for(i in 1:n){
    intgamma[i]=intBerns_eval(0,t[i],gamma,0,0,0,t,d,X, boundaryknots)

  }
  exp(-intgamma*exp(betaX))
}
save_data<-function(tempdir,K){
  for (i in c(1:K)){
    subdata= subset(data,group==i)
    tempdatadir = paste0(tempdir,"/Simdata/hospital",i)
    if(!dir.exists(tempdatadir)){
      dir.create(tempdatadir, recursive = TRUE)
    }
    save(subdata, file = paste(tempdatadir, "/Simdata.RData", sep=""))
  }
}
