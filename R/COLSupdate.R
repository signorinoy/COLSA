#' COLSA update function for demo
#'
#' @param B number of sites/subdatas
#' @param data a dataframe that contains all datasets
#' @param init initial beta estimates
#' @param p the total number of covaraites
#' @param npar the number of covariates in the model
#' @param boundaryknots the boundary knots for the Bernstein polynomial
#'
#' @return
#' @export
#'
#' @examples
updateCOLSA.demo <-

  function(B, data, init=NA, p,npar,boundaryknots){

    sum2<-diag(0,p,p);
    tol=1e-5;
    max_iter=10000;
    negl<-0
    if(sum(is.na(init))>1){init<-rep(0,p)}
    betahat<-init;

    for(b in 1:B){
      subdata= subset(data,group==b)
      t <- subdata$time
      d <- subdata$status
      X = as.matrix(subdata[,4:(4+npar-1)])





      betahat_old=betahat;

      H<-ddloglik_eval.hazard(betahat,t,d,X,boundaryknots)

      U=chol(sum2+H)
      L=t(U)

      for (r in 1:max_iter){

        g_0=-dloglik_eval.hazard(betahat,t,d,X,boundaryknots)
        g_1=-t(crossprod((betahat-betahat_old),sum2));

        g=g_0+g_1;

        d_beta=backsolve(U,forwardsolve(L,g))
        df_beta=crossprod(g,d_beta);


        if (abs(df_beta)<tol){
          break
        }else {
          betahat=betahat+d_beta;
        }


        H<-ddloglik_eval.hazard(betahat,t,d,X,boundaryknots) #J2

        U=chol(sum2+H)
        L=t(U)
      }

      H_new<-ddloglik_eval.hazard(betahat,t,d,X,boundaryknots)

      sum2<-sum2+H_new
      negl<-negl+density_surv_spline_eval.hazard(betahat,t,d,X,boundaryknots)
    }


      sd<-sqrt(diag(solve(sum2)));

    pvalue<-2*pnorm(-abs(betahat)/sd)
    result<-cbind(betahat=betahat,sd=sd,pvalue=pvalue,negll=negl)#
    colnames(result)<-c("Estimates","Std.Errors","p-values","neg-logll")#
    if(r == max_iter){
      stop("maximum iterations reached")
      result<-cbind(betahat=rep(NA,length(betahat)),sd=rep(NA,length(sd)),pvalue=rep(NA,length(pvalue)),negll=rep(NA,length(negl)))#

    }
    return(result)

  }
#' COLSA update function for distributed datasets
#'
#' @param subdata a dataframe that contains the data for a single site
#' @param statistics a list that contains the current estimates of beta and Hessian
#' @param p the total number of parameters
#' @param npar the number of covariates in the model
#' @param boundaryknots the boundary knots for the Bernstein polynomial
#'
#' @return
#' @export
#'
#' @examples
updateCOLSA.outloop <-
  function(subdata, statistics, p,npar,boundaryknots){


    tol=1e-5;
    max_iter=10000;
    negl = statistics$negl

    init = statistics$betahat
    sum2 = statistics$Hessian

    if(sum(is.na(init))>1){init<-rep(0,p)}
    betahat<-init;
    if(sum(is.na(sum2))>=1){sum2<-diag(0,p,p)}


      t <- subdata$time
      d <- subdata$status
      X = as.matrix(subdata[,4:(4+npar-1)])





      betahat_old=betahat;

      H<-ddloglik_eval.hazard(betahat,t,d,X,boundaryknots)

      U=chol(sum2+H)
      L=t(U)

      for (r in 1:max_iter){

        g_0=-dloglik_eval.hazard(betahat,t,d,X,boundaryknots)
        g_1=-t(crossprod((betahat-betahat_old),sum2));

        g=g_0+g_1;

        d_beta=backsolve(U,forwardsolve(L,g))
        df_beta=crossprod(g,d_beta);


        if (abs(df_beta)<tol){
          break
        }else {
          betahat=betahat+d_beta;
        }


        H<-ddloglik_eval.hazard(betahat,t,d,X,boundaryknots) #J2

        U=chol(sum2+H)
        L=t(U)
      }

      H_new<-ddloglik_eval.hazard(betahat,t,d,X,boundaryknots)

      sum2<-sum2+H_new
      negl<-negl+density_surv_spline_eval.hazard(betahat,t,d,X,boundaryknots)



    sd<-sqrt(diag(solve(sum2)));

    pvalue<-2*pnorm(-abs(betahat)/sd)
    result<-cbind(betahat=betahat,sd=sd,pvalue=pvalue,negll=negl)#
    colnames(result)<-c("Estimates","Std.Errors","p-values","neg-logll")#
    if(r == max_iter){
      stop("maximum iterations reached")
      result<-cbind(betahat=rep(NA,length(betahat)),sd=rep(NA,length(sd)),pvalue=rep(NA,length(pvalue)),negll=rep(NA,length(negl)))#

    }
    statistics = list(betahat = betahat,Hessian = sum2,negl=negl)
    return(list(result=result,statistics = statistics))

  }
