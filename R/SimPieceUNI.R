#'SimUnivPiecewise
#'
#' This function simulates univariate survival data from a piecewise exponential model with a proportional hazards assumption given a covariate matrix, true beta vector, baseline hazard splits, baseline hazard heights and a right censoring time.
#'@param x1 - Matrix of patient covariates for hazard 1 simulation
#'@param beta1 - vector of size ncol(x1) that is the true regression coefficient vector for the baseline hazard function
#'@param s1 - vector of size at least length 2, where the first entry is 0. This characterizes the split point locations of baseline hazard
#'@param lam1 - vector of the same size as s1. This vector is the true baseline hazard heights and the last entry represents the height on the interval [max(s1),infinity)
#'@param cens - This is the administrative right censoring time of the study. All patients who have survival outcomes after cens have survival times set to cens.
#'@return Returns a list of size 4 containing the semi-competing risks simulated data. Entry 1 contains the non-terminal event times for the patients.
#'Entry 2 contains the terminal event times for the patients. Entry 3 contains the patient indicators for whether or not a patient experienced a non-terminal event prior to death.
#'Entry 4 contains the patient indicators for whether or not they experienced a terminal event.
#'@import stats
#'@examples
#' ##Set number of patients and covariate matrices
#'n=100
#'x1=matrix(rnorm(n*10,0,1),nrow=n)
#'##Sets up true covariate vector
#'beta1=rnorm(10,0,1)
#'##Sets up true baseline hazard split locations
#'s1=c(0,7,30,100,1000)
#'##Sets up baseline hazard heights
#'lam1=c(.1,.1,.3,.1,.1)
#'##Runs Function and returns a list of simulated data
#'X=SimUNIVPiecewise(x1,beta1,s1,lam1,1000)
#'X
#'
#'@references
#' Lee, K. H., Haneuse, S., Schrag, D. and Dominici, F. (2015), Bayesian semiparametric analysis of semicompeting risks data: investigating hospital readmission after a pancreatic cancer diagnosis. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 253-273.
#'
#'@export
SimUNIVPiecewise=function(x1,beta1,s1,lam1,cens){

  Y1=rep(NA,nrow(x1))
  I1=Y1

  bounds1=rep(NA,(length(s1)-1))


  diff1=diff(s1)



  for(b in 1:length(Y1)){
    sum1=0

    for(k in 1:length(bounds1)){
      sum1=sum1+exp(x1[b,]%*%beta1)*diff1[k]*lam1[k]

      bounds1[k]=exp(-sum1)
    }



    bounds1=1-bounds1


    if(length(s1)==2){

      U1=runif(1,0,1)
      spot=0

      A=exp(x1[b,]%*%beta1)


      if(U1>bounds1){
        spot=1
      }

      if(spot==0){
        Y1[b]=-(log(1-U1))/(A*lam1[1])
      }else{
        Y1[b]=-(log(1-U1)+diff1*A*lam1[1])/(A*lam1[2])+s1[2]

      }


    }else{

      U1=runif(1,0,1)

      if(U1<bounds1[1]){
        spot=0



      }else{

        for(k in 1:(length(bounds1)-1)){
          if(U1>bounds1[k] && U1 < bounds1[k+1]){
            spot=k+1
          }

        }

      }


      A=exp(x1[b,]%*%beta1)

      if(spot==0){
        Y1[b]=-(log(1-U1))/(A*lam1[1])
      }else{
        Y1[b]=-(log(1-U1)+A*diff1[1:(spot-1)]%*%lam1[1:(spot-1)])/(A*lam1[spot])
        Y1[b]=Y1[b]+s1[spot]

      }


    }







if(Y1[b]>cens){
  Y1[b]=cens
  I1[b]=0
}else{
  I1[b]=1
}


  }

  z=list(Y1,I1)

  return(z)

}

