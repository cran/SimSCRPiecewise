#'SimSCRPiecewise
#'This function simulates semi-competing risks data based on three piecewise exponential hazards, three true regression vectors and three matrices of patient covariates (which can be different or the same).
#'This simulates from the semi-markov model of Lee et al (2015) given patient covariates, regression parameters and baseline hazard functions.
#'@param x1 - Matrix of patient covariates for hazard 1 simulation
#'@param x2 - Matrix of patient covariates for hazard 2 simulation
#'@param x3 - Matrix of patient covariates for hazard 3 simulation
#'@param beta1 - vector of size ncol(x1) that is the true regression coefficient vector for hazard 1
#'@param beta2 - vector of size ncol(x2) that is the true regression coefficient vector for hazard 2
#'@param beta3 - vector of size ncol(x3) that is the true regression coefficient vector for hazard 3
#'@param s1 - vector of size at least length 2, where the first entry is 0. This characterizes the split point locations of baseline hazard 1
#'@param s2 - vector of size at least length 2, where the first entry is 0. This characterizes the split point locations of baseline hazard 2
#'@param s3 - vector of size at least length 2, where the first entry is 0. This characterizes the split point locations of baseline hazard 3
#'@param lam1 - vector of the same size as s1. This vector is the true baseline hazard 1 heights and the last entry represents the height on the interval [max(s1),infinity)
#'@param lam2 - vector of the same size as s2. This vector is the true baseline hazard 2 heights and the last entry represents the height on the interval [max(s2),infinity)
#'@param lam3 - vector of the same size as s3. This vector is the true baseline hazard 3 heights and the last entry represents the height on the interval [max(s3),infinity)
#'@param gamma - vector containing patient frailties.
#'@param cens - This is the administrative right censoring time of the study. All patients who have survival outcomes after cens have survival times set to cens.
#'@return Returns a list of size 4 containing the semi-competing risks simulated data. Entry 1 contains the non-terminal event times for the patients.
#'Entry 2 contains the terminal event times for the patients. Entry 3 contains the patient indicators for whether or not a patient experienced a non-terminal event prior to death.
#'Entry 4 contains the patient indicators for whether or not they experienced a terminal event.
#'@import stats
#'@examples
#' ##Set number of patients and covariate matrices
#'n=100
#'x1=matrix(rnorm(n*10,0,1),nrow=n)
#'x2=x1
#'x3=x1
#'##Sets up true covariate vectors
#'beta1=rnorm(10,0,1)
#'beta2=rnorm(10,0,1)
#'beta3=c(3,rep(0,9))
#'##Sets up three baseline hazard split locations
#'s1=c(0,7,30,100,1000)
#'s2=c(0,50,100,2000)
#'s3=c(0,10,40,50,500)
#'##Sets up baseline hazard heights
#'lam1=c(.1,.1,.3,.1,.1)
#'lam2=c(.2,.3,.1,.1)
#'lam3=c(.1,.3,.2,.2,.1)
#'gamma=rgamma(100,1,1)
#'##Runs Function and returns a list of simulated data
#'X=SimSCRPiecewise(x1,x2,x3,beta1,beta2,beta3,s1,s2,s3,lam1,lam2,lam3,gamma,1000)
#'X
#'
#'@references
#' Lee, K. H., Haneuse, S., Schrag, D. and Dominici, F. (2015), Bayesian semiparametric analysis of semicompeting risks data: investigating hospital readmission after a pancreatic cancer diagnosis. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 253-273.
#'
#'@export
SimSCRPiecewise=function(x1,x2,x3,beta1,beta2,beta3,s1,s2,s3,lam1,lam2,lam3,gamma,cens){

Y1=rep(NA,nrow(x1))
Y2=Y1
I1=Y1
I2=Y1

bounds1=rep(NA,(length(s1)-1))
bounds2=rep(NA,(length(s2)-1))
bounds3=rep(NA,(length(s3)-1))

diff1=diff(s1)
diff2=diff(s2)
diff3=diff(s3)


for(b in 1:length(Y1)){
  sum1=0

  for(k in 1:length(bounds1)){
    sum1=sum1+gamma[b]*exp(x1[b,]%*%beta1)*diff1[k]*lam1[k]

    bounds1[k]=exp(-sum1)
  }



  sum2=0

  for(k in 1:length(bounds2)){
    sum2=sum2+gamma[b]*exp(x2[b,]%*%beta2)*diff2[k]*lam2[k]

    bounds2[k]=exp(-sum2)
  }

  bounds1=1-bounds1
  bounds2=1-bounds2


  if(length(s1)==2){

    U1=runif(1,0,1)
    spot=0

    A=exp(x1[b,]%*%beta1)*gamma[b]


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


    A=exp(x1[b,]%*%beta1)*gamma[b]

    if(spot==0){
      Y1[b]=-(log(1-U1))/(A*lam1[1])
    }else{
      Y1[b]=-(log(1-U1)+A*diff1[1:(spot-1)]%*%lam1[1:(spot-1)])/(A*lam1[spot])
      Y1[b]=Y1[b]+s1[spot]

    }


  }






  U1=runif(1,0,1)

  if(length(s2)==2){

    spot=0

    A=exp(x2[b,]%*%beta2)*gamma[b]


    if(U1>bounds2){
      spot=1
    }

    if(spot==0){
      Y2[b]=-(log(1-U1))/(A*lam2[1])
    }else{
      Y2[b]=-(log(1-U1)+diff2*A*lam2[1])/(A*lam2[2])+s2[2]

    }


  }else{



    if(U1<bounds2[1]){
      spot=0
    }else{

      for(k in 1:(length(bounds2)-1)){
        if(U1>bounds2[k] && U1 < bounds2[k+1]){
          spot=k+1
        }

      }

    }


    A=exp(x2[b,]%*%beta2)*gamma[b]

    if(spot==0){
      Y2[b]=-(log(1-U1))/(A*lam2[1])
    }else{
      Y2[b]=-(log(1-U1)+A*diff2[1:(spot-1)]%*%lam2[1:(spot-1)])/(A*lam2[spot])+s2[spot]

    }

  }



  if(Y1[b]>cens &&Y2[b]>cens){
    Y1[b]=cens
    Y2[b]=cens
    I1[b]=0
    I2[b]=0
  }else{



    if(Y2[b]<Y1[b]){
      if(Y2[b]<cens){
        Y1[b]=Y2[b]
        I1[b]=0
        I2[b]=1
      }else{
        Y1[b]=cens
        Y2[b]=cens
        I1[b]=0
        I2[b]=0
      }
    }else{

      if(Y1[b]<cens && Y2[b]>cens){
        I1[b]=1
        I2[b]=0
        Y2[b]=cens

      }else{


        sum3=0

        for(k in 1:length(bounds3)){
          sum3=sum3+gamma[b]*exp(x3[b,]%*%beta3)*diff3[k]*lam3[k]

          bounds3[k]=exp(-sum3)
        }

        bounds3=1-bounds3


        U1=runif(1,0,1)

        if(length(bounds3)==1){

          spot=0

          A=exp(x3[b,]%*%beta3)*gamma[b]


          if(U1>bounds3){
            spot=1
          }

          if(spot==0){
            M=-(log(1-U1))/(A*lam3[1])
          }else{
            M=-(log(1-U1)+diff3*A*lam3[1])/(A*lam3[2])+s3[2]

          }



        }else{

          if(U1<bounds3[1]){
            spot=0
          }else{

            for(k in 1:(length(bounds3)-1)){
              if(U1>bounds3[k] && U1 < bounds3[k+1]){
                spot=k+1
              }

            }

          }


          A=exp(x3[b,]%*%beta3)*gamma[b]

          if(spot==0){
            M=-(log(1-U1))/(A*lam3[1])
          }else{

            M=-(log(1-U1)+A*diff3[1:(spot-1)]%*%lam3[1:(spot-1)])/(A*lam3[spot])+s3[spot]
          }

        }

        Y2[b]=Y1[b]+M



        if(Y2[b]>cens){
          Y2[b]=cens
          I2[b]=0
          I1[b]=1
        }else{
          I2[b]=1
          I1[b]=1
        }



      }



    }


  }


















}

z=list(Y1,Y2,I1,I2)

return(z)

}




