#####################################################################
## Jared Studyvin
## 13 June 2016
## Weighted likelihood
#####################################################################


weightedLikelihood <- function(fatDist,fatW,distribution,...){


    havevgam <- require(VGAM)
    if(!havevgam){
        install.packages('VGAM')
    }
    havefadist <- require(FAdist)
    if(!havefadist){
        install.packages('FAdist')
    }



    distn <- tolower(distribution)
    allowedDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    if(!distn%in%allowedDist){
        stop(paste0('distribution must be one of the following: ',paste(allowedDist,collapse=', ')))
    }


    if(length(fatW)!=length(fatDist)){
        fatW <- rep(1,length(fatDist))
        message('The length of fatW is not same as the length of fatDist and will be ignored in calculating the weight likelihood parameter estimates.')
    }




    wLogLik <- function(parm,x,w,distn,low=-Inf,up=Inf){
        ## lik <- switch(distn,
        ##              rayleigh=dtrunc(x=x,spec=distn,a=low,b=up,scale=parm[1]),
        ##                gamma=dtrunc(x=x,spec=distn,a=low,b=up,shape=parm[1],rate=parm[2]),
        ##                weibull=dtrunc(x=x,spec=distn,a=low,b=up,shape=parm[1],scale=parm[2]),
        ##                llog=dtrunc(x=x,spec=distn,a=low,b=up,shape=parm[1],scale=parm[2]),
        ##                norm=dtrunc(x=x,spec=distn,a=low,b=up,mean=parm[1],sd=parm[2]),
        ##                gompertz=dtrunc(x=x,spec=distn,a=low,b=up,scale=parm[1],shape=parm[2])
        ##              )# end switch

        lik <- switch(distn,
                      rayleigh=drayleigh(x=x,scale=parm[1]),
                      gamma=dgamma(x=x,shape=parm[1],rate=parm[2]),
                      weibull=dweibull(x=x,shape=parm[1],scale=parm[2]),
                      llog=dllog(x=x,shape=parm[1],scale=parm[2]),
                      norm=dnorm(x=x,mean=parm[1],sd=parm[2]),
                      gompertz=dgompertz(x=x,scale=parm[1],shape=parm[2])
                      )# end switch

        wLL <- sum(w*log(lik))
        return(-wLL)
    }



    ## starting values for the optimization
    (startValue <- getStartValue(x=fatDist,distribution=distn,w=fatW))


    ## lower limits of the parameters
    ## Only the mean from the normal distribution is unconstrained the rest have to be greater than 0.
    lowLim <- ifelse(distn=='norm',c(-Inf,1e-10),1e-10)
    ##fit <- nlminb(start=startValue,objective=wLogLik,x=fatDist,distn=distn,w=fatW,lower=lowLim)
    fit <- nlminb(start=startValue,objective=wLogLik,x=fatDist,distn=distn,w=fatW,...,lower=lowLim)

    ## estimated values,make of length 2, if not already
    (fitParm <- c(fit$par,NA)[1:2])


    ## weighted empirical cdf
    ecdf <- function(q,data,w){
        require(plyr)
        out <- aaply(q,1,function(t,d,w){
                  sum(w*(d<=t))/sum(w)
              },d=data,w=w)
        return(out)
    }

    ## cdf of fit distn, requiring the parameter estimates
    CDF <- function(q,parm,distn,low=-Inf,up=Inf){
        ## cdf <- switch(distn,
        ##               rayleigh=ptrunc(q=q,spec=distn,a=low,b=up,scale=parm[1]),
        ##               gamma=ptrunc(q=q,spec=distn,a=low,b=up,shape=parm[1],rate=parm[2]),
        ##               weibull=ptrunc(q=q,spec=distn,a=low,b=up,shape=parm[1],scale=parm[2]),
        ##               llog=ptrunc(q=q,spec=distn,a=low,b=up,shape=parm[1],scale=parm[2]),
        ##               norm=ptrunc(q=q,spec=distn,a=low,b=up,mean=parm[1],sd=parm[2]),
        ##               gompertz=ptrunc(q=q,spec=distn,a=low,b=up,scale=parm[1],shape=parm[2])
        ##               )# end switch

        ## non truncated distributions
        cdf <- switch(distn,
                      rayleigh=prayleigh(q=q,scale=parm[1]),
                      gamma=pgamma(q=q,shape=parm[1],rate=parm[2]),
                      weibull=pweibull(q=q,shape=parm[1],scale=parm[2]),
                      llog=pllog(q=q,shape=parm[1],scale=parm[2]),
                      norm=pnorm(q=q,mean=parm[1],sd=parm[2]),
                      gompertz=pgompertz(q=q,scale=parm[1],shape=parm[2])
                      )# end switch


        return(cdf)
    }

    ## For the Kolmogorov¡VSmirnov
    ksStat <- function(q,parm,distn,data,w){
        -abs(CDF(q=q,parm=parm,distn=distn)-ecdf(q=q,data=data,w=w))
    }

    ## for debugging
    ##plot(CDF(1:100,parm=fitParm,distn=distn),type='l');points(ecdf(1:100,data=fatDist,w=fatW))


    ksFit <- nlminb(start=mean(fatDist),ksStat,parm=fitParm,distn=distn,data=fatDist,w=fatW)

    (KS <- -ksFit$objective)



    ## calculate the AIC value
    k <- length(startValue)
    n <- length(fatDist)
    aic <- 2*(k+fit$objective) ## AIC value
    aicc <- aic+2*k*(k+1)/(n-k-1) ## corrected AIC value
    ## I'm not sure this is the best way to compare between distributions

    out <- data.frame(distn=distn,param1=fitParm[1],param2=fitParm[2],code=fit$convergence,message=fit$message,aic=aic,aicc=aicc,KS=KS)

    return(out)

} # end weighteLikelihood function
