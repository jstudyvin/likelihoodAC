################################################
## Jared Studyvin
## 8 June 2016
## get parameter estimates and AIC for the weighted distribution
################################################

weightedDistribution <- function(distribution,fatDist,weightFun,...){




    require(VGAM)
    require(FAdist)


    distn <- tolower(distribution)
    allowedDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    if(!distn%in%allowedDist){
        stop(paste0('distribution must be one of the following: ',paste(allowedDist,collapse=', ')))
    }

    ## for debugging
    ##print(distn)


    wFat <- weightFun(fatDist,...)

    ## starting values for the optimization
    (startValue <- getStartValue(fatDist,distn,w=wFat))


    toIntegrateFun <- function(x,parm,w,distn,...){
        ## This function is the weight function times the pdf
        ## This function is integrated in the likelihood for the constant
        out <- switch(distn,
                      rayleigh=w(x,...)*drayleigh(x=x,scale=parm[1]),
                      gamma=w(x,...)*dgamma(x=x,shape=parm[1],rate=parm[2]),
                      weibull=w(x,...)*dweibull(x=x,shape=parm[1],scale=parm[2]),
                      llog=w(x,...)*dllog(x=x,shape=parm[1],scale=parm[2]),
                      norm=w(x,...)*dnorm(x=x,mean=parm[1],sd=parm[2]),
                      gompertz=w(x,...)*dgompertz(x=x,scale=parm[1],shape=parm[2])
                      ) # end switch
        return(out)
    }





    ## for testing
    ##plot(fatDist,toIntegrateFun(x=fatDist,startValue,weightFun,distn))


    loglik <- function(parm,x,distn,w,...){

        ## for debugging
        ##print(parm)

        ## all dist have a lower bound of zero except the normal distribution
        lowLim <- ifelse(distn=='norm',-Inf,0)
        upLim <- Inf ## all distribution have an upper bound of infinity

        ## integral value
        c <- tryCatch({
        integrate(toIntegrateFun,lower=lowLim,upper=upLim,parm=parm,w=w,distn=distn,...)$value
        },error=function(cond){
            message('Error with integration in the weighted distribution function.')
            message('The distribution is ',distn,', with parameter values: ',paste(parm,collapse=','))
            return(NA)

        },finally={}) #end tryCatch


        ## get the correct log of f
        logf <- switch(distn,
                       rayleigh=drayleigh(x=x,scale=parm[1],log=TRUE),
                       gamma=dgamma(x=x,shape=parm[1],rate=parm[2],log=TRUE),
                       weibull=dweibull(x=x,shape=parm[1],scale=parm[2],log=TRUE),
                       llog=dllog(x=x,shape=parm[1],scale=parm[2],log=TRUE),
                       norm=dnorm(x=x,mean=parm[1],sd=parm[2],log=TRUE),
                       gompertz=dgompertz(x=x,scale=parm[1],shape=parm[2],log=TRUE)
               ) # end switch

        ## sum the liklihood
        loglik <- sum(logf-log(c)) + sum(log(w(x,...)))
        return(-loglik)
    } # end function loglik

    ## for testing
    ##loglik(startValue,fatDist,distn,weightFun,subdivisions=1000)

    lowLim <- ifelse(distn=='norm',c(-Inf,1e-10),1e-10)

    fit <- nlminb(start=startValue,objective=loglik,x=fatDist,distn=distn,w=weightFun,...,lower=lowLim)

    ## estimated values,make of length 2, if not already
    (fitParm <- c(fit$par,NA)[1:2])


    ## calculate the AIC value
    k <- length(startValue)
    n <- length(fatDist)
    aic <- 2*(k+fit$objective) ## AIC value
    aicc <- aic+2*k*(k+1)/(n-k-1) ## corrected AIC value
    ## I'm not sure this is the best way to compare between distributions



    ## I would like to calculate the kolmogorov-smirnov test statistic
    ## I'm not sure how to calculate the emipirical CDF

    ## ## cdf of fit distn, requiring the parameter estimates
    ## CDF <- function(q,parm,distn,low=-Inf,up=Inf){
    ##     cdf <- switch(distn,
    ##                   rayleigh=prayleigh(q=q,scale=parm[1]),
    ##                   gamma=pgamma(q=q,shape=parm[1],rate=parm[2]),
    ##                   weibull=pweibull(q=q,shape=parm[1],scale=parm[2]),
    ##                   llog=pllog(q=q,shape=parm[1],scale=parm[2]),
    ##                   norm=pnorm(q=q,mean=parm[1],sd=parm[2]),
    ##                   gompertz=pgompertz(q=q,scale=parm[1],shape=parm[2])
    ##                   )# end switch

    ##     return(cdf)
    ## }



    out <- data.frame(distn=distn,param1=fitParm[1],param2=fitParm[2],code=fit$convergence,message=fit$message,aic=aic,aicc=aicc)


    return(out)

} # end weightedDistribution function
