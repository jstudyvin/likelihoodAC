##########################################
## Jared Studyvin
## 26 July 2016
## 998.03
## generic cdf function for the distributions of interest for the AC
##########################################

CDF <- function(q,parm,distn){
    library(VGAM)
    library(FAdist)
    cdf <- switch(distn,
                  rayleigh=prayleigh(q=q,scale=parm[1]),
                  gamma=pgamma(q=q,shape=parm[1],rate=parm[2]),
                  weibull=pweibull(q=q,shape=parm[1],scale=parm[2]),
                  llog=pllog(q=q,shape=parm[1],scale=parm[2]),
                  norm=pnorm(q=q,mean=parm[1],sd=parm[2]), ## this might need to be truncated.
                  gompertz=pgompertz(q=q,scale=parm[1],shape=parm[2])
                  )# end switch
    return(cdf)
} ## end CDF function
