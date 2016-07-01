########################################################
## Jared Studyvin
## 14 June 2016
## get starting values for the weighted likelihood or weighted distribution
########################################################


getStartValue <- function(x,distribution,w){

    require(SDMTools) ## for weighted variance and mean

    distn <- tolower(distribution)
    allowedDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    if(!distn%in%allowedDist){
        stop(paste0('distribution must be one of the following: ',paste(allowedDist,collapse=', ')))
    }

    if(length(w)!=length(x)){
        w <- rep(1,length(x))
        message('The length of w is not same as the length of x and will be ignored in calculated the start values in the getStarValue function.')
    }
    ## weighted mean and variance
    (wmean <- wt.mean(x,w))
    (wvar <- wt.var(x,w))
    (maxx <- max(w*x))




    gompertzStart <- function(xbar,varx){

        eulerConst <- -digamma(1) ## Euler-Mascheroni Constant
        g <- function(b,xbar,varx,eulerConst){
            ## approximate solution see Lenart 2012
            a.b <- pi^2/12 - b^2*varx/2
            1/b * exp(a.b)*(a.b - log(a.b) - eulerConst)-xbar
        }
        (maxg <- sqrt(pi^2/(6*varx))) ## max value

        low <- 1e-5;up <- maxg-1e-5

        ##y <- seq(low,up,length=100)
        ##g(y,eulerConst=eulerConst,xbar=xbar,varx=varx)


suppressWarnings(b <- nlminb(low,g,eulerConst=eulerConst,xbar=xbar,varx=varx,lower=1e-5)$par)

        ## find the root
        ##b <- uniroot(g,interval=c(low,up),eulerConst=eulerConst,xbar=xbar,varx=varx)$root
        ## use b to find a
        a <- b*pi^2/12 - b^3/2 * varx
        c(b,a)
        return(c(b,a)) ## scale and shape
    }


    weibullStart <- function(xbar,varx,maxx){
        ## method of moments for weibull parameters
        wei <- function(b,varx,xbar){
            out <- exp(lgamma(1+2/b)-2*lgamma(1+1/b)) - 1 - varx/xbar^2
            return(out)
            }
        b <- uniroot(wei,interval=c(1e-1,maxx),varx=varx,xbar=xbar)$root
        a <- xbar/gamma(1+1/b)
        return(c(b,a))
    }


    ## This function gets the starting values for the distribution of interest.
    out <- switch(distn,
                  rayleigh=c(wmean/sqrt(pi/2)), ## scale
                  gamma=c(wmean^2/wvar,mean(x)/var(x)), ## shape and rate
                  weibull=weibullStart(xbar=wmean,varx=wvar,maxx=max(x)), ## shape and scale
                  llog=c(pi/sqrt(3*wvar),log(wmean)), ## shape and scale
                  norm=c(wmean,sqrt(wvar)), ## mean and sd
                  gompertz=gompertzStart(xbar=wmean,varx=wvar) ## scale and shape
                  ) #end switch
    return(out)
}

