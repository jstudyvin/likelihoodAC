###################################################
## Jared Studyvin
## 23 June 2016
## summary functions for the AC simulation
###################################################



extractDist <- function(resultList,distribution=NULL,...){

    require(plyr)

    if(is.null(distribution)){
        distn <- NULL
    }else{
        distn <- tolower(distribution)
        ##allowedDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
        allowedDist <- c('gamma','weibull','llog','norm')
        if(!distn%in%allowedDist){
            stop(paste0('distribution must be one of the following: ',paste(allowedDist,collapse=', ')))
        }
    }# end else


    extract <- function(el,distn=NULL,criteria='aicc'){

        if(!(criteria%in%names(el))){
            warning('The criteria specified: ', criteria,' is available changing to aicc.')
            criteria <- 'aicc'
        }

        if(is.null(distn)){
            out <- eval(parse(text=paste0('subset(el,code==0&',criteria,'==min(',criteria,'))')))
        }else{
            out <- eval(parse(text=paste0("subset(el,distn=='",distn,"')")))
        }
        return(out)
    } # end extract

    elem <- resultList[[2]]

    extractClass <- function(i,alist,...){
        elem <- alist[[i]]
        if(is.list(elem)){
            out <- ldply(elem,extract,...)
        }else{
            out <- extract(el=elem,...)
        }
        out$rep <- i
        return(out)
    }# end extractClass

    out <- ldply(seq_along(resultList),extractClass,alist=resultList,distn=distn,...)
    ##(out <- ldply(seq_along(resultList),extractClass,alist=resultList,distn='llog'))

    return(out)
} #end extractDist




CDF <- function(q,parm,distn){
    cdf <- switch(distn,
                  rayleigh=prayleigh(q=q,scale=parm[1]),
                  gamma=pgamma(q=q,shape=parm[1],rate=parm[2]),
                  weibull=pweibull(q=q,shape=parm[1],scale=parm[2]),
                  llog=pllog(q=q,shape=parm[1],scale=parm[2]),
                  norm=pnorm(q=q,mean=parm[1],sd=parm[2]),
                  gompertz=pgompertz(q=q,scale=parm[1],shape=parm[2])
                  )# end switch
    return(cdf)
} ## end CDF function




getAC <- function(resultList,wFun,bands=1:100,...){

    dat <- extractDist(resultList=resultList,...)

    out <- adply(dat,1,function(row,q){
                     q <- sort(q)
                     q <- c(head(q,1)-1,q)
                     par <- c(row$param1,row$param2)
                     CDF(q=q,parm=par,distn=row$distn)
                 },q=sort(bands))
    return(out)
} # end getAC
