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

    ## for debugging
    ##el <- elem[[4]]

    extract <- function(el,distn=NULL,criteria='aicc'){

        if(!(criteria%in%names(el))){
            warning('The criteria specified: ', criteria,' is unavailable changing to aicc.')
            criteria <- 'aicc'
        }

        ## Infinity aicc values happen when the sample is close to the number of parameters.
        ## This isn't an issue for aic, so switch over if two or more value are infinity
        if(sum(el[,criteria] == Inf,na.rm=TRUE)>1){
            criteria <- 'aic'
        }




        if(is.null(distn)){
            ## negative aic values imply a bad optimization and should be removed
            elConverge <- subset(el,code==0&aic>0)
            out <- eval(parse(text=paste0('subset(elConverge,code==0&',criteria,'==min(',criteria,'))')))
        }else{
            out <- eval(parse(text=paste0("subset(el,distn=='",distn,"')")))
        }

        if(nrow(out)==0){
            out <- el[1,]
            out[,c('distn','param1','param2')] <- NA
        }

        out$distn <- as.character(out$distn)

        return(out)
    } # end extract

    ## for debugging
    ##elem <- resultList[[206]]

    extractClass <- function(i,alist,...){
        elem <- alist[[i]]
        if(class(elem)=='list'){
            out <- ldply(elem,extract,...)
        }else{
            out <- extract(el=elem,...)
        }
        out$rep <- i
        return(out)
    }# end extractClass

    out <- ldply(seq_along(resultList),extractClass,alist=resultList,distn=distn,...)
    ##out <- ldply(132,extractClass,alist=resultList,distn=NULL)

    ## for debugging
    ## out <- NULL
    ## for(i in seq_along(resultList)){
    ##     out <- rbind(out,extractClass(i,resultList))
    ## }

    return(out)
} #end extractDist




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




getAC <- function(resultList,wFun,bands=1:100,plotType=NULL,...){



    dat <- extractDist(resultList=resultList,...)


    if(is.null(plotType)){
        out <- adply(dat,1,function(row,q,wF,...){
                         if(is.na(row$distn)){
                             row$AC <- NA
                             return(row)
                         }
                         q <- sort(q)
                         q <- c(head(q,1)-1,q)
                         par <- c(row$param1,row$param2)
                         fatDen <- diff(CDF(q=q,parm=par,distn=as.character(row$distn)))
                         ##
                         row$AC <- sum(wF(q[-1],type=row$plotType,...)*fatDen)

                         return(row)
                     },q=sort(bands),wF=wFun)

        out$nHat <- with(out,n/(piHat*AC))
    }else{

        names(plotType) <- plotType

        w <- wFun(sort(bands))

        out <- adply(dat,1,function(row,q,wF,pt,...){
                         q <- sort(q)
                         q <- c(head(q,1)-1,q)
                         par <- c(row$param1,row$param2)
                         fatDen <- diff(CDF(q=q,parm=par,distn=as.character(row$distn)))
                         ##
                         ac <- laply(pt,function(x,wF,q,fatDen,...){
                                         ## q was made 1 element longer for the diffence
                                         ## This is why the first element needs removed
                                         sum(wF(q[-1],type=x,...)*fatDen)
                                     },wF=wF,fatDen=fatDen,q=q)
                         for(i in 1:length(pt)){
                             row[,paste0('areaCor',pt[i])] <- ac[i]
                         }

                         return(row)
                     },q=sort(bands),wF=wFun,pt=plotType)
    } # end else
    return(out)
} # end getAC


getMHat <- function(resultList,wFun,datList=NULL,...){

    acDat <- getAC(resultList=resultList,wFun=wFun,...)

    ##acDat <- getAC(resultList=resultList,wFun=wFun,plotType=plotType)

    getM <- function(acRow,datList){
        library(tidyr)
        library(dplyr)
        if(is.null(datList)){
            stop('add code for weighted distribution results')
        }
        dat <- datList[[acRow$rep]]
        ## assumes plotType and piHat are column names
        m <- ddply(dat,~plotType+piHat,summarize,mHat=sum(1/piHat))
        vGather <- paste0(names(acRow)[grepl('areaCor',names(acRow))],collapse=',')
        gatherAC <- eval(parse(text=paste0('gather(acRow,key=plotType,value=AC,',vGather,')')))[,c('plotType','AC')]
        gatherAC$plotType <- gsub('areaCor','',gatherAC$plotType)

        out <- suppressWarnings(full_join(m,gatherAC,by='plotType'))
        out$nHat <- with(out,mHat/AC)

        return(out)

    }# end getM function

    out <- adply(acDat,1,getM,datList=datList)
    return(out)
} # end getMHat function
