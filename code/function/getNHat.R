##########################################
## Jared Studyvin
## 26 July 2016
## 998.03
## calculate the adjusted fatality count (NHat) from the simulation results
##########################################

getNHat <- function(resultList,wFun,datList=NULL,...){

    acDat <- getAC(resultList=resultList,wFun=wFun,...)

    ## for debugging
    ##acDat <- getAC(resultList=resultList,wFun=wFun,plotType=plotType)

    getN <- function(acRow,datList){
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

    out <- adply(acDat,1,getN,datList=datList)
    return(out)
} # end getMHat function
