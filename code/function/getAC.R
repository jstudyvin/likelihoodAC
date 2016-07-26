##########################################
## Jared Studyvin
## 26 July 2016
## 998.03
## calculate the area correction value from the simulation results
##########################################

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
