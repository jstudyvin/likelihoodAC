##########################################
## Jared Studyvin
## 26 July 2016
## 998.03
## extract the distribution from each simulation iteration
##########################################


## Several distribution are fit is each iteration and this function extracts the best based on some critera or will always extract the named distribution

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
