######################################################################
## Jared Studyvin
## 13 June 2016
## Generate data for simulation
######################################################################


createFatData <- function(distribution,trueParam,piTypeComb,nRep,weightFun,...){
### distribution: indicate the truth
### trueParam: vector of the parameters for distribution
### piTypeComb: dataframe with columns plotType, piHat, n
### Needs these three column names exactly
### The values of plotType are passed as the second argument to the weight function
### different values of piHat are implying seasonal differences
### n is the sample size for the true number of fatalities
### nRep: indicates the number of data sets to be return in the list
### weightFun: first argument should be the distance, second argument should take the values of plotType in piTypeComb
### ...: additional arguments to pass to weightFun


    require(VGAM)
    require(FAdist)
    require(plyr)
    require(truncdist)

    ## check is the distribution is allowed or not
    distn <- tolower(distribution)
    allowedDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    if(!distn%in%allowedDist){
        stop(paste0('distribution must be one of the following: ',paste(allowedDist,collapse=', ')))
    }


    ## generates observed data for a single combination of piTypeComb
    getData <- function(row,parm,distn,w,...){
        n <- row$n[1]
        type <- as.character(row$plotType[1])
        piHat <- row$piHat[1]
        ## get a random sample from the distribution of interest
        out <- switch(distn,
                      rayleigh=rrayleigh(n=n,scale=parm[1]),
                      gamma=rgamma(n=n,shape=parm[1],rate=parm[2]),
                      weibull=rweibull(n=n,shape=parm[1],scale=parm[2]),
                      llog=rllog(n=n,shape=parm[1],scale=parm[2]),
                      norm=rtrunc(n=n,spec='norm',a=0,mean=parm[1],sd=parm[2]),
                      gompertz=rgompertz(n=n,scale=parm[1],shape=parm[2])
                      ) # end switch

        out <- ceiling(out) # observed data comes as whole numbers
        ## for debugging
        ## cat('w:',w(75,type,...),'\n')
        ## cat('type =',type,'\n')
        ## cat('piHat =',piHat,'\n')
        ## cat('n =',n,'\n')


        ## which observations to keep based on the weight function w and the value of piHat
        keep <- rbinom(n=n,size=1,w(out,type,...)*piHat) ==1

        if(sum(keep)==0){
            return(NA)
        }

        return(out[keep])

    } #end getData function

    ## combined plotType and piHat into a single indexing column
    piTypeComb$comb <- unlist(alply(piTypeComb,1,function(row){paste0(row$plotType[1],row$piHat[1],collapse='')}))



    getRep <- function(index,piTypeComb,parm,distn,w,...){
        ## apply the getData function across the rows of piTypeComb
        obsDatList <- dlply(piTypeComb,~comb,getData,parm=parm,distn=distn,w=w,...)
        ##(obsDatList <- dlply(piTypeComb,~comb,getData,parm=parm,distn=distn,w=w))
        ##print(index)
        ## reformat list into a nice data frame
        obsDat <- adply(names(obsDatList),1,function(name,obsDatList,piTypeComb){
                            row <- eval(parse(text=paste0("subset(piTypeComb,comb=='",name,"')"))) ## subset piTypeComb this combination
                            data.frame(distance=obsDatList[[name]],plotType=row$plotType[1],piHat=row$piHat[1])
          },obsDatList=obsDatList,piTypeComb=piTypeComb)

        obsDat$X1 <- NULL

        ## NA happens when no observation happen in a strata
        outDat <- subset(obsDat,!is.na(distance))

        return(outDat)
    }# end getRep function

    ## repeat the above process nRep times
    dataList <- llply(1:nRep,getRep,piTypeComb=piTypeComb,parm=trueParam,distn=distn,w=weightFun,...)

    return(dataList)

} # end createFatData function




