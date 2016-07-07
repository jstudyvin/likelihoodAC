##########################################################
## Jared Studyvin
## 5 July 2016
## estimated the weighted likelihood on the simulated data
##########################################################

rm(list=ls())

library(plyr)

userPath <- '~/GoogleDrive/wind/fatality/areaCorrection/likelihoodAC/'
dataPath <- paste0(userPath,'data/')
codePath <- paste0(userPath,'code/')
outPath <- paste0(userPath,'output/')

source(paste0(codePath,'weightFun.R'))
source(paste0(codePath,'weightedLikelihood.R'))
source(paste0(codePath,'getStartValue.R'))



simListWL <- function(dat,weightFun,...){

    require(plyr)

    new <- adply(dat,1,function(row,w,...){
              area <- w(row$distance,type=as.character(row$plotType[1]),...)
              data.frame(row,area)
          },w=weightFun)

    ## the weights need to have a multiplicative effect on the likelihood
    fatW <- with(new,1/(piHat*area)) ## invert the weights
    fatDist <- new$distance


    ##allDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    allDist <- c('gamma','weibull','llog','norm')
    out <- ldply(allDist,weightedLikelihood,fatDist=fatDist,fatW=fatW)
        ##print(out)
    out$n <- length(fatDist)


    return(out)


}# end simListWL



library(parallel)
library('doSNOW')
detectCores(all.tests = TRUE, logical = TRUE)




#######################################################################
## gamma distribution
#######################################################################
## small

load(paste0(dataPath,'gammaSmall.Rdata'))



cl <- makeCluster(4)
registerDoSNOW(cl)
Sys.time()
system.time(
    gammaSmallResult <- llply(gammaListSmall,simListWL,weightFun=weightFun,.parallel=TRUE,.paropts=list(.export=c('weightedLikelihood','getStartValue')))
    )
stopCluster(cl)

save(gammaSmallResult,gammaParam,file=paste0(outPath,'gammaSmallWL.Rdata'))


load(paste0(dataPath,'gammaMid.Rdata'))
cl <- makeCluster(4)
registerDoSNOW(cl)
Sys.time()
system.time(
    gammaMidResult <- llply(gammaListMid,simListWL,weightFun=weightFun,.parallel=TRUE,.paropts=list(.export=c('weightedLikelihood','getStartValue')))
    )
stopCluster(cl)

save(gammaMidResult,gammaParam,file=paste0(outPath,'gammaMidWL.Rdata'))


load(paste0(dataPath,'gammaBig.Rdata'))
cl <- makeCluster(4)
registerDoSNOW(cl)
Sys.time()
system.time(
    gammaBigResult <- llply(gammaListBig,simListWL,weightFun=weightFun,.parallel=TRUE,.paropts=list(.export=c('weightedLikelihood','getStartValue')))
    )
stopCluster(cl)

save(gammaBigResult,gammaParam,file=paste0(outPath,'gammaBigWL.Rdata'))



#######################################################################
## weibull distribution
#######################################################################
## small

load(paste0(dataPath,'weibullSmall.Rdata'))


cl <- makeCluster(4)
registerDoSNOW(cl)
Sys.time()
system.time(
    weibullSmallResult <- llply(weibullListSmall,simListWL,weightFun=weightFun,.parallel=TRUE,.paropts=list(.export=c('weightedLikelihood','getStartValue')))
    )
stopCluster(cl)

save(weibullSmallResult,weibullParam,file=paste0(outPath,'weibullSmallWL.Rdata'))


load(paste0(dataPath,'weibullMid.Rdata'))
cl <- makeCluster(4)
registerDoSNOW(cl)
Sys.time()
system.time(
    weibullMidResult <- llply(weibullListMid,simListWL,weightFun=weightFun,.parallel=TRUE,.paropts=list(.export=c('weightedLikelihood','getStartValue')))
    )
stopCluster(cl)

save(weibullMidResult,weibullParam,file=paste0(outPath,'weibullMidWL.Rdata'))


load(paste0(dataPath,'weibullBig.Rdata'))
cl <- makeCluster(4)
registerDoSNOW(cl)
Sys.time()
system.time(
    weibullBigResult <- llply(weibullListBig,simListWL,weightFun=weightFun,.parallel=TRUE,.paropts=list(.export=c('weightedLikelihood','getStartValue')))
    )
stopCluster(cl)

save(weibullBigResult,weibullParam,file=paste0(outPath,'weibullBigWL.Rdata'))




