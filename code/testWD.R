##########################################################
## Jared Studyvin
## 9 June 2016
## Test the weighted distribution
##########################################################
rm(list=ls())

library(plyr)

source('weightFun.R')
source('weightedDistribution.R')
source('getStartValue.R')
source('summaryFun.R')

userPath <- '~/GoogleDrive/wind/fatality/areaCorrection/weightedDistance/'
dataPath <- paste0(userPath,'data/')
outPath <- paste0(userPath,'output/')


load(paste0(dataPath,'llogSmall.Rdata'))

ls()

##listEl <- llogListSmall[[1]]
##dat <- subset(listEl,plotType=='RP'&piHat==.7)

estWD <- function(listEl,...){

    outList <- dlply(listEl,~plotType+piHat,function(subDat,...){
                       ##allDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
                       allDist <- c('gamma','weibull','llog','norm')
                       plotType <- as.character(subDat$plotType[1])
                       piHat <- subDat$piHat[1]
                       message('The plot type is ',plotType,' and piHat is ',piHat)
                       fatDist <- subDat$distance
                       out <- ldply(allDist,weightedDistribution,fatDist=fatDist,...,type=plotType)
                       out$n <- length(fatDist)
                       out$plotType <- plotType
                       out$piHat <- piHat
                       return(out)
                   },...)

    return(outList)
}


system.time(
    catch <- llply(llogListSmall[1:2],estWD,weightFun=weightFun,subdivisions=10000)
)


extractDist(catch,criteria='aic')

getAC(catch,weightFun)



catch[1]
weightedDistribution('weibull',dat$distance,weightFun=weightFun,subdivisions=1000)

distribution <- 'weibull';fatDist <- dat$distance




source('weightedDistribution.R')
simWD <- function(fatDist,w,...){
    allDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    out <- ldply(allDist,weightedDistribution,fatDist=fatDist,weightFun=w,...)
    ##print(out)
    out$n <- length(fatDist)
    return(out)
}


system.time(
    WDResult <- llply(weibullList,simWD,w=weightFun,subdivisions=1000)
)


save(WDResult,file=paste0(outPath,'weightedDistFromWeibullData.Rdata'))


simWD(weibullList[[1]],weightFun,subdivisions=1000)

weightedDistribution('gompertz',fatDist=weibullList[[1]],weightFun=weightFun,subdivisions=1000)


weightedDistribution('norm',fatDist=fatDist,weightFun=w,subdivisions=1000)

