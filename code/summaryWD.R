#################################################
## Jared Studyvin
## 15 July 2016
## summarize the results from the weighted distribution
#################################################
rm(list=ls())

library(plyr)

userPath <- '~/GoogleDrive/wind/fatality/areaCorrection/likelihoodAC/'
dataPath <- paste0(userPath,'data/')
outPath <- paste0(userPath,'output/')
codePath <- paste0(userPath,'code/')



source(paste0(codePath,'weightFun.R'))
##source(paste0(codePath,'weightedDistribution.R'))
##source(paste0(codePath,'getStartValue.R'))
source(paste0(codePath,'summaryFun.R'))


########################################
## llog
########################################
## small
## load(paste0(outPath,'llogSmallWD.Rdata'))
## llogSmallAC <- getAC(resultList=llogSmallResult,wFun=weightFun)
## write.csv(llogSmallAC,paste0(outPath,'areaCorResult/llogSmallWD.csv'),row.names=FALSE)

## mid
load(paste0(outPath,'llogMidWD.Rdata'))
llogMidAC <- getAC(resultList=llogMidResult,wFun=weightFun)
write.csv(llogMidAC,paste0(outPath,'areaCorResult/llogMidWD.csv'),row.names=FALSE)

## big
load(paste0(outPath,'llogBigWD.Rdata'))
llogBigAC <- getAC(resultList=llogBigResult,wFun=weightFun)
write.csv(llogBigAC,paste0(outPath,'areaCorResult/llogBigWD.csv'),row.names=FALSE)

########################################
## norm
########################################
## small
load(paste0(outPath,'normSmallWD.Rdata'))
normSmallAC <- getAC(resultList=normSmallResult,wFun=weightFun)
write.csv(normSmallAC,paste0(outPath,'areaCorResult/normSmallWD.csv'),row.names=FALSE)

## mid
load(paste0(outPath,'normMidWD.Rdata'))
normMidAC <- getAC(resultList=normMidResult,wFun=weightFun)
write.csv(normMidAC,paste0(outPath,'areaCorResult/normMidWD.csv'),row.names=FALSE)

## big
load(paste0(outPath,'normBigWD.Rdata'))
normBigAC <- getAC(resultList=normBigResult,wFun=weightFun)
write.csv(normBigAC,paste0(outPath,'areaCorResult/normBigWD.csv'),row.names=FALSE)

########################################
## weibull
########################################
## small
load(paste0(outPath,'weibullSmallWD.Rdata'))
weibullSmallAC <- getAC(resultList=weibullSmallResult,wFun=weightFun)
write.csv(weibullSmallAC,paste0(outPath,'areaCorResult/weibullSmallWD.csv'),row.names=FALSE)

## mid
load(paste0(outPath,'weibullMidWD.Rdata'))
weibullMidAC <- getAC(resultList=weibullMidResult,wFun=weightFun)
write.csv(weibullMidAC,paste0(outPath,'areaCorResult/weibullMidWD.csv'),row.names=FALSE)

## big
load(paste0(outPath,'weibullBigWD.Rdata'))
weibullBigAC <- getAC(resultList=weibullBigResult,wFun=weightFun)
write.csv(weibullBigAC,paste0(outPath,'areaCorResult/weibullBigWD.csv'),row.names=FALSE)

########################################
## gamma
########################################
## small
load(paste0(outPath,'gammaSmallWD.Rdata'))
gammaSmallAC <- getAC(resultList=gammaSmallResult,wFun=weightFun)
write.csv(gammaSmallAC,paste0(outPath,'areaCorResult/gammaSmallWD.csv'),row.names=FALSE)

## mid
load(paste0(outPath,'gammaMidWD.Rdata'))
gammaMidAC <- getAC(resultList=gammaMidResult,wFun=weightFun)
write.csv(gammaMidAC,paste0(outPath,'areaCorResult/gammaMidWD.csv'),row.names=FALSE)

## big
load(paste0(outPath,'gammaBigWD.Rdata'))
gammaBigAC <- getAC(resultList=gammaBigResult,wFun=weightFun)
write.csv(gammaBigAC,paste0(outPath,'areaCorResult/gammaBigWD.csv'),row.names=FALSE)

