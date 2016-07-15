#################################################
## Jared Studyvin
## 15 July 2016
## summarize the results from the weighted likelihood
#################################################

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
## gamma
########################################
