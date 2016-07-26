#########################################################
## Jared Studyvin
## 18 July 2016
## 998.030 area correction simulation
## plot the results
#########################################################
rm(list=ls())

library(plyr);library(tidyr);library(dplyr)

userPath <- '~/GoogleDrive/wind/fatality/areaCorrection/likelihoodAC/'
dataPath <- paste0(userPath,'data/')
outPath <- paste0(userPath,'output/')
codePath <- paste0(userPath,'code/function/')

lapply(list.files(codePath,pattern='*.R',full.names=TRUE),source)




acFile <- list.files(paste0(outPath,'areaCorResult/'),pattern='*.csv')

bigFile <- acFile[grepl('Big',acFile)]
midFile <- acFile[grepl('Mid',acFile)]
smallFile <- acFile[grepl('Small',acFile)]
trueFile <- acFile[grepl('True',acFile)]






plotN(files=smallFile,path=paste0(outPath,'areaCorResult/'))
plotN(files=midFile,path=paste0(outPath,'areaCorResult/'))
plotN(files=bigFile,path=paste0(outPath,'areaCorResult/'))

write.csv(data.frame(sampleSize=c('small','mid','big'),gamma=c(7,2,0),llog=c(37,0,1),norm=c(20,8,4),weibull=c(43,2,1)),paste0(outPath,'areaCorResult/nonConvergeWD.csv'),row.names=FALSE)



##path=paste0(outPath,'areaCorResult/')
##files=bigFile
##trueFile

plotAC(files=smallFile,path=paste0(outPath,'areaCorResult/'),trueFile=trueFile)
plotAC(files=midFile,path=paste0(outPath,'areaCorResult/'),trueFile=trueFile)
plotAC(files=bigFile,path=paste0(outPath,'areaCorResult/'),trueFile=trueFile)





## table the best distributions from the simulation
tabDist <- tableDist(file=c(smallFile,midFile,bigFile),path=paste0(outPath,'areaCorResult/'))


write.csv(tabDist,file=paste0(outPath,'tableDistn/tabDistn.csv'),row.names=FALSE,na='0')



tabDist


