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
codePath <- paste0(userPath,'code/')


acFile <- list.files(paste0(outPath,'areaCorResult/'))

bigFile <- acFile[grepl('Big',acFile)]
midFile <- acFile[grepl('Mid',acFile)]
smallFile <- acFile[grepl('Small',acFile)]
trueFile <- acFile[grepl('True',acFile)]




plotN <- function(files,path){

    names(files) <- gsub('.csv','',files)

    datList <- llply(files,function(x,p){
                         read.csv(paste0(p,x))
                     },p=path)

    ## trueList <- llply(trueFile,function(x,p){
    ##                      read.csv(paste0(p,x))
    ##                  },p=path)


    ##   names(datList)
    ##   d <- head(subset(datList[[1]],plotType=='RP'))
    ## spread(d,piHat,value=AC)

    ##d <- head(datList[[1]],100)
    resultList <- llply(datList,function(d){ddply(d,~rep,summarize, nhat=sum(nHat,na.rm=TRUE),bad=sum(is.na(nHat)))})

    nDF <- ldply(resultList,rbind,.id='distLike')
    subset(nDF,bad!=0)
    colSums(nSpread <- spread(nDF,key=distLike,value=nhat))
   ## return(colSums(spread(nDF,key=distLike,value=bad,fill=0)))

    if(any(grepl('Big',files))){
        hline <- 5000
        y.lim <- c(0,2*hline)
        plotName <- 'bigN.pdf'
    }else if(any(grepl('Mid',files))){
        hline <- 2000
        y.lim <- c(0,2.5*hline)
        plotName <- 'midN.pdf'
    }else{
        hline <- 500
        y.lim <- c(0,3*hline)
        plotName <- 'smallN.pdf'
    }



    pdf(paste0(path,plotName))
    par(mfrow=c(2,2))
    boxplot(nSpread[,grepl('gamma',names(nSpread))],ylim=y.lim,ylab='N')
    abline(h=hline,col='red')
    boxplot(nSpread[,grepl('llog',names(nSpread))],ylim=y.lim,ylab='N')
    abline(h=hline,col='red')
    boxplot(nSpread[,grepl('weibull',names(nSpread))],ylim=y.lim,ylab='N')
    abline(h=hline,col='red')
    boxplot(nSpread[,grepl('norm',names(nSpread))],ylim=y.lim,ylab='N')
    abline(h=hline,col='red')
    graphics.off()

    return(NA)

} # end plotN function



plotN(files=smallFile,path=paste0(outPath,'areaCorResult/'))
plotN(files=midFile,path=paste0(outPath,'areaCorResult/'))
plotN(files=bigFile,path=paste0(outPath,'areaCorResult/'))


##path=paste0(outPath,'areaCorResult/')
##files=bigFile
##trueFile
plotAC <- function(files,path,trueFile){

    names(files) <- gsub('.csv','',files)

    datList <- llply(files,function(x,p){
                         read.csv(paste0(p,x))
                     },p=path)

    names(trueFile) <- gsub('.csv','',trueFile)
    trueDF <- ldply(trueFile,function(x,p){
                         read.csv(paste0(p,x))
                     },p=path,.id='distn')

    trueDFRP <- subset(trueDF,plotType=='RP')


    ##head(datList[[1]]);head(datList[[2]])

    ##nameList <- names(datList)

    ##(nameEl <- nameList[2])


    nrow(dat <- datList[[4]])

    getACData <- function(dat){
            ac7 <- subset(dat,plotType=='RP'&piHat==.7,select=c(rep,AC))
            names(ac7) <- c('rep','RP7')
            ac9 <- subset(dat,plotType=='RP'&piHat==.9,select=c(rep,AC))
            names(ac9) <- c('rep','RP9')
            out <- full_join(ac7,ac9)
            return(out)
        }


    resultList <- llply(datList,getACData)


    rList <- llply(names(resultList),function(y,r){
                       x <- r[[y]]
                       names(x) <- paste0(y,names(x))
                       return(x[,!grepl('rep',names(x))])
                   },r=resultList)



    acDF <- do.call(cbind,rList)

    acDF <- acDF[,!grepl('WLRP9',names(acDF))]


    if(any(grepl('Big',names(acDF)))){
        plotName <- 'areaCorBigN.pdf'
    }else if(any(grepl('Mid',names(acDF)))){
        plotName <- 'areaCorMidN.pdf'
    }else{
        plotName <- 'areaCorSmallN.pdf'
    }

    y.lim <- NULL#c(0,.15)
 pdf(paste0(path,plotName),width=13)
##    windows(width=13)
    par(mfrow=c(2,2))
    boxplot(acDF[,grepl('gamma',names(acDF))],ylim=y.lim)
    abline(h=trueDFRP[with(trueDFRP,grepl('gamma',distn)),'A'],col='red')
    boxplot(acDF[,grepl('llog',names(acDF))],ylim=y.lim)
    abline(h=trueDFRP[with(trueDFRP,grepl('llog',distn)),'A'],col='red')
    boxplot(acDF[,grepl('weibull',names(acDF))],ylim=y.lim)
    abline(h=trueDFRP[with(trueDFRP,grepl('weibull',distn)),'A'],col='red')
    boxplot(acDF[,grepl('norm',names(acDF))],ylim=y.lim)
    abline(h=trueDFRP[with(trueDFRP,grepl('norm',distn)),'A'],col='red')
    graphics.off()

    return(NA)

} # end plotN function



plotAC(files=smallFile,path=paste0(outPath,'areaCorResult/'),trueFile=trueFile)
plotAC(files=midFile,path=paste0(outPath,'areaCorResult/'),trueFile=trueFile)
plotAC(files=bigFile,path=paste0(outPath,'areaCorResult/'),trueFile=trueFile)
