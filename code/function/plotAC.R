###################################################
## Jared Studyvin
## 26 July 2016
## 998.03
## plot the area correction value from the simulation
###################################################


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

    y.lim <- c(0,.2)
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
