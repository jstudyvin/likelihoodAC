###################################################
## Jared Studyvin
## 26 July 2016
## 998.03
## plot the adjusted fatality count from the simulation
###################################################

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
  print(colSums(spread(nDF,key=distLike,value=bad,fill=0)))

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
