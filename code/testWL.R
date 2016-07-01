##########################################################
## Jared Studyvin
## 9 June 2016
## Test the weighted likelihood
##########################################################

library(plyr)

userPath <- '~/GoogleDrive/wind/fatality/areaCorrection/weightedDistance/'
dataPath <- paste0(userPath,'data/')
codePath <- paste0(userPath,'code/')
outPath <- paste0(userPath,'output/')

source(paste0(codePath,'weightFun.R'))
source(paste0(codePath,'weightedLikelihood.R'))
source(paste0(codePath,'getStartValue.R'))



simListWL <- function(dat,weightFun,...){

    ##weightedLikelihood(fatDist=fatDist,fatW=fatW,distribution='gompertz')

    ##simWL(fatDist=fatDist,fatW)

    new <- adply(dat,1,function(row,w,...){
              area <- w(row$distance,type=as.character(row$plotType[1]),...)
              data.frame(row,area)
          },w=weightFun)


    fatW <- with(new,1/(piHat*area)) ## invert the weights
    fatDist <- new$distance


    allDist <- c('rayleigh','gamma','weibull','llog','norm','gompertz')
    out <- ldply(allDist,weightedLikelihood,fatDist=fatDist,fatW=fatW)
        ##print(out)
    out$n <- length(fatDist)


    return(out)


}# end simListWL




########################################################
## gompertz
load(paste0(dataPath,'gompertzList.Rdata'))
gompertzParam

system.time(
    gompertzResult <- llply(gompertzList[1:50],simListWL,weightFun=weightFun)
)
gompertzResult




dat <- extractDist(gompertzResult)





getWeibull <- function(result){
    subset(result,distn=='weibull')
}

weibullEst <- ldply(WLResult,getWeibull)

hist(weibullEst$param1);abline(v=weibullParam[1],col='red')
hist(weibullEst$param2);abline(v=weibullParam[2],col='red')


WLResult[[1]]

table(ldply(WLResult,getBestDist,type='aic'))
table(ldply(WLResult,getBestDist,type='aicc'))
table(ldply(WLResult,getBestDist,type='KS'))


getProb <- function(q,parm,distn){


}


CDF <- function(q,parm,distn){
    cdf <- switch(distn,
                  rayleigh=prayleigh(q=q,scale=parm[1]),
                  gamma=pgamma(q=q,shape=parm[1],rate=parm[2]),
                  weibull=pweibull(q=q,shape=parm[1],scale=parm[2]),
                  llog=pllog(q=q,shape=parm[1],scale=parm[2]),
                  norm=pnorm(q=q,mean=parm[1],sd=parm[2]),
                  gompertz=pgompertz(q=q,scale=parm[1],shape=parm[2])
                  )# end switch
    return(cdf)
}



getAC <- function(result,bestType,wFun){
    bestDist <- getBestDist(result=result,bestType)
    if(bestDist=='FAIL'){
        return(NA)
    }
    est <- eval(parse(text=paste0("subset(result,distn=='",bestDist,"')")))
    estParm <- with(est,c(param1,param2))


    AC <- sum(wFun(1:100)*diff(CDF(0:100,estParm,bestDist)))
    return(AC)
}



ks <- ldply(WLResult,getAC,bestType='KS',wFun=weightFun)
aic <- ldply(WLResult,getAC,bestType='aic',wFun=weightFun)
aicc <- ldply(WLResult,getAC,bestType='aicc',wFun=weightFun)

names(ks) <- 'ks'
names(aic) <- 'aic'
names(aicc) <- 'aicc'

estAC <- cbind(ks,aic,aicc)

boxplot(estAC)

(trueAC <- sum(weightFun(1:100)*diff(CDF(0:100,weibullParam,'weibull'))))
abline(h=trueAC,col='red')





save(WLResult,file=paste0(outPath,'weightedLikeFromWeibullData.Rdata'))


simWD(weibullList[[1]],weightFun,subdivisions=1000)

weightedDistribution('gompertz',fatDist=weibullList[[1]],weightFun=weightFun,subdivisions=1000)


weightedDistribution('norm',fatDist=fatDist,weightFun=w,subdivisions=1000)

