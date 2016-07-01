##################################################
## Jared Studyvin
## 22 June 2016
## Test the generateData function
##################################################
rm(list=ls());graphics.off()


userPath <- '~/GoogleDrive/wind/fatality/areaCorrection/weightedDistance/'
codePath <- paste0(userPath,'code/')
dataPath <- paste0(userPath,'data/')


## This is the weight function
## The second argument corresponds to the plot type
source(paste0(codePath,'weightFun.R'))

## This is the function to generate the data
source(paste0(codePath,'createFatData.R'))


plotParm <- function(parm,distn,x=0:100){
    require(truncdist)
    d <- switch(distn,
                       rayleigh=drayleigh(x=x,scale=parm[1],log=FALSE),
                       gamma=dgamma(x=x,shape=parm[1],rate=parm[2],log=FALSE),
                       weibull=dweibull(x=x,shape=parm[1],scale=parm[2],log=FALSE),
                       llog=dllog(x=x,shape=parm[1],scale=parm[2],log=FALSE),
                       norm=dtrunc(x=x,spec='norm',mean=parm[1],sd=parm[2],log=FALSE),
                       gompertz=dgompertz(x=x,scale=parm[1],shape=parm[2],log=FALSE)
               ) # end switch

    plot(x,d,main=distn)
}


plotN <- function(l,...){
    require(plyr)

    sampleSize <- laply(l,nrow)
    hist(sampleSize,xlab='sample size',...)
    print(summary(sampleSize))
    return(NA)
}





########################### set true parameter values
## gamma distribution
gammaParam <- c(2,.05) #shape and rate
## weibull distribution
weibullParam <- c(1.5,50) #shape and scale
## llog distribution
llogParam <- c(.4,3.5) #shape and scale
## norm distribution
normParam <- c(35,30) #mean and sd

## plot the distribution with the true parameter values
windows();plotParm(gammaParam,'gamma')
windows();plotParm(weibullParam,'weibull')
windows();plotParm(llogParam,'llog')
windows();plotParm(normParam,'norm')

propRP <- .9;propArrive <- .7


###################################################################################
## everything set above needs to remain fixed for the simulation
###################################################################################






## input data frame for the generateData function
## Needs these three column names exactly
## The values of plotType are passed as the second argument to the weight function
## different values of piHat are implying seasonal differences
## n is the sample size for the true number of fatalities


#################### mid sample size
midN <- 2000;

(piTypeCombMid <- data.frame(plotType=c(rep('RP',2),rep('FULL',2)),piHat=c(.7,.9,.4,.6),n=c(propRP*midN*propArrive,propRP*midN*(1-propArrive),(1-propRP)*midN*propArrive,(1-propRP)*midN*(1-propArrive))))

################# ## generate data
## gamma distribution
set.seed(45);gammaListMid <- createFatData(distribution='gamma',trueParam=gammaParam,piTypeComb=piTypeCombMid,nRep=1000,weightFun=weightFun)
## weibull distribution
set.seed(637);weibullListMid <- createFatData(distribution='weibull',trueParam=weibullParam,piTypeComb=piTypeCombMid,nRep=1000,weightFun=weightFun)
## llog distribution
set.seed(98);llogListMid <- createFatData(distribution='llog',trueParam=llogParam,piTypeComb=piTypeCombMid,nRep=1000,weightFun=weightFun)
## norm distribution
set.seed(23);normListMid <- createFatData(distribution='norm',trueParam=normParam,piTypeComb=piTypeCombMid,nRep=1000,weightFun=weightFun)

################################
## plot the sample sizes for each
windows();plotN(gammaListMid,n=100,main='mid gamma')
windows();plotN(weibullListMid,n=100,main='mid weibull')
windows();plotN(llogListMid,n=100,main='mid llog')
windows();plotN(normListMid,n=100,main='mid norm')



####################################### small sample size
smallN <- 500;

(piTypeCombSmall <- data.frame(plotType=c(rep('RP',2),rep('FULL',2)),piHat=c(.7,.9,.4,.6),n=c(propRP*smallN*propArrive,propRP*smallN*(1-propArrive),(1-propRP)*smallN*propArrive,(1-propRP)*smallN*(1-propArrive))))

################# ## generate data
## gamma distribution
set.seed(458);gammaListSmall <- createFatData(distribution='gamma',trueParam=gammaParam,piTypeComb=piTypeCombSmall,nRep=1000,weightFun=weightFun)
## weibull distribution
set.seed(6374);weibullListSmall <- createFatData(distribution='weibull',trueParam=weibullParam,piTypeComb=piTypeCombSmall,nRep=1000,weightFun=weightFun)
## llog distribution
set.seed(986);llogListSmall <- createFatData(distribution='llog',trueParam=llogParam,piTypeComb=piTypeCombSmall,nRep=1000,weightFun=weightFun)
## norm distribution
set.seed(3);normListSmall <- createFatData(distribution='norm',trueParam=normParam,piTypeComb=piTypeCombSmall,nRep=1000,weightFun=weightFun)

################################
## plot the sample sizes for each
windows();plotN(gammaListSmall,n=100,main='small gamma')
windows();plotN(weibullListSmall,n=100,main='small weibull')
windows();plotN(llogListSmall,n=100,main='small llog')
windows();plotN(normListSmall,n=100,main='small norm')




####################################### big sample size
bigN <- 5000;

(piTypeCombBig <- data.frame(plotType=c(rep('RP',2),rep('FULL',2)),piHat=c(.7,.9,.4,.6),n=c(propRP*bigN*propArrive,propRP*bigN*(1-propArrive),(1-propRP)*bigN*propArrive,(1-propRP)*bigN*(1-propArrive))))

################# ## generate data
## gamma distribution
set.seed(458);gammaListBig <- createFatData(distribution='gamma',trueParam=gammaParam,piTypeComb=piTypeCombBig,nRep=1000,weightFun=weightFun)
## weibull distribution
set.seed(6374);weibullListBig <- createFatData(distribution='weibull',trueParam=weibullParam,piTypeComb=piTypeCombBig,nRep=1000,weightFun=weightFun)
## llog distribution
set.seed(986);llogListBig <- createFatData(distribution='llog',trueParam=llogParam,piTypeComb=piTypeCombBig,nRep=1000,weightFun=weightFun)
## norm distribution
set.seed(3);normListBig <- createFatData(distribution='norm',trueParam=normParam,piTypeComb=piTypeCombBig,nRep=1000,weightFun=weightFun)

################################
## plot the sample sizes for each
windows();plotN(gammaListBig,n=100,main='big gamma')
windows();plotN(weibullListBig,n=100,main='big weibull')
windows();plotN(llogListBig,n=100,main='big llog')
windows();plotN(normListBig,n=100,main='big norm')





#############################################################
## save the data

save(gammaListBig,piTypeCombBig,gammaParam,file=paste0(dataPath,'gammaBig.Rdata'))
save(weibullListBig,piTypeCombBig,weibullParam,file=paste0(dataPath,'weibullBig.Rdata'))
save(llogListBig,piTypeCombBig,llogParam,file=paste0(dataPath,'llogBig.Rdata'))
save(normListBig,piTypeCombBig,normParam,file=paste0(dataPath,'normBig.Rdata'))


save(gammaListMid,piTypeCombMid,gammaParam,file=paste0(dataPath,'gammaMid.Rdata'))
save(weibullListMid,piTypeCombMid,weibullParam,file=paste0(dataPath,'weibullMid.Rdata'))
save(llogListMid,piTypeCombMid,llogParam,file=paste0(dataPath,'llogMid.Rdata'))
save(normListMid,piTypeCombMid,normParam,file=paste0(dataPath,'normMid.Rdata'))


save(gammaListSmall,piTypeCombSmall,gammaParam,file=paste0(dataPath,'gammaSmall.Rdata'))
save(weibullListSmall,piTypeCombSmall,weibullParam,file=paste0(dataPath,'weibullSmall.Rdata'))
save(llogListSmall,piTypeCombSmall,llogParam,file=paste0(dataPath,'llogSmall.Rdata'))
save(normListSmall,piTypeCombSmall,normParam,file=paste0(dataPath,'normSmall.Rdata'))







