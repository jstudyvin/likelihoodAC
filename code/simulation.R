###########################################
## Jared Studyvin
## 22 May 2016
## test script
###########################################

w <- function(x){
    require(dplyr)
    x <- ceiling(x)

    prop <- data.frame(d=1:100,prop=c(0.994938718149504,0.994624021920978,0.97823470461618,
                               0.903887671489263,0.781446817754437,0.660080011358875,
                                   0.536081699880172,0.409099290100651,0.29840175576812,
                                   0.189740873334674,0.111264971851594,0.0816700493858669,
                                   0.067700278957639,0.059958175937043,0.0545341895120893,
                                   0.0505828469268522,0.0475734802409185,0.0450099381748921,
                                   0.0423378182198427,0.0401874841234347,0.0385724876256597,
                                   0.0375685607674148,0.0379679064582935,0.0401752489037758,
                                   0.0424131471715669,0.0438726225903946,0.0447347583214495,
                                   0.0433469339543116,0.0412068218256184,0.0385248279209736,
                                   0.0365510108255672,0.0345598836947837,0.0324967362727387,
                                   0.0306221406330178,0.0293123303643584,0.0280691902067884,
                                   0.0269997951054828,0.0262282943445054,0.0254393320872947,
                                   0.0248405591169834, 0.0244960828525636,0.0242036187013642,
                                   0.0233586683041639,0.0226035346599462,0.0222185158984456,
                                   0.0214069613363088,0.020683810236012,0.0200167922196213,
                                   0.019411577915632,0.0187682283245286, 0.0182352956331142,
                                   0.0179184182081402,0.0178280820006544,0.0175513843953515,
                                   0.0173705435405511,0.016963286718304,0.0165774119064628,
                                   0.0161678252991839,0.0157737313237982,0.0153829454765745,
                                   0.0150233199416864,0.0147224581429039,0.014358043163838,
                                   0.0142845815296265,0.0141121743926314,0.0139173993754729,
                                   0.0136788604784495,0.013534509906191,0.0132257057163886,
                                   0.013059391311518,0.01309305396692,0.0130782401190816,
                                   0.0128978650101229,0.0127050541932945,0.0125310864539905,
                                   0.012462809770892,0.0125492535518906,0.0130525249770394,
                                   0.0132565365366079,0.0131972031553393,0.0131075146272192,
                                   0.0126218287678565,0.0120599362155872,0.0117114471331286,
                                   0.0114894197137408,0.011247687525456,0.0109519229419088,
                                   0.0109130541927116,0.0109879441900024,0.0110231274119703,
                                   0.0110033252854863,0.0108244349124521,0.0103993318123575,
                                   0.0102663216566701,0.0101761125809184,0.0100221948109964,
                                   0.00959313241747452,0.00889708968216072,0.00766000805175065,
                                   0.00598333373939794))

    out <- left_join(data.frame(d=x),prop,by=c('d'='d'))$prop

    out[is.na(out)] <- 1e-10

    return(out)


}



############## create data
alpha <- 2;beta <- 20
set.seed(29);y <- rgamma(1000,alpha,scale=beta)

plot(y,dgamma(y,alpha,1/beta))

allFat <- data.frame(dist=round(y),prob=w(round(y)))

## remove some data based on prob detection
set.seed(736);allFat$keep <- with(allFat,rbinom(length(dist),1,prob))
obsFat <- subset(allFat,keep==1,select=c(dist,prob));nrow(obsFat)
summary(obsFat)

distance <- obsFat$dist ## distance data

###############################################################
## maximum weighted likelihood


wloglik <- function(theta,x,inverseW=FALSE){

    if(any(theta<=0)){
        return(0)
    }

    if(inverseW){
        loglik <- sum(dgamma(x,shape=theta[1],scale=theta[2],log=TRUE)/w(x))
    }else{
        loglik <- sum(w(x)*dgamma(x,shape=theta[1],scale=theta[2],log=TRUE))
    }

    return(-loglik)
}

alpha*beta


wloglik(c(1,4),x=distance)

(start <- c(mean(distance)^2/var(distance),var(distance)/mean(distance)))

fit1 <- nlminb(start,wloglik,x=distance,lower=1e-100)

fit2 <- nlminb(start,wloglik,x=distance,inverseW=TRUE,lower=1e-100)

windows()
hist(distance,freq=FALSE,main='Maximum Weighted Likelihood',ylim=c(0,.15),n=50)
plotseq <- seq(1,100,by=.5)
points(plotseq,dgamma(plotseq,alpha,scale=beta),type='l')
points(plotseq,dgamma(plotseq,fit1$par[1],scale=fit1$par[2]),type='l',col='red')
points(plotseq,dgamma(plotseq,fit2$par[1],scale=fit2$par[2]),type='l',col='green')
legend('topright',legend=c('Truth','w: prob det','w: inverse prob det'),col=c('black','red','green'), lty=1)

mwle <- data.frame(rbind(c(alpha,beta),fit1$par,fit2$par))
names(mwle) <- c('alpha','beta')
rownames(mwle) <- c('truth','probDet','invProbDet')
mwle

###################################################################
## weighted distribution

h <- function(x,theta,inverseW=FALSE){

    if(inverseW){
        ##w <- w(x)
        ##w[w==0] <- -99
        out <- dgamma(x,shape=theta[1],scale=theta[1])/w(x)
        ##out[out<0] <- 0
    }else{
        out <- w(x)*dgamma(x,shape=theta[1],scale=theta[1])
    }
    return(out)
}

wDistLogLik <- function(theta,x,iW){

    c <- integrate(h,lower=0,upper=Inf,theta=theta,inverseW=iW,subdivisions=1000)$value

    log.f <- sum(dgamma(x=x,shape=theta[1],scale=theta[2],log=TRUE))

    loglik <- log.f-length(x)*log(c)

    return(-loglik)

}


wfit1 <- nlminb(start,wDistLogLik,x=distance,iW=FALSE)
wfit1
wfit2 <- nlminb(start,wDistLogLik,x=distance,iW=TRUE)
wfit2


wDistEst <- data.frame(rbind(c(alpha,beta),wfit1$par,wfit2$par))
names(wDistEst) <- c('alpha','beta')
rownames(wDistEst) <- c('truth','probDet','invProbDet')
wDistEst



wDist <- function(x,theta,inverseW=FALSE,noWeight=FALSE){
    c <- integrate(h,lower=0,upper=Inf,theta=theta,inverseW=inverseW,subdivisions=1000)$value


    if(noWeight){
        return(dgamma(x,shape=theta[1],scale=theta[2])/c)
    }

    if(inverseW){
        num <- dgamma(x,shape=theta[1],scale=theta[2])/w(x)
    }else{
        num <- dgamma(x,shape=theta[1],scale=theta[2])*w(x)
    }

    fstar <- num/c
    return(fstar)
}

windows()
hist(distance,freq=FALSE,main='Weighted Distribution',ylim=c(0,.15),n=50)
plotseq <- seq(1,100,by=.5)
points(plotseq,dgamma(plotseq,alpha,scale=beta),type='l')
points(plotseq,wDist(plotseq,wfit1$par),type='l',col='red')
points(plotseq,wDist(plotseq,wfit2$par,inverseW=TRUE),type='l',col='green')


points(plotseq,wDist(plotseq,wfit1$par,noWeight=TRUE),type='l',col='blue')
points(plotseq,wDist(plotseq,wfit2$par,inverseW=TRUE,noWeight=TRUE),type='l',col='blue')



points(plotseq,dgamma(plotseq,wfit1$par[1],scale=wfit1$par[2]),type='l',lty=2,col='red')
points(plotseq,dgamma(plotseq,wfit2$par[1],scale=wfit2$par[2]),type='l',lty=2,col='green')
legend('topright',legend=c('Truth','w dist, prob det','w dist, inverse prob det','original dist, prob det', 'original dist, inverse prob det'),col=c('black','red','green','red','green'), lty=c(1,1,1,2,2))




library(MASS)
mle <- fitdistr(distance,'gamma')
mle
1/mle$estimate[2]

mwle
wDistEst


###################
## A for weighted likelihood

integrate(h,lower=0,upper=Inf,theta=fit2$par,inverseW=FALSE,subdivisions=1000)$value

length(y)

sum(y<=100)


#######################

1/sum(wDist(1:100,wfit2$par,inverseW=TRUE,noWeight=TRUE))
