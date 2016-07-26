### A simulation comparison of weighted likelihood vs. wMLE
#1. fits
#2. estimates of dwp

### need for installing gpclib
## pathRtools <- paste(c("C:\\Rtools\\bin",
##  "C:\\Rtools\\mingw_32\\bin",
##  "C:\\Program Files (x86)\\MiKTeX 2.9\\miktex\\bin",
##  "C:\\Program Files\\R\\R-3.3.0\\bin\\i386",
##  "C:\\Windows",
##  "C:\\Windows\\System32"), collapse=";")
## Sys.setenv(PATH=paste(pathRtools,Sys.getenv("PATH"),sep=";"))
## install.packages('gpclib',type='source')

library(splancs) # library with functions for determining whether or not points lie in a polygon
library(gpclib) # ...polygon clippling library can calculate the area of intersecting polygons

#### process:
### STEP 1: define the ground configuration and search area for a turbine
RPpoly<-cbind( #polygon defining RP out to 200m at one turbine (#3 at Fowler)
c(9.528,11.04,12.249,13.122,13.636,13.776,13.329,12.563,11.496,10.156,8.578,6.801,4.871,2.838,0.754,-1.329,-3.355,-5.275,-12.998,-113.914,-115.552,-117.152,-118.674,-120.078,-121.328,-122.391,-123.24,-123.852,-124.212,-123.655,-127.485,-129.102,-129.875,-130.52,-131.34,-132.317,-133.427,-134.643,-135.936,-137.275,-138.628,-199.092,-199.605,-199.618,-134.509,-131.363,-127.924,-124.212,-118.548,-44.406,-44.406,-16.087,-10.939,-9.325,-7.466,-5.411,-3.216,-0.939,1.357,3.611,5.762,7.752,9.528),
c(8.912,7.183,5.229,3.104,0.865,-1.428,-3.47,-5.414,-7.211,-8.816,-10.185,-11.286,-12.089,-12.574,-12.728,-12.548,-12.038,-11.211,-10.181,-8.121,-8.219,-8.578,-9.189,-10.037,-11.099,-12.348,-13.752,-15.274,-16.874,-157.073,-154.103,-152.677,-26.657,-25.467,-24.391,-23.454,-22.68,-22.086,-21.688,-21.493,-21.508,-17.992,-12.558,-12.152,-16.359,-12.169,-8.215,-4.517,-1.428,-2.458,10.929,11.444,7.84,9.476,10.825,11.852,12.529,12.837,12.769,12.325,11.518,10.37,8.912)
)
RPpoly

maxd<-200
#win.graph()
par(mar=c(.5,.5,.5,.5))
plot(0,0,xlim=maxd*c(-1,1),ylim=maxd*c(-1,1),xlab="x",ylab='y')
theta<-seq(0,2*pi,length=1000)
lines(maxd*sin(theta),maxd*cos(theta))
polygon(RPpoly,col=colors()[124],density=NA)

### STEP 2: calculate weights as fraction of area in RP in 0.5m annuli (ignoring turbine base)
## The algorithm:
#1. write circles as polygons
cply<-function(nside,r){ # function to format a circle with radius r as a polygon with nside sides
   theta<-seq(0,2*pi,length=nside+1)
   cbind(r*cos(theta),r*sin(theta))
}
#2. calculate areas of intersection of disks of radius r = seq(rincr, 200, by=rincr) and RP polygon (for Fowler #3)
rincr<-0.5; nside<-1000
wx<-array(dim=c(maxd/rincr,2))
wx[,1]<-seq(rincr,maxd,by=rincr) # rincr m rings with outer boundary at wx[,1]
wc<-numeric(dim(wx)[1]) # area of RP within radius = rincr * i
RP<-as(RPpoly,'gpc.poly') # formatting a RP polygon as a gpc polygon to use with gpc functions
for (i in 1:(maxd/rincr)){
   wc[i]<-area.poly(intersect(as(cply(nside,i*rincr),'gpc.poly'),RP)) # area in RP within a radius of r = i * rincr
}
#3. calculate differences in disk areas to get ring areas in RP and divide by total area in the ring
wx[,2]<-(wc-c(0,wc[1:(length(wc)-1)]))/(pi*rincr*(wx[,1]+c(0,wx[1:(dim(wx)[1]-1)])))
wx<-rbind(wx,c(Inf,0)) # no chance of observing carcasses past 200m

### STEP 3: generate true carcass distribution and determine which lands on RP
## One example...
# In preliminary model, assume f ~ gamma(shape, scale) with family known and parameters unknown
#  reasonable gamma parameters are (correspond to y0 = 0.05 in dwp paper)
a<-2.04459; b<-19.15065

# toss the carcasses
nfat<-500 # total number of carcasses
set.seed(2827);fatd<-rgamma(nfat,shape=a,scale=b)
set.seed(2798);ang<-2*runif(nfat)*pi # directions carcasses fly
points(fatd*cbind(cos(ang),sin(ang)))

# which carcasses land on RP and are observed?
ocarc<-inout(fatd*cbind(cos(ang),sin(ang)), RPpoly, bound=T)
# distances of observed carcasses (supposing searcher efficiency = 1 on RP and 0 off RP)
fato<-fatd[ocarc] # the actual locations are not tracked; only the distance and whether or not they land on RP

points((fatd*cbind(cos(ang),sin(ang)))[ocarc,],col='red')
# check:
check.carcasses<-function(fato = fato, fatd = fatd, RPpoly = RPpoly, maxd = 200){
  plot(0,0,xlim=maxd*c(-1,1),ylim=maxd*c(-1,1),xlab="x",ylab='y')
  theta<-seq(0,2*pi,length=1000)
  lines(maxd*sin(theta),maxd*cos(theta))
  polygon(RPpoly,col=colors()[124],density=NA)
  points(fatd*cos(ang),fatd*sin(ang),pch='*',col=1+4*ocarc)
}


## library(MASS)
## fitdistr(fatd,'gamma');1/fitdistr(fatd,'gamma')$estimate[2]
## fitdistr(fato,'gamma');1/fitdistr(fato,'gamma')$estimate[2]


## cons<-function(gparm,fato,wx){
##     log(sum(wx[,2]*dgamma(wx[,1],shape=gparm[1],scale=gparm[2])))
## }

## cons(c(a,b),fato,wx)

w <- function(x,wx,small=1e-10){
    require(plyr)
    require(dplyr)
    xw <- data.frame(wx)
    names(xw) <- c('y','prop')

    (y <- round_any(x,.5))
    df <- data.frame(y)
    out <- left_join(df,xw,by=c('y'='y'))$prop
    out[is.na(out)] <- small
    return(out)
}

set.seed(2828);ocarc <- rbinom(nfat,size=1,prob=w(fatd,wx,0))==1
fato<-fatd[ocarc]

## dont <- rnorm(10)
## dont
## w(dont,wx)


## cons2 <- function(x,gparm,wx){
##     w(x,wx,0)*dgamma(x,shape=gparm[1],scale=gparm[2])
## }


## cc <- integrate(cons2,lower=0,upper=Inf,gparm=c(a,b),wx=wx,subdivisions=1000)
## log(cc$value)


## llkfun3 <- function(gparm,fato,wx){
##     cc <- integrate(cons2,lower=0,upper=Inf,gparm=gparm,wx=wx,subdivisions=1000)
##     dem <- log(cc$value)
##     num <- sum(log(dgamma(fato,shape=gparm[1],scale=gparm[2])))
##     return(num-length(fato)*dem)
## }

## llkfun3(c(.01,.1),fato,wx)

# fitting weighted likelihood (using the true parameter values as starting points for estimation)
llkfun<-function(gparm,fato,wx){
   sum(log(dgamma(fato,shape=gparm[1],scale=gparm[2]))-log(sum(wx[,2]*dgamma(wx[,1],shape=gparm[1],scale=gparm[2]))))
}
llkfun2<-function(gparm, fato, wx){
  sum(1/wx[pmax(1,findInterval(fato,wx[,1])),2]*log(dgamma(fato, shape=gparm[1], scale=gparm[2])))
}

## llkfun4<-function(gparm, fato, wx){
##     sum((1/w(fato,wx))*log(dgamma(fato, shape=gparm[1], scale=gparm[2])))
## }



##fab3<-optim(par=c(2,20),llkfun3,fato=fato,wx=wx,control=list(parscale=c(1,1),fnscale=-1),lower=c(.01,.15),upper=c(50,40),method="L-BFGS-B")$par

## mean(fato)
## prod(fab)
## prod(fab3)

# fit wMLE
fab<-optim(par=c(2,20),llkfun,fato=fato,wx=wx,control=list(parscale=c(1,1),fnscale=-1),lower=c(.01,.1),upper=c(50,40),method="L-BFGS-B")$par
# fit wL
fab2<-optim(par=c(2,20),llkfun2,fato=fato,wx=wx,control=list(parscale=c(1,1),fnscale=-1),lower=c(.01,.1),upper=c(50,40),method="L-BFGS-B")$par

##fab4<-optim(par=c(2,20),llkfun4,fato=fato,wx=wx,control=list(parscale=c(1,1),fnscale=-1),lower=c(.01,.1),upper=c(50,40),method="L-BFGS-B")$par

xx <- seq(0,200,by=.1)

plot(xx,pgamma(xx,shape=a,scale=b), col=2, typ='l',xlab='Distance from turbine', ylab='Fraction of carcasses within r meters of turbine',lwd=2)
Twl<-wx[,2]*dgamma(wx[,1],shape=a,scale=b)
lines(wx[,1], cumsum(Twl/sum(Twl)), col=2)
# the distribution of simulated carcass distances
lines(sort(fatd),1/nfat*(1:nfat), lty=3, lwd=2)
# the empirical distribution of observed counts
lines(sort(fato),(1:length(fato))/length(fato), lty=3)
# the fitted distribution of counts based on RP sample
lines(xx,pgamma(xx,shape=fab[1],scale=fab[2]),col=4)
lines(xx,pgamma(xx,shape=fab2[1],scale=fab2[2]),col=3)
legend(x='bottomright',.2,legend=c(
  "true distribution",
  "true RP distribution",
  "empirical distribution",
  "empirical RP distribution",
  "fitted distribution (wD)",
  "fitted distribution (wL)"),
  col=c(2,2,1,1,4,3),lty=c(1,1,3,3,1,1),lwd=c(2,1,2,11,1))
dwp<-list()
dwp$true<-sum(wx[1:200,2]*diff(pgamma(wx[,1],shape=a, scale=b)))
dwp$wD<-sum(wx[1:200,2]*diff(pgamma(wx[,1],shape=fab[1], scale=fab[2])))
dwp$wL<-sum(wx[1:200,2]*diff(pgamma(wx[,1],shape=fab2[1], scale=fab2[2])))
dwp

### slow boat simulation:
a<-2.04459; b<-19.15065
nsim<-1000
simres<-data.frame(array(dim=c(nsim,8)))
names(simres)<-c('wD.a','wD.b','wD.dwp','wD.pdiff', 'wL.a','wL.b','wL.dwp','wL.pdiff')
nfat<-500 # total number of carcasses
simi<-634
for (simi in simi:nsim){
# NOTE: there is no guarantee that models will converge for given data, so sim may crash
# if it crashes, just restart where it left off (already built into the loop... simi:nsim)
  fatd<-rgamma(nfat,shape=a,scale=b)
##  ang<-2*runif(nfat)*pi # directions carcasses fly
##  ocarc<-inout(fatd*cbind(cos(ang),sin(ang)), RPpoly, bound=T)
  ocarc <- rbinom(nfat,size=1,prob=w(fatd,wx,0))==1
  print(sum(ocarc))
  fato<-fatd[ocarc] # the actual locations are not tracked; only the distance and whether or not they land on RP
  fab<-optim(par=c(2,20),llkfun,fato=fato,wx=wx,control=list(parscale=c(1,1),fnscale=-1),lower=c(.01,.1),upper=c(50,40),method="L-BFGS-B")$par
  fab2<-optim(par=c(2,20),llkfun2,fato=fato,wx=wx,control=list(parscale=c(1,1),fnscale=-1),lower=c(.01,.1),upper=c(50,40),method="L-BFGS-B")$par
  simres[simi,1:3]<-c(fab, sum(wx[1:200,2]*diff(pgamma(wx[,1],shape=fab[1], scale=fab[2]))))
  diffCDF<-pgamma(wx[,1],shape=fab[1], scale=fab[2])-pgamma(wx[,1],shape=a, scale=b)
  simres[simi,4]<-diffCDF[which(abs(diffCDF)==max(abs(diffCDF)))]
  simres[simi,5:7]<-c(fab2, sum(wx[1:200,2]*diff(pgamma(wx[,1],shape=fab2[1], scale=fab2[2]))))
  diffCDF<-pgamma(wx[,1],shape=fab2[1], scale=fab2[2])-pgamma(wx[,1],shape=a, scale=b)
  simres[simi,8]<-diffCDF[which(abs(diffCDF)==max(abs(diffCDF)))]
  if (simi%%100==0){cat(paste0(simi,'\n'));flush.console()}
}
simi

windows()
## two measures of model performance:
#1. how well does it estimate dwp?
boxplot(simres$wD.dwp,simres$wL.dwp, xlab='model', ylab='estimated dwp', names=c('wD', 'wL'))
axis(2)
axis(1, at=1:2, lab=c('wD', 'wL'))
box()
lines(par('usr')[1:2],rep(dwp$true,2), col=2)
legend(x='bottomleft', bty='n',legend='actual dwp', lty=1, col=2)
#2. how close does it match actual distribution?
boxplot(simres$wD.pdiff, simres$wL.pdiff, names=c('wD', 'wL'))
hist(abs(simres$wL.pdiff)-abs(simres$wD.pdiff))
sum(abs(simres$wL.pdiff)-abs(simres$wD.pdiff)>0)/nsim





