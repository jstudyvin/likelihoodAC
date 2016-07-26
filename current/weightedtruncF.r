
# updated  truncated MLE function for a weighted distribution; If you supply weights as relative 
## probabilities of detection the function works as D. Dalthorp describes and fits a distribution
## with the understanding that the data are weighted in such a way as to distort the distribution
## If you supply weights as the inverse fo the relative probabilities of detection then the fitted
## function describes the observed distribution (it is w(x)f(x) / integral(w(y)f(y))) and 
## it needs to be "un-weighted" to get the area correction out.
## The function should be adaptable to any distribution 
# with out to 2 parameters
weightedtruncF <- structure( 
function(
#-------------------------------------------------------------
##############################################################
# DESCRIPTION: log likelihood of weighted truncated distribution
# AUTHOR: Daniel Riser-Espinoza
# LAST UPDATED: 5-21-2016 by PR
# REQUIRED PACKAGES: truncdist
##############################
theta,
# vector of parameters, shape and rate
fun,
# character; function of interest to be used in MLE fit of weighted truncated distribution
x,
# numeric vector of input values 
minRad,
# minimum radius; defaults to 0
maxRad,
# truncation radius
wX = NULL, # Weights associated with each carcass.  These are the relative probabilities 
# that a carcass within an annulus will enter the sample and each carcass needs a wX.  
# If both wX and wAnnulus are supplied
# wX is used by default.  If wX is not supplied, the weight is assumed to be the weight
# from the appropriate annulus in wAnnulus.  If neither is supplied, 
wAnnulus = NULL, #vector of inverse weights for each distance minRad:(maxRad - 1).  
# Should all be between 0 and 1. 
# vector of weights = inverse proportion of area searched in a distance band, or probability of 
# carcass being found in ring, given it was in that ring; must include 0-m = 1.0
inverseWeights = FALSE
#-------------------------------------------------------------
) {

##check that weights are supplied as relative probabilities:
if(!inverseWeights & any(wAnnulus > 1 | wX > 1 | wAnnulus ==0 | wX == 0)){
  stop('weightedtruncF: Weights must be between 0 and 1')
  }
 
##if  wAnnulus is supplied, check that it is the correct length
if(!is.null(wAnnulus) & length(wAnnulus) != length(minRad:maxRad)){
  stop('weightedtruncF: wAnnulus weight vector must have one weight for each annulus')
  }
  
##if  wX is supplied, check that it is the correct length
if(!is.null(wX) & length(wX) != length(x)){
  stop('weightedtruncF: wX weight vector must have one weight for each carcass')
  }   
  
##derive wX from wAnnulus if necessary
if(!is.null(wAnnulus) & is.null(wX)){
  #Associate the appropriate distance-based weight with each observation
  wX = wAnnulus[x - minRad + 1]
}

##derive wAnnulus from wX if necessary
if(is.null(wAnnulus) & !is.null(wX) ){
  #by distance mean of wX and interpolation if necessary
  wAnnulus = aggregate(wX, by = list(x), mean)
  wAnnulus = approx(x = wAnnulus$Group.1, y = wAnnulus$x, xout = minRad:maxRad, rule = 2)$y
}


# sample size
n <- length( x )

 arguments = list( x, fun, minRad, maxRad, theta[1], theta[2] )
 arguments = arguments[!is.na(arguments)]
# likelihood term 1: density log-likelihood term
llik1 <- -sum( log( do.call('dtrunc', arguments) ) )

# likelihood term 2: weighted log-likelihood term
# llik2 <- n * log( simpIntQuad( fun, lower = 0, upper = maxRad, w = w, theta[1], theta[2]  ) )
 arguments = list( fun, minRad, maxRad, w = wAnnulus, theta[1], theta[2] )
 arguments = arguments[!is.na(arguments)]
llik2 <- n *  log( do.call('simpIntQuad', arguments  ) )  
# sum likelihood terms
ll <- llik1 + llik2
# return( list(llik1, llik2) )

return(ll)

}

)
