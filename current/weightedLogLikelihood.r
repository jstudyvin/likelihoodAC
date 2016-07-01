###Density for the log-likelihood of a truncated distribution with weighted data
## Author Rabie
## last modified 5-21-2016 by PR

###This function is a weighted log likelihood that does not perform as well as a weighted likelihood.
###use weightedtruncF instead.

wtruncLL = function(
  theta, #parameters associated with the distribution function
  
  x, #vector of distances assumed to be integer and within minRad:maxRad (inclusive)
  
  minRad, #Lower truncation bound on search plot(eg radius of turbine).  If there are
  # multiple plot types contributing to the sample this should be the minimum possible
  # distance from the center of the turbine at which a carcass could occur.
  # Variation in minRad among plots should be addressed in the 
  # weights (wAnnulus or wX, below)

  maxRad, #Upper truncation bound on search plot (eg max plot radius)  If there are
  # multiple plot types contributing to the sample this should be the maximum possible
  # distance from the center of the turbine at which a carcass could occur. 
  # Variation in maxRad among plots should be addressed in the 
  # weights (wAnnulus or wX, below)

  fun, #underlying probability distribution.  Needs to be a probability distribution
  #  for which a dfun() function is defined
  #  wtruncLL can accommodate probability distributions with up to 5 parameters, which 
  #  sounds excessive unless you are considering mixture distributions.
  
# only one of wX and wAnnulus needs to be supplied.  Read below.  
  
  wX = NULL, # Weights associated with each carcass.  These are the relative probabilities 
  # that a carcass within an annulus will enter the sample and each carcass needs a wX.  
  # If both wX and wAnnulus are supplied
  # wX is used by default.  If wX is not supplied, the weight is assumed to be the weight
  # from the appropriate annulus in wAnnulus.  If neither is supplied, 
  
  wAnnulus = NULL #vector of inverse weights for each distance minRad:(maxRad - 1).  
  # Should all be between 0 and 1. 
  
  #####################################################################################
  #####################################################################################
  ##WEIGHTS:  Usually these will be proportion of area searched within the annulus
  #  (wAnnulus), or associated with the carcass (wX)  but if
  # multiple plot types are pooled, then other factors that vary as a function of distance
  # should be considered, including variable SEEF, variable CR, variable underlying 
  # CR rates & plots of differing sizes/configurations. 
  # For complicated data sets it may be simplst to associate one weight with
  # each carcass (i.e. supply wX) rather than try to calculate a defensible 
  # weighted average within annuli (i.e. supply wAnnulus).
  
  # If wX is supplied it will be used and wAnnulus will be ignored.  If only wAnnulus
  # is supplied, wAnnulus will be used.  If neither is supplied, all weights are 
  # assumed to be 1.0 (i.e. assumption that all plots are circular, completely searched, 
  # and with a common radius)   
  ){ 
  require(truncdist)
  require(VGAM) #some of the density functions
  require(FAdist) #log logistic
  
##check that weights are supplied as relative probabilities:
if(any(wAnnulus > 1 | wX > 1 | wAnnulus ==0 | wX == 0)){
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


#if neither wX nor wAnnulus is supplied all weights are 1, otherwise invert the weights
if(!is.null(wX)) {
  wX = 1 / wX
  } else {
  wX = rep(1, length(x))
}

  ##This formulation permits the use of distributions with up to 5 parameters
  arguments = list(x, fun, minRad, maxRad, theta[1], theta[2], theta[3], theta[4], theta[5])
  arguments = arguments[!is.na(arguments)]

  ll = - sum(wX * log(do.call('dtrunc', arguments)))
  
return( ll )
}
