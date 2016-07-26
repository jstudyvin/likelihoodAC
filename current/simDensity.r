
##simulate carcass observations given a density distribution and weights for each annulus

##function to simulate a truncated, weighted probability density
simDensity = function(x, fun, theta, minRad, maxRad, wAnnulus = NULL, wX = NULL, simN = 10000){
##If weights are not given for the annuli, interpolate them from the carcass distances	
  interp = F
	if(is.null(wAnnulus)){
		#Throw a flag to indicate interpolation
		interp = T
		wAnnulus = approx(x = x,  y = wX, xout = minRad:maxRad, rule = 2)$y
  } # if(is.null(wAnnulus))
	simX = NULL
  rtruncArgs = list(mean(wAnnulus)^-1 * simN, fun, minRad, maxRad-1, theta[1], theta[2],
    theta[3], theta[4], theta[5])
  rtruncArgs = rtruncArgs[!is.na(rtruncArgs)]
	#Iterate the next few steps until there are enough samples
	while(length(simX) < simN){
	#Generate according to the underlying distribution (truncated)
	simX0 = round(do.call(rtrunc, rtruncArgs))
	#Weights are probability in a binomial trial for each carcass
  tProbs =  wAnnulus[simX0 - minRad + 1]
	test = rbinom(rep(1, length(simX0)), rep(1, length(simX0)), tProbs)
	simX = c(simX, simX0[test == 1])
	}
	return(list(simX = simX, interp = interp))
	}
