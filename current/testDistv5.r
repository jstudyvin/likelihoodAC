
########################
#####                       
########################

testDistv5 <- structure( 
function(
#-------------------------------------------------------------
##############################################################
# DESCRIPTION: test candidate distributions to find best fit, 
# via AIC, for distance data set; distributions are
# gamma
# gompertz
# log-logistic
# rayleigh
# truncated normal
# weibull

##adding distributions should be straightforward

# AUTHOR: DRE, PAR
# LAST UPDATED: 5-21-2016
# REQUIRES: VGAM, FAdist, truncdist,  weightedtruncF (sourced function)
################################################
# BEGIN FUNCTION ARGUMENTS:
x, 
# vector: input carcass distance data;
minRad,
# Scalar: minimum radius at which carcasses can occur (about 2m for the radius of a turbine)
maxRad, 
# Scalar: maximum distance at which carcasses can occur (largest plot radius)
wX = NULL,
# weights for individual carcasses.  takes precedence over wAnnulus. Default value is 1.
# see commentary in the function weightedtruncF (weightedTruncatedLogLikelihood.r)
wAnnulus = NULL,
# weights for individual annuli.  Ignored if wX is given.  Default value is 1
# see commentary in the function weightedtruncF (weightedTruncatedLogLikelihood.r)
writeOut = FALSE, 
# generate .csv file with model fit results?
project = NULL,
# project name
figures = FALSE , 
# should any figures of fits be generated?  options are 
  # TRUE (one figure per distribution)
  # FALSE (no figures)
  # 'best' (only the best fitting distribution is plotted)
simN = 10000,
# minimmum number of simulated data points for figures
inverseWeights = F ##If T, weights are 1/relative probability of detection a la Jared's
## weighted distribution method
# END FUNCTION ARGUMENTS

) { 
require(VGAM)  # needed for erf(x), the error function, and for Rayleigh dist.
require(FAdist) # needed for some distributions?
require( truncdist ) # needed for truncated data functions (dtrunc etc.)

 # add small buffer to plot radii so fatalities on the edge are included
 minRad = max(0, minRad - 0.01)
 maxRad = maxRad + 0.01
 
 # QA checks on data     
 x <- as.numeric(x) 
 
 if(any(x > maxRad ) || any(x < minRad)){
  stop('testDistv5: Carcasses are outside of plot boundaries')
}
  
 # set up output dataframe
 output <- data.frame( 
  distribution = c('rayleigh', 'gamma', 'norm', 'llog', 'gompertz', 'weibull'),   
  nParm = c(            1,        2,      2,       2,        2,         2),
  parameter1 = NA, parameter2 = NA,   aic = NA, convergence = NA  )
  
  starts = list( mean(x)/sqrt(pi/2),  #rayleigh starting parameter values
            c(1, 1 / mean(x)),  #gamma starting parameters
            c(mean(x), sd(x)),  #normal starting parameters
            c(1, mean(log( x ))) , #log-logistic starting parameters
            c(0.001, 0.001),  #gompertz starting parameters
            c(1,10)) #weibull starting parameters
            
# for(i in 1:nrow(output)){
 for(i in c(1,2,3,4,5,6)){
 # for(i in 1){

fit = tryCatch(nlminb(starts[[i]], objective = weightedtruncF, lower = c(0.001, 0.001),  
  x = x,  minRad = minRad,  maxRad = maxRad,   fun = output$distribution[i],   
    wX = wX,   wAnnulus = wAnnulus, inverseWeights = inverseWeights), error = function(e) 'noFit')
	

if(fit[1] != 'noFit'){
output$parameter1[i] = fit$par[1]
output$parameter2[i] = fit$par[2]
output$convergence[i] = fit$convergence
output$aic[i] = weightedtruncF(theta = fit$par, x = x,  minRad = minRad,  
  maxRad = maxRad, fun = output$distribution[i], wX = wX,  
  wAnnulus = wAnnulus) * 2 + 2 * output$nParm[i] 
  } else {
output$parameter1[i] = NA
output$parameter2[i] = NA
output$convergence[i] = 'noFit'
output$aic[i] = NA
}  #noFit
}  #cycle through distributions
  
##Write AIC table
output$deltaAIC = output$aic - min( output$aic )
output <- output[order( output$deltaAIC ),] 

if( writeOut == TRUE ) { 
  write.csv(output, paste( project, 'DistriubtionFitOutput.csv', sep = "_"),
  row.names = FALSE) 
  } #write out
  
if(figures == TRUE | figures == 'best'){
  for(i in 1:nrow(output)){
  jpeg(paste0(project, '_', output$distribution[i], '.jpg'), 
  height = 6, width = 8, res = 300, units = 'in')
  plotIt(x = x, fun = output$distribution[i], 
    theta = c(output$parameter1[i], output$parameter2[i]), minRad = minRad, maxRad = maxRad, 
    wAnnulus = wAnnulus, wX = wX, simN = simN, rank = i, field = nrow(output), figures = figures)
  dev.off()
  if(figures == 'best') break()
  }
}  #figures
return(output)
}  #function
)  #structure

