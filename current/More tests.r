
library(truncdist)
options(stringsAsFactors = F)
setwd('Y:/_R code/weightedMLE')

source('testDistv5.r')
source('testDistv6.r')
source('simpIntQuad.r')
source('simDensity.r')
source('plotIt.r')
source('weightedtruncF.r')
source('weightedLogLikelihood.r')
##source('testDistWeightedLikelihood.r')
##source('calcAreaCorrection.r')



xbar = 40
sdx = 18
godSprinkles = round(rtrunc(100, 'norm', a = 1, b = 99, mean = xbar, sd = sdx))

hist(godSprinkles)

r = 1:100
w = dexp(r, 1/20) / dexp(1, 1/20)
w = dgamma(r, 3, 1/20) / max(dgamma(r, 3, 1/20))

probs = w[godSprinkles]
found =  sapply(probs, function(p) rbinom(1, 1, p) == 1)
x = godSprinkles[found]

# plot(density(x))
# lines(dgamma(r, 3, 1/20))


##DRE & DD method
wtdDist = testDistv5(x, minRad = 1, maxRad = 100, wX = NULL, wAnnulus = w)
##PR method
wtdLik = testDistWeightedLikelihood(x, minRad = 1, maxRad = 100, wX = NULL, wAnnulus = w)
##JS method
invWtdDist =  testDistv6(x, minRad = 1, maxRad = 100, wX = NULL, wAnnulus = 1/w, inverseWeights = T)

# plot(invWtdDist$fStarTable[,3])
# lines(density(x))


minRad = 1
 maxRad = 100
 wX = NULL
 wAnnulus = 1/w
 inverseWeights = T

calcAreaCorrection(fun = wtdDist[1, 1], propSrchd = matrix(w, nrow = 1), plotRad = max(r),
	bandwidth = 1, maxRad = max(r), wtdDist[1, 3], wtdDist[1, 4])

calcAreaCorrection(fun = wtdLik[1, 1], propSrchd = matrix(w, nrow = 1), plotRad = max(r),
	bandwidth = 1, maxRad = max(r), wtdLik[1, 3], wtdLik[1, 4])


sum(invWtdDist$fStarTable[,'prop'] / invWtdDist$fStarTable[, 'inputWeight'])

calcAreaCorrection('norm', propSrchd = matrix(w, nrow = 1), plotRad = max(r),
	bandwidth = 1, maxRad = max(r), xbar, sdx)

