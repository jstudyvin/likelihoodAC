
##plots output from testDistv4.r

plotIt = function(
  x, #observed carcass locations
  fun, #name of underlying probability density function
  theta, #parameters associated with fun
  minRad, #turbine diameter
  maxRad, #maximum plot radius
  ##wAnnulus
  wAnnulus = NULL, #weights (as probabilities) for each annulus
  wX = NULL,  #weights for each carcass
  simN = 10000,
  rank = NULL, 
  field = NULL,
  figures = F){
  #simulate data
theta = theta[!is.na(theta)]
f = match.fun(paste0('d', fun))
pfun = match.fun(paste0('p', fun))
simX = simDensity(x, fun, theta, minRad, maxRad, wAnnulus = wAnnulus, wX = wX, simN = simN)
interp = simX$interp
simX = simX$simX
densX = density(x, from = minRad, to = maxRad)
densSimX = density(simX, from = minRad, to = maxRad)
##underlying density function
parametricArgs = minRad:(maxRad-1)
parametricArgs = list(parametricArgs, theta[1], theta[2], theta[3], theta[4], theta[5])
parametricArgs = parametricArgs[!is.na(parametricArgs)]
parametricDens = do.call(f, parametricArgs)
##Carcasses expected in cleared plots with radius = maxrad
carcsInClearedPlotArgs1 = list(maxRad, theta[1], theta[2], theta[3], theta[4], theta[5])
carcsInClearedPlotArgs1 = carcsInClearedPlotArgs1[!is.na(carcsInClearedPlotArgs1)]
carcsInClearedPlotArgs2 = list(minRad, theta[1], theta[2], theta[3], theta[4], theta[5])
carcsInClearedPlotArgs2 = carcsInClearedPlotArgs2[!is.na(carcsInClearedPlotArgs2)]
carcsInClearedPlot = 
  do.call(pfun, carcsInClearedPlotArgs1) - do.call(pfun, carcsInClearedPlotArgs2) 
##plot accoutrements
main = paste0('Distance distribution fit: truncated ', fun, ' density\n')
m3 = 4  #adjust figure margins
if(interp) {
  m3 = 5  #adjust figure margins
	main = paste0(main, 
    'Relative detection probabilities were interpolated between data points')
}
pred = paste0('Model predicts ', 100 * round(carcsInClearedPlot,2), 
  '% of carcasses within ', round(maxRad), ' m')
ylim = ylim = c(0, max(c(densX$y, densSimX$y), parametricDens))
##plot
par(mar = c(5, 4, m3, 2) + 0.1) #adjust figure margins
plot(densX, ylim = ylim, xlab = 'Distance from turbine', 
	ylab = 'Relative density', main = '')
lines(densSimX, lty = 'dashed')
lines(parametricDens, col = 'pink')
points(x = x, y = rep(ylim[2]/75, length(x)), pch = 124, col = adjustcolor('black', 0.15))
legend('topright', c(paste0('Observed density; n = ', length(x)), 
	paste0('Simulated density: n = ', length(simX)), 'Fitted density', 'Observed carcasses'), 
	lty = c('solid', 'dashed', 'solid', 'solid'), lwd = c(1, 1, 1, -1), pch = c(-1, -1, -1, 124), 
	col = c('black', 'black', 'pink', 'black'))
 mtext(pred, 3, 0.8)
 if(figures == TRUE) mtext(paste('Model rank is', rank, 'of', field), 1, 4)
 mtext(main, 3, m3-2.5)
}
