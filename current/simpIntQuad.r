

simpIntQuad <- structure( 
function(
#-------------------------------------------------------------
##############################################################
# DESCRIPTION: SIMPSON'S RULE QUADRATURE W/ INTERPOLATION BETWEEN SEARCHED-AREA BANDS TO GENERATE MORE DIVISIONS
# AUTHOR: Daniel Riser-Espinoza
# LAST UPDATED: 02-05-2015
# REQUIRED PACKAGES:
##############################
fun,
# function, to be evaluated using simpson's quadrature formula
minRad,
# numeric vector of input values 
maxRad,
# truncation radius
w,
# vector of weights = inverse proportion of area searched in a distance band, or probability of 
# carcass being found in ring_i, given it was anywhere in that ring_i; must include 0-m = 1.0
...
# other arguments to be passed to input function ( e.g. mu, sigma, shape, scale, etc. )
#-----------------------------------------------------------------------------------------
) {

# do linear interpolation on average area searched values; NB: area searched vector needs a 0 entry;
# interpolates 3 points between every one supplied
areaInterp <- approx( minRad : maxRad, w, xout = seq( minRad, maxRad, 
  by = maxRad / ( ( 4 * maxRad ) ) ) )
w <- areaInterp$y
# set new number of divisions
d <- length( w )

# calculate interval length
delta = ( maxRad - minRad ) / ( d - 1 )
# evaluate function at interval end points
endpts = areaInterp$y *  dtrunc( areaInterp$x, fun, a = minRad, b = maxRad, ... ) 
# separate and combine even and odd components as per Simpson's Rule
evens = 4 * sum( endpts[seq( 2, d - 2, by = 2 )] ) + endpts[maxRad]
odds = 2 * sum( endpts[seq( 3, d - 1, by = 2 )] ) + endpts[1]
# approximate integral
A = ( delta / 3 ) * sum( c( evens, odds ) )

return( A )

}

)
