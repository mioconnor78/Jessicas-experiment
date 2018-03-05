## script for estimating equilibrium oxygen concentration at nonstandard pressure. Build based on equations at [add link]

## We're going to assume standard pressure for experiments at UBC, near sea level.

C.p # oxygen concentration at equilibrium, mg/L, at nonstandard pressure. go get the rest of the equation if you ever need it.

C.star # equilibrium concentration at standard pressure, mg/L

C.star <- function(t) exp(7.7117 - 1.31303*log(t + 45.93)) # for t in celcius. these constants must be in the mortimer references; find them.

C.star(25)
C.star(0)
