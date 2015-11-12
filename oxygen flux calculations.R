### Oxygen data file for correcting for temperature dependence of oxygen exchange with the atmosphere.

ln(P) = 5.25 x ln(1-h/44.3)

#atmospheric pressure at our ponds:
h <- 0.3 #km above sea level
pressure.mod <- 5.25 x log(1-h/44.3) # = log(P) in atm at altitude (h, km) relative to standard partial pressre (PST) at 760 mm HG or 101.325 kpa at sea level.

#equilibrium oxygen exhange
T <- 273 # T in Kelvin
P.wv <- function(T) exp((11.88571 - *3840.70/T) - (216,961/T^2))
theta <- fuction(tc) 0.000975 - (tc * 1.426 x 10^-5) + ((Tc^2) * 6.436 x 10^-8)

Cp <- c.star * P*[((1-P.wv/P)(1-theta*P))/((1-P.wv)(1-theta))]

O.sat <- function(DO) [100*DO]/Cp