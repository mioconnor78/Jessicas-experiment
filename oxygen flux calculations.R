### Oxygen data file for correcting for temperature dependence of oxygen exchange with the atmosphere.

ln(P) = 5.25 x ln(1-h/44.3)

#atmospheric pressure at our ponds:
h <- 0.3 #km above sea level
pressure.mod <- 5.25 x log(1-h/44.3) # = log(P) in atm at altitude (h, km) relative to standard partial pressre (PST) at 760 mm HG or 101.325 kpa at sea level.
c.star <- 8.6  #equilibrium oxygen concentration at standard pressure of 1 atm, (mg/L)

#equilibrium oxygen exhange
# T in Kelvin
P.wv <- function(T) exp((11.88571 - *3840.70/(T+273) - (216,961/(T + 274)^2))
theta <- fuction(T) 0.000975 - (T * 1.426 x 10^-5) + ((T^2) * 6.436 x 10^-8)

Cp <- c.star * P*[((1-P.wv/P)(1-theta*P))/((1-P.wv)(1-theta))]

O.sat <- function(DO) [100*DO]/Cp


### other examples:
C.star = function(Pt, p) c.760 * (Pt - p) / (760 - p) # C.star is teh oxygen solubility at 15C
# c*760 = saturation value at 760 mm HG = 10.08 mg/L
# Pt = barometric pressure
# p = vapour pressure of water (get this from a table, this is temperature dependent)

# so if 
Pt = 29.33
p = 12.79
c.760 = 10.08

C.star(Pt, p)

### another version, works well, is simplest: using this one.
C.star <- function(T) exp(7.7117 - 1.31403*log(T+45.93)) - 0.035
# the -0.035 is an approximate adjustment for elevation
C.star(T)

# So the change in oxygen resulting from a change in temperature, following this equation, tells me approximately how much o2 moved for physical reasons. I think I should just subtract this difference from the observed difference, right?

data$NPP2 <- (((data$dusk[1] - data$dawn1[1]) - (C.star(data2$temp.dusk)[1] - C.star(data2$temp.dawn1)[1]))*1000)/(32) 

o2.obs <- (data$dusk[1] - data$dawn1[1]) #this is how much oxygen changed in the water, as we saw it. some of that change was due to physical processes (o2.phys)
o2.phys <- (C.star(data2$temp.dusk)[1] - C.star(data2$temp.dawn1)[1]) #over that same time, there was this change in oxygen due to physical processes. negative number means less oxygen. So the NPP should be the observation plus the absolute value (or minus?) the physical processes. 
o2.npp <- o2.obs - o2.phys


## bring in full data file
data <- read.csv("./temporal_dataFEB12.csv")

## bring in data file with temperatures at each sampling time
o.data <- read.csv("./oxygen_temp_temporal.csv")

data2 <- merge(data, o.data, by.x = c("week", "Tank"), by.y = c("week", "Tank"))
data2 <- data2[,-(38:49)]


http://www.env.gov.bc.ca/wat/wq/BCguidelines/do/do-01.htm
