### Garzke et al temperature experiment
### this is the current working code file
### MO made a new file on Aug29 from June file when I decided to take week out as a fixed effect, and instead model autocorrelation. 

### load libraries
#library(qpcR)
library(nlme)
library(MuMIn)
library(plyr)
library(tidyverse) 
library(broom)

### set working directory and load data
data <- read.csv("./data/temporal_dataFEB12.csv")
dim(data)
head(data)
tail(data)
data <- data[-(241:255),]


### load temperature data
tdata <- read.csv("./data/avgtemps.csv")
head(tdata)
temp.data <- ddply(tdata, .(Week, Tank), summarise, mean(Temperature)) 
head(temp.data)
names(temp.data) <- c('Week', 'tank', 'wklyTemp')

tank.means <- ddply(tdata, .(Tank), summarise, mean(Temperature)) 
head(tank.means)
names(tank.means) <- c("Tank", "TankTemp")

## bring in data file with temperatures at each sampling time
o.data <- read.csv("./data/oxygen_temp_temporal.csv")
o.data <- o.data[,-(4:14)]
o.data <- o.data[,-2]
head(o.data)
dim(o.data)
data2 <- merge(data, o.data, by.x = c("week", "Tank"), by.y = c("week", "Tank"))
data <- data2

data3 <- merge(data, tank.means, by.x = "Tank", by.y = "Tank") #tank.means from temperature analysis file

head(data3)
data <- data3

### load libraries, define variables and add columns
k <- 8.617342*10^-5  # eV/K
data$invTi <-  1/((data$average.temp + 273)*k) # average temp of the week
data$invTT <-  1/((data$TankTemp + 273)*k) # average temp of the tank over all weeks
data$invTi <- as.numeric(as.character(data$invTi))
data$invTT <- as.numeric(as.character(data$invTT))

## some plots for the two temperature terms:
plot((data$invTi - data$invTT) ~ data$week)
plot((data$invTi - data$invTT) ~ data$invTT)
plot((data$invTi) ~ data$invTT)
plot((data$invTi - data$invTT) ~ data$Tank)

### estimate biomass from chlorophyll concentration
data$PP.biomass <- (data$chla*55) #chla (ug/L)* 55 C in PP / 1 chla = ugPPC/L
data$ZP.carbon1 <- ifelse(data$trophic.level=='P',  0, data$zoo.ug.carbon.liter) # for adding
data$total.carbon <- data$PP.biomass + data$ZP.carbon1 #I'm assuming 0 for ZP in P treatments here.

### estimating NPP and ER from the raw data: 
## correct of water-air oxygen flux
### adjusting for effects of temperature on physical exchange
C.star <- function(T) exp(7.7117 - 1.31403*log(T+45.93)) - 0.035
# the -0.035 is an approximate adjustment for elevation
C.star(T) # yields the oxygen concentration expected at a given temperature (T in C).


#calculate NPP and ER (hourly), in terms of umol O2 / l / hr, following yvon durochers 2010.
# O2 has molar mass of 32g/mol. so 1 umol = 32 ug. so take ug/32
data$NPP2 <- (((data$dusk - data$dawn1) - (C.star(data$temp.dusk) - C.star(data$temp.dawn1)))*1000)/(32)  # oxygen produced umol / L /day, net all respiration. raw data is mg/L. Subtract o2 water/atm flux due to change in temperature 
data$ER2 <- -(24/data$hours2)*(((data$dawn2 - data$dusk) - (C.star(data$temp.dawn2) - C.star(data$temp.dusk)))*1000)/(32)  # amount of oxygen consumed per day via respiaration. negative to get the change in oxygen umol / L /day; oxygen used in the dark and daylight. MeanER can be greater than meanNPP, because NPP reflects ER already.
data$GPP <- data$NPP2+(data$ER2/24)*data$hours1 # daily oxygen production (NPP2) + estimated daytime community respiration (daily R / 24 * hours daylight)
data$NEM <- data$ER2/data$GPP  # following Yvon Durochers 2010. NEM > 1 means the system is respiring more than it's fixing per day. This does not need to be logged.
data$NPP.mass <- data$NPP2 / (data$PP.biomass)  # NPP on ummol 02/L/day/ugCPP
data$ER.mass <- data$ER2/(data$total.carbon) # ER on ummol 02/L/day/ugTPP

data <- data[data$week >= '4',]


### data prep complete

### ANALYSIS
## trying the method suggested by Van de pol and Wright 2009
# center by within-tank temperature

data1 <- data

## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invTi - invTT), with Tank as a random intercept effect.  
modNPP0 <- lme(log(NPP2) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modNPP1 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")
modNPP2 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)  
modNPP4 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP5 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP6 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP7 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modNPP0, modNPP2, modNPP4,modNPP1, modNPP5, modNPP6, modNPP7)

## Best model: create fitted values to use later for plotting
intervals(modNPP7r, which = "fixed")
modNPP7r <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
#predw <- predict(modNPP7, level = 0:1)
## alternatively (yes, this works): 
mod.coefs <- augment(modNPP7r, effect = "random") # puts fitted values back in the original dataset.

### SOME BASIC PLOTS
#data1 <- data[(data$NPP2 >= 0.5),] # three negative values and one very small value now, not sure what to do about them.
hist(data[(data$NPP2 >= 0.5),]$NPP2)
hist(log(data1$NPP2))

plot(log(data1$NPP2)~data1$Tank, pch = 19, col = data1$trophic.level)
plot(log(data1$NPP2)~data1$week, pch = 19, col = data1$trophic.level)
plot(log(data1$NPP2)~I(data1$invTT-mean(data1$invTT)), pch = data1$Tank, col = data1$Tank, ylim = c(0,6))
plot(log(data1$NPP2)~I(data1$invTi-(data1$invTT)), pch = data1$Tank, col = data1$Tank, ylim = c(0,6))
abline(3.5435395, -0.6159088, lwd = 2, col = 1)
plot(log(data1$NPP2)~data1$invTi, pch = data1$Tank, col = data1$Tank)

### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

NPP.plot <- ggplot(data = data1, aes(x = invTi, y = log(NPP2))) + 
  theme_bw() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  xlab("Temperature 1/kTi") +
  ylab("ln(NPPi)")

## PLOT 1: Individual regression lines fitted within and among groups
NPP.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black')

## PLOT 2: Use fitted lines from the model. Create a dataframe with those model output coefficients... or add the predictions of the model to the original dataset and plot those here, then fit lines to those using linear regressions
NPP.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black')

## PLOT 3: Now try to plot the model results directly, as lines. Following van de pol and wright, we can plot all this on one temperature axis, with one slope for between group change, and another for within group change. So, we just have to figure out what are those coefficients...
# B0 = intercept varies within group (fixed + random)
# B1 = effect of temperature within groups (fixed) reflects main and interactive effect with the among-groups term. I'll plot it with the interaction term in the within group slopes
# B2 = the among groups term

## define NPP.func for effect of weekly temperature on NPP, the linear model relating temperature to NPP, where T is mean tank temperature and m is the tank temp in week i  

NPP.func <- function(Tw) fixef(modNPP7)[1] + (fixef(modNPP7)[3])*Tw 

# i thought i needed to include coefficients from all terms here, but this just isn't working
z <- 0.5 #invTi - invTT for each tank; i'm not sure if this should be the average devation? # ok if we make z = 0 the line is close. ideally, we set z as a function to estimate the mean deviation from the average for each tank, and then use that in the formula. come back to this.
#z <- function(Ti) { (Ti - mean(Ti))}
NPP.func2 <- function(x, Ti) { (fixef(modNPP7r)[1] - fixef(modNPP7r)[3]*mean(mod.coefs$invTT) - fixef(modNPP7r)[4]*(z)*mean(mod.coefs$invTT)) + (fixef(modNPP7r)[3] + fixef(modNPP7r)[4]*(z))*x }

yvals <- NPP.func2(mod.coefs$invTT)

## some attempts to write the z function; didn't get this working
z <- function (Ti) {
  diffs <- Ti - mean(Ti)
  return(data.frame(diffs))
}

devs <- data1 %>%
  group_by(Tank) %>%
  do(devs = z(.$invTi)) 

devs[[2]]

## this plot works if z = 0 or is small.
NPP.plot +
  #geom_smooth(data = mod.coefs, aes(x = (invTi), y = log(NPP2), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE, lwd = 0.75) +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvals, ymin = yvals - 0.3, ymax = yvals + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvals), lwd = 2) +
  geom_text(label = "B2i = -0.62", x = 40.0, y = 4.5) + 
  geom_abline(slope = (fixef(modNPP7r)[3] + fixef(modNPP7r)[4]*(z)), intercept = (fixef(modNPP7r)[1] - fixef(modNPP7r)[3]*mean(mod.coefs$invTT) - fixef(modNPP7r)[4]*(z)*mean(mod.coefs$invTT)), colour = "red") # this plots the model coefs, use it to check that the above method worked (and it did). this red line should be exactly on top of the black line.

ggsave("NPPplot.png", device = "png")

################################################
## Does mass-specific NPP vary with temperature?  
################################################
### NPP in umol/L/day per ugC


## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invTi - invTT), with Tank as a random intercept effect.  
modNPPm0 <- lme(log(NPP.mass) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modNPPm1 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")
modNPPm2 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)  
modNPPm4 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm5 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm6 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm7 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modNPPm0, modNPPm2, modNPPm4, modNPPm1, modNPPm5, modNPPm6, modNPPm7)

## Best model: create fitted values to use later for plotting
modNPPm6r <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
intervals(modNPPm6r, which = "fixed") 
mod.coefs <- augment(modNPPm6r, effect = "random") # puts fitted values back in the original dataset.

## figures  
### NPP in umol/L/day per ugC
data1 <- data[(data$NPP2 >= 0.5),]
hist(data1$NPP.mass)
hist(log(data1$NPP.mass))
# hist(log(data$NPP.mass+.001)) # no reason to add a number, there are no -inf values.

plot(log(data$NPP.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$NPP.mass)~I(data1$invTT-mean(data1$invTT)), pch = 19, col = data1$trophic.level)
plot(log(data1$NPP.mass)~I(data1$invT-(data1$invTT)), pch = 19, col = data1$trophic.level)
abline(-2.62, -1.83, lwd = 2, col = 1)
abline((-2.62+1.53), (-1.83-1.40), lwd = 2, col = 2)
abline((-2.62+0.13), (-1.83-0.05), lwd = 2, col = 3)


### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

NPPm.plot <- ggplot(data = data1, aes(x = invTi, y = log(NPP.mass))) + 
  theme_bw() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  xlab("Temperature 1/kTi") +
  ylab("ln(NPPi)")

## PLOT 1: Individual regression lines fitted within and among groups
NPPm.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black')

## PLOT 2: Fit lines to fitted data from the model [close here - go get the coefficients as estimated for NPP, and write them out here with a z value. should plot up well.]
NPPm.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black')

## PLOT 3: Plot lines from model
# after some algebra, i can isolate the slope for data modeled as centered
NPPmPP.func <- function(x) {(fixef(modNPPm6r)[1] - fixef(modNPPm6r)[3]*mean(data1$invTT)) + (fixef(modNPPm6r)[3] - fixef(modNPPm6r)[2] - fixef(modNPPm6r)[6]*x + fixef(modNPPm6r)[6]*mean(data1$invTT))*x} #x = invTT
# use this function to compute yvals for plotting.
yvals <- NPPmPP.func(data1$invTT)

NPP.plot +
  #geom_smooth(data = mod.coefs, aes(x = (invTi), y = log(NPP2), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE, lwd = 0.75) +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvals, ymin = yvals - 0.3, ymax = yvals + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvals), lwd = 2)
#stat_function(data = mod.coefs, aes(x = (mod.coefs$invTT)), fun = NPP.func2, geom = 'line')

NPPm.plot +
  #geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = interaction(Tank, trophic.level), color = trophic.level), method = "lm", se = FALSE) +
  geom_line(aes(x = (data1$invTT), y = yvals))
  
  stat_function(data = mod.coefs, aes(x = invTT), fun = NPPmPP.func, geom = 'line')
  
  geom_abline(slope = (fixef(modNPPm6r)[3] + fixef(modNPPm6r)[2]*mean(mod.coefs$invTT)), intercept = (fixef(modNPPm6r)[1]- fixef(modNPPm6r)[3]*(mean(mod.coefs$invTT))), colour = "pink") +
  geom_abline(slope = (fixef(modNPPm6r)[3] + fixef(modNPPm6r)[7] + fixef(modNPPm6r)[2]*mean(mod.coefs$invTT)), intercept = (fixef(modNPPm6r)[1] + fixef(modNPPm6r)[4] - fixef(modNPPm6r)[3]*(mean(mod.coefs$invTT)) - fixef(modNPPm6r)[4]*(mean(mod.coefs$invTT))), colour = "seagreen") + 
  geom_abline(slope = (fixef(modNPPm6r)[3] + fixef(modNPPm6r)[8] + fixef(modNPPm6r)[2]*mean(mod.coefs$invTT)), intercept = (fixef(modNPPm6r)[1] + fixef(modNPPm6r)[5] - fixef(modNPPm6r)[3]*(mean(mod.coefs$invTT)) - fixef(modNPPm6r)[5]*(mean(mod.coefs$invTT))), colour = "blue")
  
  
    geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = interaction(Tank, trophic.level), color = trophic.level), method = "lm", se = FALSE) +

  geom_abline(slope = fixef(modER2r)[3], intercept = (fixef(modER2r)[1])- fixef(modER2r)[3]*(mean(mod.coefs$invTT))) +
  geom_abline(slope = fixef(modER2r)[3], intercept = (fixef(modER2r)[1] + fixef(modER2r)[4] - fixef(modER2r)[3]*mean(mod.coefs$invTT))) +
  geom_abline(slope = fixef(modER2r)[3], intercept = (fixef(modER2r)[1] + fixef(modER2r)[5] - fixef(modER2r)[3]*mean(mod.coefs$invTT)))


geom_abline(slope = fixef(modNPPm6r)[3], intercept = (fixef(modNPPm6r)[1]- fixef(modNPPm6r)[3]*(mean(mod.coefs$invTT))), colour = "red")

ggsave("ERplot.png", device = "png")

######################################
## Does ER vary with temperature?  
#####################################
hist(data$ER2)
data1 <- data[(data$ER2 >= 0),]
hist(log(data1$ER2))

#analysis
modER0<-lme(log(ER2) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modER1<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modER2<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER4<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER5 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)
modER6 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER7 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modER0, modER1, modER2, modER4, modER5, modER6, modER7)

## Best model: create fitted values to use later for plotting 
modER2r<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
summary(modER2r)
intervals(modER2r, which = "fixed")
mod.coefs <- augment(modER2r, effect = "random") 

### SOME BASIC PLOTS
## figures 
plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$ER2)~I(data1$invTT-mean(data1$invTT)), pch = 19, col = data1$trophic.level, ylim = c(0,6))
abline(4.27, -0.43, lwd = 2, col = 1)
abline((4.27+0.35), (-0.43), lwd = 2, col = 2)
abline((3.88-0.04), (-0.43), lwd = 2, col = 3)

### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs
ER.plot <- ggplot(data = data1, aes(x = invTi, y = log(ER2))) + 
  theme_bw() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  xlab("Temperature 1/kTi") +
  ylab("ln(ERi)")

ER.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23) +
  geom_smooth(method = "lm", se = FALSE, aes(group = trophic.level), formula = y ~ x, color = 'black')

#1. the linear model for within groups variation
#NPPwg.fun <- function(invTi) 3.5435395 + (0.028891)*(invTi - mean(invTi))
#NPPwg.fun <- function(invTi) 3.5435395 + (0.028891)*(invTi) - (0.028891)*(mean(invTi))

#2. the linear model for among groups variation
#NPP.func <- function(invTT) 3.5435395 + (-0.6159088)*(invTT - mean(invTT))

summary(modER2)
#within groups, for each TG (different intercept)
#1. ERwg.fun <- function(invTi) fixef(modER2r)[1] + fixef(modER2r)[2]*(invTi - mean(invTi))
#   ERwg.fun <- function(invTi) fixef(modER2r)[1] + fixef(modER2r)[2]*(invTi) - fixef(modER2r)[2]*(mean(invTi))

#2. ER.func <- function(invTT) {fixef(modER2r)[1] + fixef(modER2r)[3]*(invTT - mean(invTT))}
    ER.func <- function(invTT) {fixef(modER2r)[1] + fixef(modER2r)[3]*(invTT) - fixef(modER2r)[3]*(mean(invTT))} # for trophic level 1

yvals <- ER.func(mod.coefs$invTT)

ER.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = interaction(Tank, trophic.level), color = trophic.level), method = "lm", se = FALSE) +
  #geom_ribbon(aes(x = (mod.coefs$invTT), y = yvals, ymin = yvals - 0.3, ymax = yvals + 0.3), fill = "grey70", alpha = 0.6) +
  geom_abline(slope = fixef(modER2r)[3], intercept = (fixef(modER2r)[1])- fixef(modER2r)[3]*(mean(mod.coefs$invTT))) +
  geom_abline(slope = fixef(modER2r)[3], intercept = (fixef(modER2r)[1] + fixef(modER2r)[4] - fixef(modER2r)[3]*mean(mod.coefs$invTT))) +
  geom_abline(slope = fixef(modER2r)[3], intercept = (fixef(modER2r)[1] + fixef(modER2r)[5] - fixef(modER2r)[3]*mean(mod.coefs$invTT)))
  
ggsave("ERplot.png", device = "png")




##################################################
## Does mass specific ER vary with temperature? 
##################################################
## figures 
data1 <- data[(data$ER2 >= 0),]
hist(data1$ER.mass)
hist(log(data1$ER.mass))

plot(log(data$ER.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~I(data$invTT-mean(data$invTT)), pch = 19, col = data$trophic.level)
abline(-1.89, -1.49, lwd = 2, col = 1)
abline((-1.89+1.54), (-1.49-0.74), lwd = 2, col = 2)
abline((-1.89+0.21), (-1.49+0.43), lwd = 2, col = 3)

## analysis
modERm0<-lme(log(ER.mass) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modERm1<-lme(log(ER.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modERm2<-lme(log(ER.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modERm4<-lme(log(ER.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modERm5 <- lme(log(ER.mass) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

model.sel(modERm0, modERm1, modERm2, modERm4, modERm5)

# for model fitting: 
modERm4r <- lme(log(ER.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 

modERm5r <- lme(log(ER.mass) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit)
summary(modERm4r)

m.avg <- model.avg(modERm4r, modERm5r)
summary(m.avg)

### plotting

Tb <- function(x) -0.61-1.03*x
x <- seq(37, 42, 0.001)
ggplot(data = mod.coefs, aes(x = invT, y = log(NPP2))) + 
  theme_minimal() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  geom_abline(intercept = 3.54, slope = Tb(x)) +
  #geom_ribbon(aes(x = x, ymin = Tb(x) ), fill = "grey70")
  geom_smooth(data = mod.coefs, aes(x = invT, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black') +
  xlab("Temperature 1/kT") +
  ylab("ln(NPP)")


##### BIOMASS RESULTS #####

#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))

plot(log(data$PP.biomass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~I(data$invTT-mean(data$invTT)), pch = 19, col = data$trophic.level)

abline(6.11, 0.87, col = 1, lwd = 2)
abline((6.11 -1.23), (.87 + 1.47), col = 2, lwd = 2)
abline((6.11 - 0.01), (0.87 + 0.25), col = 3, lwd = 2)

## analysis
modPP0<-lme(log(PP.biomass) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modPP1<-lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modPP2<-lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modPP4<-lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modPP5 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)
modPP6 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modPP7 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modPP0, modPP1, modPP2, modPP4, modPP5, modPP6, modPP7) 

## Best model: create fitted values to use later for plotting  
modPP6r <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
summary(modPP6r)
intervals(modPP6r, which = "fixed")
mod.coefs <- augment(modPP6r, effect = "random") 


### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

PP.plot <- ggplot(data = data1, aes(x = invTi, y = log(PP.biomass))) + 
  theme_bw() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  xlab("Temperature 1/kTi") +
  ylab("ln(PPi)")

## PLOT 1: Individual regression lines fitted within and among groups
PP.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black')

## PLOT 2: Use fitted lines from the model. 
PP.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black')

## PLOT 3: 
## define PP.func for effect of weekly temperature on PP, the linear model relating temperature to PP, where T is mean tank temperature and m is the tank temp in week i  

## next steps: follow template of ER model or NPP model; figure out the terms here.

# after some algebra, i can isolate the slope for data modeled as centered: UPDATE THIS FOR PP
NPP.func2 <- function(x) {fixef(modNPP7)[1] - fixef(modNPP7)[3]*mean(x) + (fixef(modNPP7)[3])*(x)} #x = invTT
# use this function to compute yvals for plotting.
yvals <- PP.func(mod.coefs$invTT)

## FOR PLOT, BORROW ABLINE APPROACH USED IN ER PLOT
NPP.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvals, ymin = yvals - 0.3, ymax = yvals + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvals), lwd = 2)


ggsave("PPplot.png", device = "png")






## Does zooplankton carbon vary with temperature?  ### leaving this for now. 
## figures 
dataz <- data[(which(data$trophic.level != 'P')),] # because we don't have observations for the P treatments (we're just assuming they are 0) I don't htink we should analyze those tanks here.
hist(dataz$zoo.ug.carbon.liter)
hist(log(dataz$zoo.ug.carbon.liter))
hist(log(dataz$zoo.ug.carbon.liter+1))
plot(log(dataz$zoo.ug.carbon.liter)~dataz$Tank, pch = 19, col = data$trophic.level)

plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$Tank, pch = 19, col = dataz$trophic.level)
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$week, pch = 19, col = dataz$trophic.level)
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$invT, pch = 19, col = dataz$trophic.level)

## analysis
## determine need for random effects in the full model: 
modZa<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modZb<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 
anova(modZa, modZb)

modTCc<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modTCc, modTCb)

modTC0<-lme(log(total.carbon) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modTC1<-lme(log(total.carbon)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modTC2<-lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modTC4<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

model.sel(modTC0, modTC1, modTC2, modTC4)

# for model fitting: 
modTC2r <- lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

summary(modTC2r)




### TOTAL CARBON #### 
## Does total biomass vary with temperature and trophic structure?  
## have to analyze this for weeks 4 and later, because we have no ZP data from weeks 2 and 3. 
## figures 
hist(data$total.carbon)
hist(log(data$total.carbon))

plot(log(data$total.carbon)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$total.carbon)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$total.carbon)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modTCa<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTCb<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 
anova(modTCa, modTCb)

modTC0<-lme(log(total.carbon) ~ 1, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTC1<-lme(log(total.carbon)~1+I(invT-mean(invT)), random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTC2<-lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTC4<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  

model.sel(modTC0, modTC1, modTC2, modTC4)
anova(modTC4, modTC2)
anova(modTC1, modTC2)
anova(modTC1, modTC0)

# for model fitting: 
modTC4r <- lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="REML", na.action=na.omit)  

summary(modTC4r)


