### MO made a new file on Aug29 from June file when I decided to take week out as a fixed effect, and instead model autocorrelation. 

### revised this on Nov 13 to include correction fo oxygen flux

### load libraries
#library(qpcR)
library(nlme)
library(MuMIn)

### set working directory and load data
data <- read.csv("./data/temporal_dataFEB12.csv")
dim(data)
head(data)
tail(data)
data <- data[-(241:255),]

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
data$invT <-  1/((data$average.temp + 273)*k)
data$invTT <-  1/((data$TankTemp + 273)*k)
data$invT <- as.numeric(as.character(data$invT))
data$invTT <- as.numeric(as.character(data$invTT))

## some plots for the two temperature terms:
plot((data1$invT - data1$invTT) ~ data1$week)
plot((data1$invT - data1$invTT) ~ data1$invTT)
plot((data1$invT) ~ data1$invTT)
plot((data1$invT - data1$invTT) ~ data1$Tank)


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

## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invT - invTT), with Tank as a random intercept effect.  
modNPP0 <- lme(log(NPP2) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modNPP1 <- lme(log(NPP2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")
modNPP2 <- lme(log(NPP2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)  
modNPP4 <- lme(log(NPP2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP5 <- lme(log(NPP2) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP6 <- lme(log(NPP2) ~ 1 + I(invT - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modNPP0, modNPP2, modNPP4,modNPP1, modNPP5, modNPP6)

## next step: Model average modNPP1 and modNPP6. Plot and interpret that.

## We need to be sure we don't have inflated degrees of freedom on the intercept and slope (invT) parameters. I think we still do:

anova(modNPP1)
numDF denDF  F-value p-value
(Intercept)                1   149 5186.790  <.0001
I(invT - invTT)            1   149    3.028  0.0839
I(invTT - mean(invTT))     1    28   43.847  <.0001

intervals(modNPP1, which = "fixed")

# so the inference for the invTT term is ok here. The invTT term is the important one.

## Does NPP vary with temperature?  
## figures on invT
data1 <- data[(data$NPP2 >= 0.5),] # three negative values and one very small value now, not sure what to do about them.
hist(data[(data$NPP2 >= 0.5),]$NPP2)
hist(log(data1$NPP2))

plot(log(data$NPP2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$NPP2)~I(data1$invTT-mean(data1$invTT)), pch = data1$Tank, col = data1$Tank, ylim = c(0,6))
plot(log(data1$NPP2)~I(data1$invT-(data1$invTT)), pch = data1$Tank, col = data1$Tank, ylim = c(0,6))
abline(3.00, -0.82, lwd = 2, col = 1)
plot(log(data1$NPP2)~data1$invT, pch = data1$Tank, col = data1$Tank)

library(tidyverse) group = Tank, color = trophic.level formula = log(data1$NPP2) ~ data1$invT, inherit.aes = FALSE

## here is a plot with basic regression lines fitted
ggplot(data = data1, aes(x = invT, y = log(NPP2))) + 
  theme_minimal() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black') +
  xlab("Temperature 1/kT") +
  ylab("ln(NPP)")

## I think the way to do this with lme results is to create a dataframe with those model output coefficients... or add the predictions of the model to the original dataset and plot those here...

plot(modNPP1)
predw <- predict(modNPP1, level = 0:1)
## ok, I think this predict is predicting each value in the dataset. by replotting this (hiding the points) but fitting the lines to these data, we'd get our modeled lines overlayed on the real data. and i can see here that the effect of time is missing. can i rewrite - in the model - the different temperature terms in terms if invT, so that i can plot it all against invT more straightforwardly?
## here is a plot with basic regression lines fitted
dim(data1)
dim(predw)
data1$Tank
data1$ID <- seq.int(nrow(data1))
predw$ID <- seq.int(nrow(predw))
data.pred <- left_join(data1, predw, by = "ID")

ggplot(data = data.pred, aes(x = invT, y = log(NPP2))) + 
  theme_minimal() +
  geom_point(aes(group = Tank.x, color = trophic.level)) +
  geom_smooth(data = data.pred, aes(x = invT, y = predict.fixed, group = Tank.x, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = data.pred, aes(x = invT, y = predict.fixed), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black') +
  xlab("Temperature 1/kT") +
  ylab("ln(NPP)")



### analysis with autocorrelation term
modNPP4a<-lm(log(NPP2)~1+I(invT-mean(invT))*trophic.level, data=data1, na.action=na.omit, correlation = corCAR1(0.2, ~ week|Tank)) 

display.brewer.all()

# for model fitting: 
modNPP1 <- lme(log(NPP2) ~ 1+ I(invT - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="REML")
summary(modNPP1)

confint(modNPP1)


################################################
## Does mass-specific NPP vary with temperature?  
################################################

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


## analysis
modNPPm0 <- lme(log(NPP.mass) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modNPPm1 <- lme(log(NPP.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modNPPm2 <- lme(log(NPP.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm4 <- lme(log(NPP.mass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm5 <- lme(log(NPP.mass) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)
model.sel(modNPPm0, modNPPm1, modNPPm2, modNPPm4, modNPPm5)

anova(modNPPm4, modNPPm2)
anova(modNPPm2, modNPPm1)
anova(modNPPm1, modNPPm2)

# for model fitting: 
modNPPm5r <- lme(log(NPP.mass) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 

summary(modNPPm5r)

######################################
## Does ER vary with temperature?  
#####################################
## figures 
hist(data$ER2)
data1 <- data[(data$ER2 >= 0),]
hist(log(data1$ER2))
#hist(log(data$calc.ER+1))

plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$ER2)~I(data1$invTT-mean(data1$invTT)), pch = 19, col = data1$trophic.level, ylim = c(0,6))
abline(4.27, -0.43, lwd = 2, col = 1)
abline((4.27+0.35), (-0.43), lwd = 2, col = 2)
abline((3.88-0.04), (-0.43), lwd = 2, col = 3)

#analysis
modER0<-lme(log(ER2) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modER1<-lme(log(ER2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modER2<-lme(log(ER2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER4<-lme(log(ER2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER5 <- lme(log(ER2) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

model.sel(modER0, modER1, modER2, modER4, modER5)

anova(modER4, modER2)
anova(modER1, modER2)
anova(modER1, modER0)

# for model fitting: 
modER2r<-lme(log(ER2) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
summary(modER2r)

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
modPP1<-lme(log(PP.biomass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modPP2<-lme(log(PP.biomass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modPP4<-lme(log(PP.biomass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modPP5 <- lme(log(PP.biomass) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

model.sel(modPP0, modPP1, modPP2, modPP4, modPP5) 

anova(modPP4, modPP5)

# for model fitting: 
modPP4r <- lme(log(PP.biomass) ~ 1 + I(invT - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
modPP5r <- lme(log(PP.biomass) ~ 1 + I(invT - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit)

m.avg <- model.avg(modPP4r, modPP5r)
summary(m.avg)


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


