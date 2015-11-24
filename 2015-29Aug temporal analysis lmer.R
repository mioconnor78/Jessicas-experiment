### MO made a new file on Aug29 from June file when I decided to take week out as a fixed effect, and instead model autocorrelation. 

### revised this on Nov 17 to redo models in lme4 so i can estimate confidence intervals on coefficients

### load libraries
#library(qpcR)
library(lme4)
library(MuMIn)

### set working directory and load data
data <- read.csv("./temporal_dataFEB12.csv")
dim(data)
head(data)
tail(data)
data <- data[-(241:255),]

## bring in data file with temperatures at each sampling time
o.data <- read.csv("./oxygen_temp_temporal.csv")
o.data <- o.data[,-(4:14)]
o.data <- o.data[,-2]
head(o.data)
dim(o.data)
data2 <- merge(data, o.data, by.x = c("week", "Tank"), by.y = c("week", "Tank"))
data <- data2

data3 <- merge(data, mass.data, by.x = c("week", "Tank"), by.y = c("week", "tank"), all=TRUE)
data3[is.na(data3)] <- 0
data <- data3

### load libraries, define variables and add columns
k <- 8.617342*10^-5  # eV/K
## create test temp
data$temp.t <- (data$average.temp + data$temp.dawn1 + data$temp.dusk)/3
data$invT <-  1/((data$temp.t + 273)*k)
data$invT <- as.numeric(as.character(data$invT))

data$PP.biomass <- (data$chla*55) #chla (ug/L)* 55 C in PP / 1 chla = ugPPC/L
data$ZP.carbon1 <- ifelse(data$trophic.level=='P',  0, data$M.uncorrz) # raw carbon
data$HeteroB <- ifelse(data$trophic.level=='PZN', data$M.uncorrz + 10.8/2, data$ZP.carbon1) # add notonectid weight
data$total.carbon <- data$PP.biomass + data$HeteroB
data$HA <- data$HeteroB/data$PP.biomass

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
data$ER.mass <- data$ER2/(data$M.corrT + data$PP.biomass) # ER on ummol 02/L/day/ugTPP

data <- data[data$week >= '4',]


### data prep complete

### ANALYSIS

## Does NPP vary with temperature?  
## figures on invT
data1 <- data[(data$NPP2 >= 0.3),] # three negative values and one very small value now, not sure what to do about them.
hist(data[(data$NPP2 >= 0.5),]$NPP2)
hist(log(data1$NPP2))
hist(log(data$NPP2))

plot(log(data$NPP2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$NPP2)~I(data1$invT-mean(data1$invT)), pch = 19, col = data1$trophic.level, ylim = c(0,6))
abline(3.00, -1.14, lwd = 2, col = 1)
abline((3.00+0.18), (-1.14-0.32), lwd = 2, col = 2)
abline((3.00+0.23), (-1.14+0.06), lwd = 2, col = 3)


## analysis
## determine need for random effects in the full model: 
modNPP4a <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPP4b<-lmer(log(NPP2)~1+I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit) 

anova(modNPP4a, modNPP4b)

modNPP0 <- lmer(log(NPP2) ~ 1 + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPP1 <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPP2 <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPP4 <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  

model.sel(modNPP0, modNPP1, modNPP2, modNPP4)
anova(modNPP4, modNPP1)
anova(modNPP1, modNPP2)
anova(modNPP2, modNPP0)

# for model fitting: 
modNPP1r <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = TRUE, na.action=na.omit)  
modNPP2r <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)  
modNPP4r <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)  

summary(modNPP4r)
confint(modNPP4r)

m.avg <- model.avg(modNPP2r, modNPP4r, modNPP1r)
summary(m.avg)

confint(m.avg)


# net ecosystem metabolism
hist(data$NEM)

data1 <- data[(data$NEM >= 0),]
plot(log(data1$NEM)~I(data1$invT-mean(data1$invT)), pch = 19, col = data1$trophic.level)
abline(-0.13, 0.04, lwd = 2, col = 1)
abline((-0.13+0.13), (0.04+0.19), lwd = 2, col = 2)
abline((-0.13-0.15), (0.04+0.29), lwd = 2, col = 3)
plot((data$NEM)~data$week, pch = 19, col = data$trophic.level)
plot((data$NEM)~data$Tank, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modNEMa <- lmer(log(NEM) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data1, REML = FALSE, na.action=na.omit)  
modNEMb <- lmer(log(NEM) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit) 

anova(modNEMa, modNEMb)

modNEM0 <- lmer(log(NEM) ~ 1 + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modNEM1 <- lmer(log(NEM) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modNEM2 <- lmer(log(NEM) ~ 1 + I(invT-mean(invT))+trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modNEM4 <- lmer(log(NEM) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)

model.sel(modNEM0, modNEM1, modNEM2, modNEM4)
anova(modNEM4, modNEM2)
anova(modNEM0, modNEM2)
anova(modNEM1, modNEM0)


# for model fitting: 
modNEM4r <- lmer(log(NEM) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)
modNEM2r <- lmer(log(NEM) ~ 1 + I(invT-mean(invT))+trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)

## the best model is modNEM4r. Use this one for making figures. 
m.avg <- model.avg(modNEM4r, modNEM2r)
summary(m.avg)
confint(m.avb)

summary(modNEM4r)
confint(modNEM4r)

## Does mass-specific NPP vary with temperature?  
## figures  
### NPP in umol/L/day per ugC
data1 <- data[(data$NPP2 >= 0.3),]
hist(data1$NPP.mass)
hist(log(data1$NPP.mass))
# hist(log(data$NPP.mass+.001)) # no reason to add a number, there are no -inf values.

plot(log(data$NPP.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$NPP.mass)~I(data1$invT-mean(data1$invT)), pch = 19, col = data1$trophic.level)
abline(-2.83, -1.89, lwd = 2, col = 1)
abline((-2.83+1.36), (-1.89-1.41), lwd = 2, col = 2)
abline((-2.83+0.46), (-1.89-0.05), lwd = 2, col = 3)


## analysis
## determine need for random effects in the full model: 
modNPPma <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPPmb <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit) 

anova(modNPPma, modNPPmb)

modNPPm0 <- lmer(log(NPP.mass) ~ 1 + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPPm1 <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPPm2 <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modNPPm4 <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  

model.sel(modNPPm0, modNPPm1, modNPPm2, modNPPm4)

anova(modNPPm4, modNPPm2)
anova(modNPPm2, modNPPm1)
anova(modNPPm1, modNPPm2)

# for model fitting: 
modNPPm4r <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)   
  
summary(modNPPm4r)
confint(modNPPm4r)


## Does ER vary with temperature?  
## figures 
hist(data$ER2)
data1 <- data[(data$ER2 >= 0),]
hist(log(data1$ER2))
#hist(log(data$calc.ER+1))

plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$ER2)~I(data1$invT-mean(data1$invT)), pch = 19, col = data1$trophic.level, ylim = c(0,6))
abline(3.88, -0.75, lwd = 2, col = 1)
abline((3.88+0.49), (-0.75+0.29), lwd = 2, col = 2)
abline((3.88+0.01), (-0.75+0.62), lwd = 2, col = 3)

#analysis
## determine need for random effects in the full model: 
modERa <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data1, REML = FALSE, na.action=na.omit)  
modERb <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit) 

anova(modERa, modERb)

modER0 <- lmer(log(ER2) ~ 1 + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modER1 <- lmer(log(ER2) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modER2 <- lmer(log(ER2) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modER4 <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)

model.sel(modER0, modER1, modER2, modER4)

anova(modER4, modER2)
anova(modER1, modER2)
anova(modER1, modER0)

# for model fitting: 
modER4r <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)
modER2r <- lmer(log(ER2) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)
m.avg <- model.avg(modER4r, modER2r)
summary(m.avg)
confint(m.avg)

summary(modER4r)
confint(modER4r)

## Does mass specific ER vary with temperature?  
## figures 
data1 <- data[(data$ER2 >= 0),]
hist(data1$ER.mass)
hist(log(data1$ER.mass))

plot(log(data$ER.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass+1)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modERma <- lmer(log(ER.mass+1) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data1, REML = FALSE, na.action=na.omit)  #, control = lmerControl(optimizer = "Nelder_Mead")
modERmb <- lmer(log(ER.mass+1) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  #, control = lmerControl(optimizer = "Nelder_Mead")

anova(modERma, modERmb)

modERm0 <- lmer(log(ER.mass+1) ~ 1 + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modERm1 <- lmer(log(ER.mass+1) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modERm2 <- lmer(log(ER.mass+1) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  
modERm4 <- lmer(log(ER.mass+1) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)  

model.sel(modERm0, modERm1, modERm2, modERm4)

# for model fitting: 
modERm4r <- lmer(log(ER.mass+1) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit) 
summary(modERm4r)
confint(modERm4r)




##### BIOMASS RESULTS #####

#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))

plot(log(data$PP.biomass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modPPa <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = FALSE, na.action=na.omit)  
modPPb <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = FALSE, na.action=na.omit) 
anova(modPPa, modPPb)

# modPPc<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit)
# anova(modPPc, modPPb)

modPP0 <- lmer(log(PP.biomass) ~ 1 + (I(invT-mean(invT))|week), data=data, REML = FALSE, na.action=na.omit)
modPP1 <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT)) + (I(invT-mean(invT))|week), data=data, REML = FALSE, na.action=na.omit)
modPP2 <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))+trophic.level + (I(invT-mean(invT))|week), data=data, REML = FALSE, na.action=na.omit)
modPP4 <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = FALSE, na.action=na.omit)

model.sel(modPP0, modPP1, modPP2, modPP4) 

anova(modPP4, modPP2)

# for model fitting: 
modPP4r <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

summary(modPP4r)
confint(modPP4r)

## Does zooplankton carbon vary with temperature?  

dataz <- data[(which(data$trophic.level != 'P')),] # because we don't have observations for the P treatments (we're just assuming they are 0) I don't htink we should analyze those tanks here.
hist(dataz$zoo.ug.carbon.liter)
hist(log(dataz$zoo.ug.carbon.liter))
hist(dataz$M.uncorrz)
hist(log(dataz$M.uncorrz+1))

plot(log(dataz$M.uncorrz+1)~dataz$invT, pch = 19, col = dataz$trophic.level)

## analysis
## determine need for random effects in the full model: 
modZa <- lmer(log(M.uncorrz0 + 1) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=dataz, REML = FALSE, na.action=na.omit)  
modZb <- lmer(log(M.uncorrz0 + 1) ~ 1 + I(invT-mean(invT))*trophic.level + 1|week, data=dataz, REML = FALSE, na.action=na.omit)  
anova(modZa, modZb)





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
### total carbon
hist(log(data$total.carbon))

modTCb<-lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = FALSE, na.action=na.omit)
modTCc<-lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = FALSE, na.action=na.omit)
anova(modTCc, modTCb)

modTC0 <- lmer(log(total.carbon) ~ 1 + (1|week), data=data, REML = FALSE, na.action=na.omit)
modTC1 <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT)) + (1|week), data=data, REML = FALSE, na.action=na.omit)
modTC2 <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data, REML = FALSE, na.action=na.omit)
modTC4 <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = FALSE, na.action=na.omit)

model.sel(modTC0, modTC1, modTC2, modTC4)
anova(modTC2, modTC4)
anova(modTC2, modTC1)

# for model fitting: 
modTC4r <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = TRUE, na.action=na.omit)

summary(modTC4r)
confint(modTC4r)

#### HA ratio
data1 <- data[which(data$HA > '0'),]
hist(log(data1$HA))

plot(log(data1$HA)~((data1$invT)), pch = 19, col = data1$trophic.level)

modHAb <- lmer(log(HA) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data1, REML = FALSE, na.action=na.omit)
modHAc <- lmer(log(HA) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)
anova(modHAb, modHAc)

modHA0 <- lmer(log(HA) ~ 1 + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modHA1 <- lmer(log(HA) ~ 1 + I(invT-mean(invT)) + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modHA2 <- lmer(log(HA) ~ 1 + I(invT-mean(invT)) + trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)
modHA4 <- lmer(log(HA) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = FALSE, na.action=na.omit)

model.sel(modHA0, modHA1, modHA2, modHA4)
anova(modHA2, modHA4)
anova(modHA2, modHA1)

# for model fitting: 
modHA4r <- lmer(log(HA) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)
modHA2r <- lmer(log(HA) ~ 1 + I(invT-mean(invT)) + trophic.level + (I(invT-mean(invT))|week), data=data1, REML = TRUE, na.action=na.omit)

m.avg <- model.avg(modHA4r, modHA2r)
summary(m.avg)
confint(m.avg)

summary(modHA4r)
confint(modHA4r)
