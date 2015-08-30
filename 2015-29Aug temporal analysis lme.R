### MO made a new file on Aug29 from June file when I decided to take week out as a fixed effect, and instead model autocorrelation. 

### load libraries
#library(qpcR)
library(nlme)
library(MuMIn)

### set working directory and load data
data <- read.csv("./temporal_dataFEB12.csv")
dim(data)
head(data)
tail(data)
data <- data[-(241:255),]

### load libraries, define variables and add columns
k <- 8.617342*10^-5  # eV/K
data$invT <-  1/((data$average.temp + 273)*k)
data$invT <- as.numeric(as.character(data$invT))
#data$kT <-  ((data$average.temp + 273)*k)
#data$kT <- as.numeric(as.character(data$kT))*100

data$PP.biomass <- (data$chla*55) #chla (ug/L)* 55 C in PP / 1 chla = ugPPC/L
data$ZP.carbon1 <- ifelse(data$trophic.level=='P',  0, data$zoo.ug.carbon.liter) # for adding
data$total.carbon <- data$PP.biomass + data$ZP.carbon1 #I'm assuming 0 for ZP in P treatments here.

### estimating NPP and ER from the raw data: 
#calculate NPP and ER (hourly), in terms of umol O2 / l / hr, following yvon durochers 2010.
# O2 has molar mass of 32g/mol. so 1 umol = 32 ug. so take ug/32
data$NPP2 <- ((data$dusk - data$dawn1)*1000)/(32)  # oxygen produced umol / L /day, net all respiration. raw data is mg/L. 
data$ER2 <- -(24/data$hours2)*((data$dawn2 - data$dusk)*1000)/(32)  # amount of oxygen consumed per day via respiaration. negative to get the change in oxygen umol / L /day; oxygen used in the dark and daylight. MeanER can be greater than meanNPP, because NPP reflects ER already.
data$GPP <- data$NPP2+(data$ER2/24)*data$hours1 # daily oxygen production (NPP2) + estimated daytime community respiration (daily R / 24 * hours daylight)
data$NEM <- data$ER2/data$GPP  # following Yvon Durochers 2010. NEM > 1 means the system is respiring more than it's fixing per day. This does not need to be logged.
data$NPP.mass <- data$NPP2 / (data$PP.biomass)  # NPP on ummol 02/L/day/ugCPP
data$ER.mass <- data$ER2/(data$total.carbon) # ER on ummol 02/L/day/ugTPP

data <- data[data$week >= '4',]

## temporal

## Does NPP vary with temperature?  
## figures on invT
hist(data$NPP2)
hist(log(data$NPP2))

plot(log(data$NPP2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$NPP2)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)
abline(3.00, -1.14, lwd = 2, col = 1)
abline((3.00+0.18), (-1.14-0.32), lwd = 2, col = 2)
abline((3.00+0.23), (-1.14+0.06), lwd = 2, col = 3)


## analysis
## determine need for random effects in the full model: 
modNPP4a<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPP4b<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 

anova(modNPP4a, modNPP4b)

modNPP0<-lme(log(NPP2) ~ 1, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPP1<-lme(log(NPP2)~1+I(invT-mean(invT)), random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPP2<-lme(log(NPP2)~1+I(invT-mean(invT))+trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPP4<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  

model.sel(modNPP0, modNPP1, modNPP2, modNPP4)
anova(modNPP2, modNPP4)
anova(modNPP4, modNPP1)
anova(modNPP1, modNPP0)

# for model fitting: 
modNPP2r <- lme(log(NPP2)~1+I(invT-mean(invT))+trophic.level, random=~I(invT-mean(invT))|week, data=data, method="REML", na.action=na.omit)
modNPP4r <- lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random=~I(invT-mean(invT))|week, data=data, method="REML", na.action=na.omit)

m.avg <- model.avg(modNPP2r, modNPP4r)
summary(m.avg)

# net ecosystem metabolism
hist(data$NEM)

data1 <- data #[(data$NEM<=14.2),]
plot(log(data$NEM)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)
abline(-0.13, 0.04, lwd = 2, col = 1)
abline((-0.13+0.13), (0.04+0.19), lwd = 2, col = 2)
abline((-0.13-0.15), (0.04+0.29), lwd = 2, col = 3)
plot((data$NEM)~data$week, pch = 19, col = data$trophic.level)
plot((data$NEM)~data$Tank, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modNEMa<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNEMb<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 

anova(modNEMa, modNEMb)

modNEM0<-lme(log(NEM) ~ 1, random=~1|week, data=data, method="ML", na.action=na.omit)
modNEM1<-lme(log(NEM)~1+I(invT-mean(invT)), random=~1|week, data=data, method="ML", na.action=na.omit)
modNEM2<-lme(log(NEM)~1+I(invT-mean(invT))+trophic.level, random=~1|week, data=data, method="ML", na.action=na.omit)
modNEM4<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random=~1|week, data=data, method="ML", na.action=na.omit)

model.sel(modNEM0, modNEM1, modNEM2, modNEM4)
anova(modNEM4, modNEM2)
anova(modNEM1, modNEM2)
anova(modNEM1, modNEM0)


# for model fitting: 
modNEM4r<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random=~1|week, data=data, method="REML", na.action=na.omit)

summary(modNEM4r)


## Does mass-specific NPP vary with temperature?  
## figures  
### NPP in umol/L/day per ugC
hist(data$NPP.mass)
hist(log(data$NPP.mass))
# hist(log(data$NPP.mass+.001)) # no reason to add a number, there are no -inf values.

plot(log(data$NPP.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)
abline(-2.83, -1.89, lwd = 2, col = 1)
abline((-2.83+1.36), (-1.89-1.41), lwd = 2, col = 2)
abline((-2.83+0.46), (-1.89-0.05), lwd = 2, col = 3)


## analysis
## determine need for random effects in the full model: 
modNPPma<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPPmb<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 

anova(modNPPma, modNPPmb)

modNPPm0<-lme(log(NPP.mass) ~ 1, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPPm1<-lme(log(NPP.mass)~1+I(invT-mean(invT)), random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPPm2<-lme(log(NPP.mass)~1+I(invT-mean(invT))+trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modNPPm4<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  

model.sel(modNPPm0, modNPPm1, modNPPm2, modNPPm4)

# for model fitting: 
modNPPm4r <- lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="REML", na.action=na.omit)  

summary(modNPPm4r)

## Does ER vary with temperature?  
## figures 
hist(data$ER2)
hist(log(data$ER2))
#hist(log(data$calc.ER+1))

plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)
abline(3.88, -0.75, lwd = 2, col = 1)
abline((3.88+0.49), (-0.75+0.29), lwd = 2, col = 2)
abline((3.88+0.01), (-0.75+0.62), lwd = 2, col = 3)

#analysis
## determine need for random effects in the full model: 
modERa<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modERb<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 

anova(modERa, modERb)

modER0<-lme(log(ER2) ~ 1, random=~1|week, data=data, method="ML", na.action=na.omit)
modER1<-lme(log(ER2)~1+I(invT-mean(invT)), random=~1|week, data=data, method="ML", na.action=na.omit)
modER2<-lme(log(ER2)~1+I(invT-mean(invT))+trophic.level, random=~1|week, data=data, method="ML", na.action=na.omit)
modER4<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random=~1|week, data=data, method="ML", na.action=na.omit)

model.sel(modER0, modER1, modER2, modER4)

anova(modER4, modER2)
anova(modER1, modER2)
anova(modER1, modER0)

# for model fitting: 
modER4r<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random=~1|week, data=data, method="REML", na.action=na.omit)
summary(modER4r)


## Does mass specific ER vary with temperature?  
## figures 
hist(data$ER.mass)
hist(log(data$ER.mass))

plot(log(data$ER.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)
abline(-1.89, -1.49, lwd = 2, col = 1)
abline((-1.89+1.54), (-1.49-0.74), lwd = 2, col = 2)
abline((-1.89+0.21), (-1.49+0.43), lwd = 2, col = 3)

## analysis
## determine need for random effects in the full model: 
modERma<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modERmb<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 

anova(modERma, modERmb)

modERm0<-lme(log(ER.mass) ~ 1, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modERm1<-lme(log(ER.mass)~1+I(invT-mean(invT)), random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modERm2<-lme(log(ER.mass)~1+I(invT-mean(invT))+trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modERm4<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  

model.sel(modERm0, modERm1, modERm2, modERm4)

# for model fitting: 
modERm4r <- lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="REML", na.action=na.omit)
summary(modERm4r)








##### BIOMASS RESULTS #####

#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))

plot(log(data$PP.biomass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)

abline(5.87, 0.79, col = 1, lwd = 2)
abline((5.87-1.19), (0.79+1.05), col = 2, lwd = 2)
abline((5.87-0.23), (0.79+0.08), col = 3, lwd = 2)

## analysis
## determine need for random effects in the full model: 
modPPa<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modPPb<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 
anova(modPPa, modPPb)

# modPPc<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit)
# anova(modPPc, modPPb)

modPP0<-lme(log(PP.biomass) ~ 1, random=~1|week, data=data, method="ML", na.action=na.omit)
modPP1<-lme(log(PP.biomass)~1+I(invT-mean(invT)), random=~1|week, data=data, method="ML", na.action=na.omit)
modPP2<-lme(log(PP.biomass)~1+I(invT-mean(invT))+trophic.level, random=~1|week, data=data, method="ML", na.action=na.omit)
modPP4<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random=~1|week, data=data, method="ML", na.action=na.omit)

model.sel(modPP0, modPP1, modPP2, modPP4) 

anova(modPP4, modPP2)

# for model fitting: 
modPP4r<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random=~1|week, data=data, method="REML", na.action=na.omit)

summary(modPP4r)


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


