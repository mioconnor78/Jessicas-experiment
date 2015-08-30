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
plot(log(data$NPP2)~data$invT, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modNPP4a<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modNPP4b<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 

anova(modNPP4a, modNPP4b)

modNPP4c<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modNPP4c, modNPP4b)

modNPP0<-lme(log(NPP2) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modNPP1<-lme(log(NPP2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modNPP2<-lme(log(NPP2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modNPP4<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

model.sel(modNPP0, modNPP1, modNPP2, modNPP4)

# for model fitting: 
modNPP2r <- lme(log(NPP2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modNPP1r <- lme(log(NPP2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modNPP4r <- lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

m.avg <- model.avg(modNPP2r, modNPP1r, modNPP4r)
summary(m.avg)

# net ecosystem metabolism
hist(data$NEM)

data1 <- data #[(data$NEM<=14.2),]
plot((data$NEM)~data$invT, pch = 19, col = data$trophic.level)
abline(h=(1))
plot((data$NEM)~data$week, pch = 19, col = data$trophic.level)
plot((data$NEM)~data$Tank, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modNEMa<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modNEMb<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 

anova(modNEMa, modNEMb)

modNEMc<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modNEMc, modNEMb)

modNEM0<-lme(log(NEM) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNEM1<-lme(log(NEM)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNEM2<-lme(log(NEM)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNEM4<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modNEM0, modNEM1, modNEM2, modNEM4)

# for model fitting: 
modNEM2r<-lme(log(NEM)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
modNEM4r<-lme(log(NEM)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)

m.avg <- model.avg(modNEM2r, modNEM4r)
summary(m.avg)


## Does mass-specific NPP vary with temperature?  
## figures  
### NPP in umol/L/day per ugC
hist(data$NPP.mass)
hist(log(data$NPP.mass))
# hist(log(data$NPP.mass+.001)) # no reason to add a number, there are no -inf values.

plot(log(data$NPP.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~data$invT, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modNPPma<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modNPPmb<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 

anova(modNPPma, modNPPmb)

modNPPmc<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modNPPmc, modNPPmb)

modNPPm0<-lme(log(NPP.mass) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm1<-lme(log(NPP.mass)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm2<-lme(log(NPP.mass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm4<-lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modNPPm0, modNPPm1, modNPPm2, modNPPm4)

# for model fitting: 
modNPPm4r <- lme(log(NPP.mass)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)

summary(modNPPm4r)

## Does ER vary with temperature?  
## figures 
hist(data$ER2)
hist(log(data$ER2))
#hist(log(data$calc.ER+1))

plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$invT, pch = 19, col = data$trophic.level)

#analysis
## determine need for random effects in the full model: 
modERa<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modERb<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 

anova(modERa, modERb)

modERc<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modERc, modERb)

modER0<-lme(log(ER2) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER1<-lme(log(ER2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER2<-lme(log(ER2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER4<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modER0, modER1, modER2, modER4)

# for model fitting: 
modER4r<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
summary(modER4r)


## Does mass specific ER vary with temperature?  
## figures 
hist(data$ER.mass)
hist(log(data$ER.mass))

plot(log(data$ER.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$invT, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modERma<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modERmb<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 

anova(modERma, modERmb)

modERmc<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modERmc, modERmb)

modERm0<-lme(log(ER.mass) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modERm1<-lme(log(ER.mass)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modERm2<-lme(log(ER.mass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modERm4<-lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

model.sel(modERm0, modERm1, modERm2, modERm4)

# for model fitting: 
modERm2r <- lme(log(ER.mass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modERm4r <- lme(log(ER.mass)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
m.avg <- model.avg(modERm2r, modERm4r)
summary(m.avg)








##### BIOMASS RESULTS #####

#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))

plot(log(data$PP.biomass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~data$invT, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modPPa<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modPPb<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 
anova(modPPma, modPPmb)

modPPc<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modPPc, modPPb)

modPP0<-lme(log(PP.biomass) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modPP1<-lme(log(PP.biomass)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modPP2<-lme(log(PP.biomass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modPP4<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

model.sel(modPP0, modPP1, modPP2, modPP4)

# for model fitting: 
modPP1r<-lme(log(PP.biomass)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modPP2r<-lme(log(PP.biomass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modPP4r<-lme(log(PP.biomass)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

m.avg <- model.avg(modPP1r, modPP2r, modPP4r)
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
modZa<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modZb<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 
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

plot(log(data[data$week >= '4',]$total.carbon)~data[data$week >= '4',]$Tank, pch = 19, col = data$trophic.level)
plot(log(data[data$week >= '4',]$total.carbon)~data[data$week >= '4',]$week, pch = 19, col = data$trophic.level)
plot(log(data[data$week >= '4',]$total.carbon)~data[data$week >= '4',]$invT, pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modTCa<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))  
modTCb<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank)) 
anova(modTCa, modTCb)

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


