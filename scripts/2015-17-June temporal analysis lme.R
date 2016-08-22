
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
modNPP0<-lme(log(NPP2) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP1<-lme(log(NPP2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP2<-lme(log(NPP2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP3<-lme(log(NPP2)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP4<-lme(log(NPP2)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP5<-lme(log(NPP2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP6<-lme(log(NPP2)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modNPP0, modNPP1, modNPP2, modNPP3, modNPP4, modNPP5, modNPP6)

# for model fitting: 
modNPP6r<-lme(log(NPP2)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)

summary(modNPP6r)


# net ecosystem metabolism
hist(data$NEM)
#hist(sqrt(data$NEM))
#hist(log(data$NEM))
#hist(log(data$NEM+1))

data1 <- data #[(data$NEM<=14.2),]
plot((data$NEM)~data$invT, pch = 19, col = data$trophic.level)
abline(h=(1))
plot((data$NEM)~data$week, pch = 19, col = data$trophic.level)
plot((data$NEM)~data$Tank, pch = 19, col = data$trophic.level)

## analysis
modNEM0<-lme(NEM ~ 1, random=~1|Tank, data=data1, method="ML", na.action=na.omit)
modNEM1<-lme(NEM~1+I(invT-mean(invT)), random=~1|Tank, data=data1, method="ML", na.action=na.omit)
modNEM2<-lme(NEM~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data1, method="ML", na.action=na.omit)
modNEM3<-lme(NEM~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data1, method="ML", na.action=na.omit)
modNEM4<-lme(NEM~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data1, method="ML", na.action=na.omit)
modNEM5<-lme(NEM~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data1, method="ML", na.action=na.omit)
modNEM6<-lme(NEM~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data1, method="ML", na.action=na.omit)

model.sel(modNEM0, modNEM1, modNEM2, modNEM3, modNEM4, modNEM5, modNEM6)

anova(modNEM6, modNEM4)
anova(modNEM4, modNEM3)
anova(modNEM3, modNEM2)
anova(modNEM2, modNEM5)
anova(modNEM5, modNEM1)
anova(modNEM1, modNEM0)

# for model fitting: 
modNEM2r<-lme(log(NEM)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data1, method="REML", na.action=na.omit)
modNEM3r<-lme(log(NEM)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data1, method="REML", na.action=na.omit)
modNEM6r<-lme(NEM~I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data1, method="ML", na.action=na.omit)

m.avg <- model.avg(modNEM2r, modNEM3r)
summary(modNEM6r)


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
modNPPm0<-lme(log(NPP.mass) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm1<-lme(log(NPP.mass)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm2<-lme(log(NPP.mass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm3<-lme(log(NPP.mass)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm4<-lme(log(NPP.mass)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm5<-lme(log(NPP.mass)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm6<-lme(log(NPP.mass)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modNPPm0, modNPPm1, modNPPm2, modNPPm3, modNPPm4, modNPPm5, modNPPm6)

# for model fitting:
modNPPm5r<-lme(log(NPP.mass)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
summary(modNPPm5r)


## Does ER vary with temperature?  
## figures 
hist(data$ER2)
hist(log(data$ER2))
#hist(log(data$calc.ER+1))

plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$invT, pch = 19, col = data$trophic.level)

## analysis
modER0<-lme(log(ER2) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER1<-lme(log(ER2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER2<-lme(log(ER2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER3<-lme(log(ER2)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER4<-lme(log(ER2)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER5<-lme(log(ER2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER6<-lme(log(ER2)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modER0, modER1, modER2, modER3, modER4, modER5, modER6)

anova(modER5, modER4)
anova(modER4, modER2)
anova(modER2, modER3)
anova(modER3, modER6)
anova(modER6, modER1)
anova(modER1, modER0)

## for fitting
modER5r<-lme(log(ER2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
summary(modER5r)



## Does mass specific ER vary with temperature?  
## figures 
hist(data$ER.mass)
hist(log(data$ER.mass))

plot(log(data$ER.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$invT, pch = 19, col = data$trophic.level)

## analysis
modERm0<-lme(log(ER.mass) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm1<-lme(log(ER.mass)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm2<-lme(log(ER.mass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm3<-lme(log(ER.mass)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm4<-lme(log(ER.mass)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm5<-lme(log(ER.mass)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm6<-lme(log(ER.mass)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modERm0, modERm1, modERm2, modERm3, modERm4, modERm5, modERm6)

anova(modERm5, modERm4)
anova(modERm5, modERm6)
anova(modERm6, modERm3)
anova(modERm3, modERm2)
anova(modERm2, modERm1)
anova(modERm1, modERm0)

## for fitting
modERm4<-lme(log(ER.mass)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
modERm5r<-lme(log(ER.mass)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)

summary(modERm5r)








##### BIOMASS RESULTS #####

## Does chla vary with temperature?  
## MO skipped this for now b/c it's redundant with PP biomass (below) and we probably shouldn't analyze both.
### MO: went through this and made minor changes on July 28, 2015
## figures 
# data1 <- data[(which(data$week!='2')),] #could redo without week2, to see if that's driving the need for week effects.
hist(data$chla)
hist(log(data$chla))
#hist(log(data$chla+1))
plot(log(data$chla)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$chla)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$chla)~data$invT, pch = 19, col = data$trophic.level)

## analysis
modchl0<-lme(log(chla)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl1<-lme(log(chla)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl2<-lme(log(chla)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl3<-lme(log(chla)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl4<-lme(log(chla)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl5<-lme(log(chla)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl6<-lme(log(chla)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modchl0, modchl1, modchl2, modchl3, modchl4, modchl5, modchl6)

anova(modchl6, modchl5)
anova(modchl5, modchl4)
anova(modchl4, modchl3)
anova(modchl3, modchl2)
anova(modchl2, modchl1)
anova(modchl1,modchl0)

## refit best model with reml
modchl6b <- lme(log(chla)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
summary(modchl6b)




#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
datap <- data[data$week >= '4',]
hist(datap$PP.biomass)
hist(log(datap$PP.biomass))

plot(log(datap$PP.biomass)~datap$Tank, pch = 19, col = datap$trophic.level)
plot(log(datap$PP.biomass)~datap$week, pch = 19, col = datap$trophic.level)
plot(log(datap$PP.biomass)~datap$invT, pch = 19, col = datap$trophic.level)

## analysis
modBp0<-lme(log(PP.biomass)~1, random=~1|Tank, data=datap, method="ML", na.action=na.omit)
modBp1<-lme(log(PP.biomass)~1+I(invT-mean(invT)), random=~1|Tank, data=datap, method="ML", na.action=na.omit)
modBp2<-lme(log(PP.biomass)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=datap, method="ML", na.action=na.omit)
modBp3<-lme(log(PP.biomass)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=datap, method="ML", na.action=na.omit)
modBp4<-lme(log(PP.biomass)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=datap, method="ML", na.action=na.omit)
modBp5<-lme(log(PP.biomass)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=datap, method="ML", na.action=na.omit)
modBp6<-lme(log(PP.biomass)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=datap, method="ML", na.action=na.omit)

model.sel(modBp0, modBp1, modBp2, modBp3, modBp4, modBp5, modBp6)

anova(modBp4, modBp3)
anova(modBp3, modBp5)
anova(modBp5, modBp6)
anova(modBp6, modBp2)
anova(modBp2, modBp1)
anova(modBp1, modBp0)

m.avg <- model.avg(modBp4, modBp3, modBp5)
summary(m.avg)


## Does zooplankton carbon vary with temperature?  
## figures 
dataz <- data[(which(data$trophic.level != 'P')),] # because we don't have observations for the P treatments (we're just assuming they are 0) I don't htink we should analyze those tanks here.
hist(data$zoo.ug.carbon.liter)
hist(log(data$zoo.ug.carbon.liter))
hist(log(data$zoo.ug.carbon.liter+1))
plot(log(data$zoo.ug.carbon.liter)~data$Tank, pch = 19, col = data$trophic.level)

plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$Tank, pch = 19, col = dataz$trophic.level)
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$week, pch = 19, col = dataz$trophic.level)
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$invT, pch = 19, col = dataz$trophic.level)

## analysis
modzpc0<-lme(log(zoo.ug.carbon.liter+1)~1, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
modzpc1<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT)), random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
modzpc2<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
modzpc3<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
modzpc4<-lme(log(zoo.ug.carbon.liter+1)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
modzpc5<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
modzpc6<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)

model.sel(modzpc0, modzpc1, modzpc2, modzpc3, modzpc4, modzpc5, modzpc6)

anova(modzpc4, modzpc5)
anova(modzpc4, modzpc1)
anova(modzpc1, modzpc6)
anova(modzpc1, modzpc2)
anova(modzpc1, modzpc3)
anova(modzpc1, modzpc0)

modzpc4r<-lme(log(zoo.ug.carbon.liter+1)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=dataz[dataz$week >= '4',], method="ML", na.action=na.omit)
summary(modzpc4r)


## Does zooplankton density vary with temperature?  ## should probably be analyzed as count data.
## figures 
hist(dataz$total.zoo.abundance.liter, breaks = 40)
hist(log(dataz$total.zoo.abundance.liter))
hist(log(dataz$total.zoo.abundance.liter+1), breaks = 40)

plot((dataz$total.zoo.abundance.liter)~dataz$Tank, pch = 19, col = dataz$trophic.level)
plot((dataz$total.zoo.abundance.liter)~dataz$week, pch = 19, col = dataz$trophic.level)
plot((dataz$total.zoo.abundance.liter)~dataz$invT, pch = 19, col = dataz$trophic.level)

## analysis
## this data is strongly right skewed, so we can't really do an linear mixed effects model any more. we want poisson distributed errors. From my quick googling, I think this is an option, but you might check a stats book. Or maybe lmer?
modzpd0<-glmmPQL((total.zoo.abundance.liter)~1, random=~1|Tank, family = poisson, data=dataz, na.action=na.omit)
modzpd1<-glmmPQL((total.zoo.abundance.liter)~1+invT, random=~1|Tank, data=dataz, family = poisson, na.action=na.omit)
modzpd2<-glmmPQL((total.zoo.abundance.liter)~1+invT+trophic.level, random=~1|Tank, data=dataz, family = poisson, na.action=na.omit)
modzpd3<-lme(log(total.zoo.abundance.liter+1)~1+invT+week+trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpd4<-lme(log(total.zoo.abundance.liter+1)~1+week+invT*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpd5<-lme(log(total.zoo.abundance.liter+1)~1+invT*week+invT*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpd6<-lme(log(total.zoo.abundance.liter+1)~1+invT*week*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)

model.sel(modzpd0, modzpd1, modzpd3, modzpd4, modzpd5, modzpd6)

anova(modzp0, modzp1)
anova(modzp1, modzp2)
anova(modzp2, modzp3)
anova(modzp3, modzp4)

AIC(modzp0, modzp1, modzp2, modzp3, modzp4)
AICs <- as.data.frame(cbind(AICc(modzp0),AICc(modzp1), AICc(modzp2), AICc(modzp3), AICc(modzp4)))
akaike.weights(AICs)

logLik(modzp0)
logLik(modzp1)
logLik(modzp2)
logLik(modzp3)
logLik(modzp4)

model2zp0 <- update(modzp0, correlation = corAR1())
model2zp1 <- update(modzp1, correlation = corAR1())
model2zp2 <- update(modzp2, correlation = corAR1())
model2zp3 <- update(modzp3, correlation = corAR1())
model2zp4 <- update(modzp4, correlation = corAR1())

model2zp0
model2zp1
model2zp2
model2zp3
model2zp4

summary(modzp4)
modzp3<-lm(log(total.zoo.abundance.liter+1)~1+invT+week+trophic.level, data=data, na.action=na.omit)

## add lines to plot
abline(-5.07736, 0.137716, lty = 2, lwd = 3, col = 'brown')
abline(-16.995649, 0.437164, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 0.1, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')



## Does daphnia:copepod  vary with temperature?  
## figures 
hist(dataz$Daphnia.Copepod.Ratio)
hist(log(dataz$Daphnia.Copepod.Ratio + 0.1))

plot(log(dataz$Daphnia.Copepod.Ratio + 0.1)~dataz$Tank, pch = 19, col = dataz$trophic.level)
plot(log(dataz$Daphnia.Copepod.Ratio + 0.1)~dataz$week, pch = 19, col = dataz$trophic.level)
plot(log(dataz$Daphnia.Copepod.Ratio + 0.1)~dataz$invT, pch = 19, col = dataz$trophic.level)

plot(log(data$Daphnia.Copepod.Ratio)~data$invT)
plot(log(data$Daphnia.Copepod.Ratio+1)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-0.5,3.5), xlim=c(38.0,41),  xlab='inv Temperature (C)', ylab='ln (Daphnia: Copepod)') 
axis(1, at=c(38.0, 38.5,39, 39.5, 40,40.5, 41), pos=-0.5, lwd=2, cex.lab=1.5)
axis(2, at=c(-0.5,  0,1.0, 2.0, 3.0, 3.5), pos=38.0, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='PZ'),]$Daphnia.Copepod.Ratio+1)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$Daphnia.Copepod.Ratio+1)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modDCR0<-lme(log(total.zoo.abundance.liter+1)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modDCR1<-lme(log(total.zoo.abundance.liter+1)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modDCR2<-lme(log(total.zoo.abundance.liter+1)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modDCR3<-lme(log(total.zoo.abundance.liter+1)~1+invT+ week + trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modDCR4<-lme(log(total.zoo.abundance.liter+1)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modDCR0, modDCR1)
anova(modDCR1, modDCR2)
anova(modDCR2, modDCR3)
anova(modDCR3, modDCR4)

AIC(modDCR0, modDCR1, modDCR2, modDCR3, modDCR4)
AICs <- as.data.frame(cbind(AICc(modDCR0),AICc(modDCR1), AICc(modDCR2), AICc(modDCR3), AICc(modDCR4)))
akaike.weights(AICs)

logLik(modDCR0)
logLik(modDCR1)
logLik(modDCR2)
logLik(modDCR3)
logLik(modDCR4)

model2DCR0 <- update(modDCR0, correlation = corAR1())
model2DCR1 <- update(modDCR1, correlation = corAR1())
model2DCR2 <- update(modDCR2, correlation = corAR1())
model2DCR3 <- update(modDCR3, correlation = corAR1())
model2DCR4 <- update(modDCR4, correlation = corAR1())

model2DCR0
model2DCR1
model2DCR2
model2DCR3
model2DCR4

summary(modDCR4)
modDCR3<-lm(log(total.zoo.abundance.liter+1)~1+invT+ week + trophic.level, data=data, na.action=na.omit)
summary(modDCR3)

## add lines to plot
abline(-5.07706, 0.137722, lty = 2, lwd = 3, col = 'brown')
abline(-16.995649, 0.4366284, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 3.5, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')

## Does zooplankton body size vary with temperature?  [week 8]
## figures 

hist(dataz$community.size, breaks = 20)
hist(1/(dataz$community.size), breaks = 20)
hist(log(dataz$community.size), breaks = 20)

plot(dataz$community.size~dataz$Tank, pch = 19, col = dataz$trophic.level)
plot(dataz$community.size~dataz$week, pch = 19, col = dataz$trophic.level)
plot(dataz$community.size~dataz$invT, pch = 19, col = dataz$trophic.level)

## analysis
modCS0<-lme(log(zoo.ug.carbon.liter+1)~1, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modCS1<-lme(log(zoo.ug.carbon.liter+1)~1+invT, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modCS2<-lme(log(zoo.ug.carbon.liter+1)~1+invT+trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modCS4<-lme(log(zoo.ug.carbon.liter+1)~1+invT*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)

model.sel(modCS0, modCS1, modCS2, modCS4)


model2CS0 <- update(modCS0, correlation = corAR1())
model2CS1 <- update(modCS1, correlation = corAR1())
model2CS2 <- update(modCS2, correlation = corAR1())
model2CS3 <- update(modCS3, correlation = corAR1())
model2CS4 <- update(modCS4, correlation = corAR1())

model2CS0
model2CS1
model2CS2
model2CS3
model2CS4

summary(modCS3)
modCS3<-lm(log(community.size)~1+invT + week + trophic.level, data=data, na.action=na.omit)
summary(modDCR3)

## add lines to plot
abline(-14.551873, 0.361328, lty = 2, lwd = 3, col = 'brown')
abline(-14.745705,0.361328 , lty = 3, lwd = 3, col = 'blue')
legend(40.1, -0.5, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')


# temporal change in mature daphnia size
#figures

hist(data$Daphnia.mature.size)
hist(log(data$Daphnia.mature.size))
hist(log(data$Daphnia.mature.size+1))

week <-data[-which(data$Tank=='18'),]
#plot(log(data$Daphnia.mature.size)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$Daphnia.mature.size)~data$invT)
plot(log(data$Daphnia.mature.size)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(0,0.8), xlim=c(38.0,41),  xlab='inv Temperature (C)', ylab='adult Daphnia size (mm)') 
axis(1, at=c(38.0, 38.5, 39.0, 39.5, 40.0,40.5, 41.0), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0, 0, 0.2, 0.4, 0.6, 0.8), pos=38.0, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='PZ'),]$Daphnia.mature.size)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$Daphnia.mature.size)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modCS0<-lme(log(community.size)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modCS1<-lme(log(community.size)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modCS2<-lme(log(community.size)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modCS3<-lme(log(community.size)~1+invT + week + trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modCS4<-lme(log(community.size)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modCS0, modCS1)
anova(modCS1, modCS2)
anova(modCS2, modCS3)
anova(modCS3, modCS4)

AIC(modDCR0, modCS1, modCS2, modCS3, modCS4)
AICs <- as.data.frame(cbind(AICc(modCS0),AICc(modCS1), AICc(modCS2), AICc(modCS3), AICc(modCS4)))
akaike.weights(AICs)

logLik(modCS0)
logLik(modCS1)
logLik(modCS2)
logLik(modCS3)
logLik(modCS4)

model2CS0 <- update(modCS0, correlation = corAR1())
model2CS1 <- update(modCS1, correlation = corAR1())
model2CS2 <- update(modCS2, correlation = corAR1())
model2CS3 <- update(modCS3, correlation = corAR1())
model2CS4 <- update(modCS4, correlation = corAR1())

model2CS0
model2CS1
model2CS2
model2CS3
model2CS4

summary(modCS3)
modCS3<-lm(log(community.size)~1+invT + week + trophic.level, data=data, na.action=na.omit)
summary(modDCR3)

## add lines to plot
abline(-14.551873, 0.361328, lty = 2, lwd = 3, col = 'brown')
abline(-14.745705,0.361328 , lty = 3, lwd = 3, col = 'blue')
legend(40.1, -0.5, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')



## Does adult zooplankton density vary with temperature?  
## figures 
hist(data$total.adults)
hist(log(data$total.adults))
hist(log(data$total.adults+1))

plot(log(data$total.adults)~data$invT)
plot(log(data$total.adults)~data$invT, cex=1.5, pch='',  axes=FALSE,  ylim=c(-2.0,1.5), xlim=c(38.0,41),  xlab='inv Temperature (C)', ylab='ZP adult density ln(ind / L)') 
axis(1, at=c(38.0, 38.5,39, 39.5, 40,40.5, 41), pos=-2.0, lwd=2, cex.lab=1.5)
axis(2, at=c(-2.0, -1.5, -1.0, -0.5, 0.0, 0.5,1.0, 1.5), pos=38.0, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='PZ'),]$total.adults)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$total.adults)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modad0<-lme(log(total.adults+1)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modad1<-lme(log(total.adults+1)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modad2<-lme(log(total.adults+1)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modad3<-lme(log(total.adults+1)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modad4<-lme(log(total.adults+1)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modad0, modad1)
anova(modad1, modad2)
anova(modad2, modad3)
anova(modad3, modad4)

AIC(modad0, modad1, modad2, modad3, modad4)
AICs <- as.data.frame(cbind(AICc(modad0),AICc(modad1), AICc(modad2), AICc(modad3), AICc(modad4)))
akaike.weights(AICs)

logLik(modad0)
logLik(modad1)
logLik(modad2)
logLik(modad3)
logLik(modad4)

model2ad0 <- update(modad0, correlation = corAR1())
model2ad1 <- update(modad1, correlation = corAR1())
model2ad2 <- update(modad2, correlation = corAR1())
model2ad3 <- update(modad3, correlation = corAR1())
model2ad4 <- update(modad4, correlation = corAR1())

model2ad0
model2ad1
model2ad2
model2ad3
model2ad4

summary(modad3)
modad3<-lm(log(total.adults+1)~1+invT+week+trophic.level, data=data, na.action=na.omit)

## add lines to plot
abline(-5.04, 0.13, lty = 2, lwd = 3, col = 'brown')
abline(-5.18, 0.13, lty = 3, lwd = 3, col = 'blue')
legend(40, -0.7, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')


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
modTC0<-lme(log(total.carbon) ~ 1, random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)
modTC1<-lme(log(total.carbon)~1+I(invT-mean(invT)), random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)
modTC2<-lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)
modTC3<-lme(log(total.carbon)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)
modTC4<-lme(log(total.carbon)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)
modTC5<-lme(log(total.carbon)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)
modTC6<-lme(log(total.carbon)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="ML", na.action=na.omit)

model.sel(modTC0, modTC1, modTC2, modTC3, modTC4, modTC5, modTC6)

anova(modTC4, modTC5)
anova(modTC5, modTC3)
anova(modTC3, modTC6)
anova(modTC6, modTC2)
anova(modTC2, modTC1)
anova(modTC1, modTC0)

### refit with REML and average
modTC3r<-lme(log(total.carbon)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="REML", na.action=na.omit)
modTC4r<-lme(log(total.carbon)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="REML", na.action=na.omit)
modTC5r<-lme(log(total.carbon)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data[data$week >= '4',], method="REML", na.action=na.omit)

m.avg <- model.avg(modTC3r, modTC4r, modTC5r)
summary(m.avg)


#### figures ####
## Biomass figure: 

par(mfrow=c(3,1))
plot(log(datap$PP.biomass)~datap$invT, pch = as.numeric(datap$trophic.level)+14, col = (as.numeric(datap$trophic.level)), ylab = 'Bpp', ylim = c(0, 7))
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$invT, pch = as.numeric(dataz$trophic.level)+14, col = (as.numeric(dataz$trophic.level)), ylab = 'Bzp', ylim = c(0, 7))
plot(log(data[data$week >= '4',]$total.carbon)~data[data$week >= '4',]$invT, pch = as.numeric(data$trophic.level)+14, col = (as.numeric(data$trophic.level)), ylab = 'Bt', ylim = c(0, 7))

par(mfrow=c(1,1)) 
plot(log(datap[datap$week==8,]$PP.biomass)~datap[datap$week==8,]$trophic.level, pch = 19, col = datap$trophic.level, ylab = 'Bpp', ylim = c(0, 7))
summary(aov(log(datap[datap$week==8,]$PP.biomass)~datap[datap$week==8,]$trophic.level))

par(mfrow=c(3,2))
plot(log(data$calc.NPP)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)
plot(log(data$NPP.mass)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)
plot(log(data$calc.ER)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)
plot(log(data$ER.mass)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)
plot(log(data$NEM)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)


plot(log(data[data$trophic.level=='PZ',]$NPP.mass)~data[data$trophic.level=='PZ',]$invT, pch = 19, col = data[data$trophic.level=='PZ',]$trophic.level)




#### messing around with ER data, after reading Yvon Durocher et al 2012.
#### Yvon Durocher and Allen 2012 Figure 1: 



hist(data$ER2)
hist(log(data$ER2))
#hist(log(data$calc.ER+1))

plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$invT, pch = 19, col = data$trophic.level)
mod1 <- lm(log(data$ER2)~data$invT)
plot(log(I(data$ER2/mean(data$ER2)))~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)


## analysis
modER0<-lme(log(ER2) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER1<-lme(log(ER2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER2<-lme(log(ER2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER3<-lme(log(ER2)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER4<-lme(log(ER2)~1+I(invT-mean(invT))*trophic.level+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER5<-lme(log(ER2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
#modER6<-lme(log(ER2)~1+I(invT-mean(invT))*week*trophic.level, random=~1+I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modER0, modER1, modER2, modER3, modER4, modER5)

modER5r<-lme(log(ER2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
summary(modER5r)



## analysis # redo, looks like don't need ranef for temp. so that means single activation energy across tanks. great. result doesn't change if week is a fixed effect.
modNPP0<-lme(log(NPP2) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corCAR1(0.2, form = ~week|Tank))
modNPP1<-lme(log(NPP2)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corCAR1(0.2, form = ~week|Tank))
modNPP2<-lme(log(NPP2)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corCAR1(0.2, form = ~week|Tank))
#modNPP3<-lme(log(NPP2)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP4<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corCAR1(0.2, form = ~week|Tank))
#modNPP5<-lme(log(NPP2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
#modNPP6<-lme(log(NPP2)~1+I(invT-mean(invT))*week*trophic.level, random=~1+I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit)

modNPP4<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|Tank, data=data, method="ML", na.action=na.omit, correlation = corCAR1(0.2, form = ~week|Tank))  #I think this is it. the repeated observations within Tank come from the weekly measurements. We're losing a temporal trend here... but maybe that's ok. 
modNPP4b<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corCAR1(0.2, form = ~week|Tank)) 
modNPP4c<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)

#random=list(~1|Tank, ~1+I(invT-mean(invT))|week)

model.sel(modNPP0, modNPP1, modNPP2, modNPP4)
modNPP3r<-lme(log(NPP2)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
modNPP4r<-lme(log(NPP2)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
modNPP6r<-lme(log(NPP2)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
modNPP5r<-lme(log(NPP2)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)
summary(modNPP5r)

m.avg <- model.avg(modNPP3r, modNPP5r)
summary(m.avg)
