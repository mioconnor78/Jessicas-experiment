
### load libraries
#library(qpcR)
library(nlme)
library(MuMIn)

### set working directory and load data
data <- read.csv("./temporal_dataFEB12.csv")
dim(data)
head(data)
tail(data)

### load libraries, define variables and add columns
k <- 8.617342*10^-5  # eV/K
data$invT <-  1/((data$average.temp + 273)*k)
data$invT <- as.numeric(as.character(data$invT))
#data$kT <-  ((data$average.temp + 273)*k)
#data$kT <- as.numeric(as.character(data$kT))*100

data$PP.biomass <- (data$chla*55) #chla (ug/L)* 55 C in PP / 1 chla = ugPPC/L
data$total.carbon <- data$PP.biomass + data$zoo.ug.carbon.liter #I'm pretty sure zp was in ugC/L
data$NPP.mass <- data$calc.NPP / (data$PP.biomass)
data$ER.mass <- data$calc.ER/(data$total.carbon)

## temporal

## Does NPP vary with temperature?  
## figures on invT
hist(data$calc.NPP)
hist(log(data$calc.NPP))
hist(log(data$calc.NPP+1))

plot(log(data$calc.NPP+1)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$calc.NPP+1)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$calc.NPP+1)~data$invT, pch = 19, col = data$trophic.level)

## analysis
modNPP0<-lme(log(calc.NPP+1) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP1<-lme(log(calc.NPP+1)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP2<-lme(log(calc.NPP+1)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP3<-lme(log(calc.NPP+1)~1+I(invT-mean(invT))+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP4<-lme(log(calc.NPP+1)~1+week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP5<-lme(log(calc.NPP+1)~1+I(invT-mean(invT))*week+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPP6<-lme(log(calc.NPP+1)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modNPP0, modNPP1, modNPP2, modNPP3, modNPP4, modNPP5, modNPP6)

# for model fitting: 
modNPP6<-lme(log(calc.NPP+1)~1+I(invT-mean(invT))*week*trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit)

model2modNPP0 <- update(modNPP0, correlation = corAR1())
model2modNPP1 <- update(modNPP1, correlation = corAR1())
model2modNPP2 <- update(modNPP2, correlation = corAR1())
model2modNPP3 <- update(modNPP3, correlation = corAR1())
model2modNPP4 <- update(modNPP4, correlation = corAR1())

model2modNPP0
model2modNPP1
model2modNPP2
model2modNPP3
model2modNPP4


summary(modNPP4)
coef(modNPP4)
confint(modNPP4)

## add lines to plot
abline(29.941985, -0.752594, lty = 1, lwd = 3, col = 'seagreen')
abline(34.32669, -0.862494, lty = 2, lwd = 3, col = 'brown')
abline(31.138628, -0.779939, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty='n')



# net ecosystem metabolism
data$NEM <- 18*data$NPP - 24*data$ER  #this is not a thing; because ER is already part of NPP. Could
# look at the two over a 24 hour period... so ER*24 but NPP*18...

hist(data$NEM)
hist(sqrt(data$NEM))
hist(log(data$NEM))
hist(log(data$NEM+1))
plot(log(data$NEM)~data$invT, pch = 19, col = data$trophic.level)

plot(log(data$NEM)~data$invT, pch = 19, col = data$trophic.level)

plot(log(data$NEM)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-0.5,2.0), xlim=c(37.5,41), xlab='inv(Temperature) 1/eV', ylab='NEM ln(mg O/L/hr)') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40,40.5, 41), pos=-0.5, lwd=2, cex.lab=1.5)
axis(2, at=c(-0.5,0,0.5,1,1.5,2.0), pos=37.5, lwd=2, cex.lab=1.5)
abline(1, 0, lwd = 3, col = 1, lty = 2)
points(log(3+data[(data$trophic.level=='P'),]$NEM)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(3+data[(data$trophic.level=='PZ'),]$NEM)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(3+data[(data$trophic.level=='PZN'),]$NEM)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')


## analysis
modNEM0<-lme(log(NEM)~1, random=~1|Tank, data=data, method="ML", na.action=na.exclude)
modNEM1<-lme(log(3+NEM)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNEM2<-lme(log(3+NEM)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNEM3<-lme(log(3+NEM)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNEM4<-lme(log(3+NEM)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modNEM0, modNEM1)
anova(modNEM0, modNEM2)
anova(modNEM2, modNEM3)
anova(modNEM3, modNEM4)

AIC(modNEM0, modNEM1, modNEM2, modNEM3, modNEM4)
AICs <- as.data.frame(cbind(AICc(modNEM0),AICc(modNEM1), AICc(modNEM2), AICc(modNEM3), AICc(modNEM4)))
akaike.weights(AICs)

logLik(modNEM0)
logLik(modNEM1)
logLik(modNEM2)
logLik(modNEM3)
logLik(modNEM4)

model2modNEM0 <- update(modNEM0, correlation = corAR1())
model2modNEM1 <- update(modNEM1, correlation = corAR1())
model2modNEM2 <- update(modNEM2, correlation = corAR1())
model2modNEM3 <- update(modNEM3, correlation = corAR1())
model2modNEM4 <- update(modNEM4, correlation = corAR1())

model2modNEM0
model2modNEM1
model2modNEM2
model2modNEM3
model2modNEM4

summary(modNEM0)
coef(modNEM3)
confint(modNEM2)

## add lines to plot
abline(-2.74, 0.10, lty = 1, lwd = 3, col = 'seagreen')
abline(-2.83, 0.10, lty = 2, lwd = 3, col = 'brown')
abline(-2.76, 0.10, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty='n')


## Does mass-specific NPP vary with temperature?  
## figures  
hist(data$NPP.mass)
hist(log(data$NPP.mass))
hist(log(data$NPP.mass+1))

plot(log(data$NPP.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$chla)~data$Tank, col = data$trophic.level)
#data1 <- data[-which(data$Tank=='30'),]

plot(log(data$NPP.mass)~data$invT)
plot(log(data$NPP.mass)~data$invT, cex=1.5, pch='', ylim=c(-10,1),  axes=FALSE, xlim=c(37.5,41), xlab='inv Temperature (C)', ylab='NPP ln(mg O/gC/L/hr)') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40,40.5, 41), pos=-10, lwd=2, cex.lab=1.5)
axis(2, at=c(-10, -8,-6,-4,-2, 0,1), pos=37.5, lwd=2, cex.lab=1.5)
abline(10, -0.32, lwd = 3, col = 2)
points(log(data[(data$trophic.level=='P'),]$NPP.mass)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$NPP.mass)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$NPP.mass)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modNPPm0<-lme(log(NPP.mass)~1, random=~1|Tank,data=data , method="ML", na.action=na.omit)
modNPPm1<-lme(log(NPP.mass)~1+invT, random=~1|Tank,data=data , method="ML", na.action=na.omit)
modNPPm2<-lme(log(NPP.mass)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm3<-lme(log(NPP.mass)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modNPPm4<-lme(log(NPP.mass)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modNPPm0, modNPPm1)
anova(modNPPm1, modNPPm2)
anova(modNPPm2, modNPPm3)
anova(modNPPm3, modNPPm4)

AIC(modNPPm0, modNPPm1, modNPPm2, modNPPm3, modNPPm4)
AICs <- as.data.frame(cbind(AICc(modNPPm0),AICc(modNPPm1), AICc(modNPPm2), AICc(modNPPm3), AICc(modNPPm4)))
akaike.weights(AICs)

logLik(modNPPm0)
logLik(modNPPm1)
logLik(modNPPm2)
logLik(modNPPm3)
logLik(modNPPm4)

model2modNPPm0 <- update(modNPPm0, correlation = corAR1())
model2modNPPm1 <- update(modNPPm1, correlation = corAR1())
model2modNPPm2 <- update(modNPPm2, correlation = corAR1())
model2modNPPm3 <- update(modNPPm3, correlation = corAR1())
model2modNPPm4 <- update(modNPPm4, correlation = corAR1())

model2modNPPm0
model2modNPPm1
model2modNPPm2
model2modNPPm3
model2modNPPm4

summary(modNPPm4)
coef(modNPPm2)
confint(modNPPm2)

## add lines to plot for modNPPm3
abline(97.78111, -2.71040, lty = 1, lwd = 3, col = 'seagreen')
abline(127.4316, -3.44018, lty = 2, lwd = 3, col = 'brown')
abline(96.45909, -2.66458, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 1.0, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty='n')

## add lines to plot for modNPPm2
abline(coef(modNPPm2)[1], coef(modNPPm2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPPm2)[1]+coef(modNPPm2)[3]), coef(modNPPm2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPPm2)[1]+coef(modNPPm2)[4]), coef(modNPPm2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.0, 0, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), b='n')



## Does ER vary with temperature?  
## figures 
hist(data$calc.ER)
hist(log(data$calc.ER))
hist(log(data$calc.ER+1))

plot(log(data$calc.ER)~data$invT)
plot(log(data$calc.ER)~data$invT, cex=1.5, pch='',  axes=FALSE,ylim=c(-3,3), xlim=c(37.5,41),  xlab='inv Temperature (C)', ylab='ER ln(mg O/L/hr)') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40, 40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2, 3), pos=37.5, lwd=2, cex.lab=1.5)
abline(26, -0.65, lwd = 3, col = 2)
points(log(data[(data$trophic.level=='P'),]$calc.ER)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$calc.ER)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$calc.ER)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modER0<-lme(log(calc.ER)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER1<-lme(log(calc.ER)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER2<-lme(log(calc.ER)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER3<-lme(log(calc.ER)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modER4<-lme(log(calc.ER)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modER0, modER1)
anova(modER1, modER2)
anova(modER2, modER3)
anova(modER3, modER4)

AIC(modER0, modER1, modER2, modER3, modER4)
AICs <- as.data.frame(cbind(AICc(modER0),AICc(modER1), AICc(modER2), AICc(modER3), AICc(modER4)))
akaike.weights(AICs)

logLik(modER0)
logLik(modER1)
logLik(modER2)
logLik(modER3)
logLik(modER4)

model2modER0 <- update(modER0, correlation = corAR1())
model2modER1 <- update(modER1, correlation = corAR1())
model2modER2 <- update(modER2, correlation = corAR1())
model2modER3 <- update(modER3, correlation = corAR1())
model2modER4 <- update(modER4, correlation = corAR1())

model2modER0
model2modER1
model2modER2
model2modER3
model2modER4


summary(modER4)
coef(modER2)
confint(modER2)

## add lines to plot
abline(52.37740,-1.34920  , lty = 1, lwd = 3, col = 'seagreen')
abline(42.48359,-1.0836, lty = 2, lwd = 3, col = 'brown')
abline(30.04742,-0.77776, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 3.5, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


## Does mass specific ER vary with temperature?  
## figures 
hist(data$ER.mass)
hist(log(data$ER.mass))
hist(log(data$ER.mass+1))

plot(log(data$ER.mass)~data$invT)
plot(log(data$ER.mass)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(37.5,41), ylim=c(-8.5,1), xlab='inv Temperature (C)', ylab='ER ln(mg O/g C/L/hr)') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40, 40.5, 41), pos=-8.5, lwd=2, cex.lab=1.5)
axis(2, at=c(-8.5,-8,-6,-4,-2, 0, 1), pos=37.5, lwd=2, cex.lab=1.5)
abline(89, -0.65, lwd = 3, col = 2)
points(log(data[(data$trophic.level=='P'),]$ER.mass)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$ER.mass)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$ER.mass)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modERm0<-lme(log(ER.mass)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm1<-lme(log(ER.mass)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm2<-lme(log(ER.mass)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm3<-lme(log(ER.mass)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modERm4<-lme(log(ER.mass)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modERm0, modERm1)
anova(modERm1, modERm2)
anova(modERm2, modERm3)
anova(modERm3, modERm4)

AIC(modERm0, modERm1, modERm2, modERm3, modERm4)
AICs <- as.data.frame(cbind(AICc(modERm0),AICc(modERm1), AICc(modERm2), AICc(modERm3), AICc(modERm4)))
akaike.weights(AICs)

logLik(modERm0)
logLik(modERm1)
logLik(modERm2)
logLik(modERm3)
logLik(modERm4)

model2ERm0 <- update(modERm0, correlation = corAR1())
model2ERm1 <- update(modERm1, correlation = corAR1())
model2ERm2 <- update(modERm2, correlation = corAR1())
model2ERm3 <- update(modERm3, correlation = corAR1())
model2ERm4 <- update(modERm4, correlation = corAR1())

model2ERm0
model2ERm1
model2ERm2
model2ERm3
model2ERm4

summary(modERm4)
coef(modERm2)
confint(modERm2)

## add lines to plot
abline(109.44036,-2.95643, lty = 1, lwd = 3, col = 'seagreen')
abline(133.83687,-3.54009, lty = 2, lwd = 3, col = 'brown')
abline(89.73642,-2.4221, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

##### BIOMASS RESULTS #####

## Does chla vary with temperature?  
### MO: went through this and made minor changes on July 28, 2015
## figures 
# data1 <- data[(which(data$week!='2')),] #could redo without week2, to see if that's driving the need for week effects.
hist(data$chla)
hist(log(data$chla))
hist(log(data$chla+1))
plot(log(data$chla+1)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$chla+1)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$chla+1)~data$invT, pch = 19, col = data$trophic.level)

## analysis
modchl0<-lme(log(chla+1)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl1<-lme(log(chla+1)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl2<-lme(log(chla+1)~1+invT+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl3<-lme(log(chla+1)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl4<-lme(log(chla+1)~1+week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl5<-lme(log(chla+1)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modchl6<-lme(log(chla+1)~1+invT*week*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

model.sel(modchl0, modchl1, modchl2, modchl3, modchl4, modchl5, modchl6)

anova(modchl6, modchl5)
anova(modchl5, modchl4)
anova(modchl4, modchl3)
anova(modchl3, modchl2)
anova(modchl2, modchl1)
anova(modchl1,modchl0)

## refit best model with reml
modchl6b<-lme(log(chla+1)~1+invT*week*trophic.level, random=~1|Tank, data=data1, method="REML", na.action=na.omit)
summary(modchl6b)


logLik(modchl0)
logLik(modchl1)
logLik(modchl2)
logLik(modchl3)
logLik(modchl4)

model2chl0 <- update(modchl0, correlation = corAR1())
model2chl1 <- update(modchl1, correlation = corAR1())
model2chl2 <- update(modchl2, correlation = corAR1())
model2chl3 <- update(modchl3, correlation = corAR1())
model2chl4 <- update(modchl4, correlation = corAR1())

model2chl0
model2chl1
model2chl2
model2chl3
model2chl4

summary(modchl4)
confint(modchl3)
coef(modchl3)

## add lines to plot
### but, best model now needs different slopes for each week and each tank... so probably not worth adding slopes.
abline(-41.74267,1.12806, lty = 1, lwd = 3, col = 'seagreen')
abline(-51.26174, 1.35398, lty = 2, lwd = 3, col = 'brown')
abline(-41.55714, 1.11951, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 1.5, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


## MO skipped this for now b/c it's repeating chla and we probably shouldn't analyze both.
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))
hist(log(data$PP.biomass+10))
plot(log(data$PP.biomass+10)~data$invT)

plot(log(data$PP.biomass+10)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(37.5,41), ylim=c(2.0, 7.5), xlab='inv Temperature (C)', ylab='PP biomass ln(ug C / L)') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40, 40.5, 41), pos=2.0, lwd=2, cex.lab=1.5)
axis(2, at=c(2.0, 3.0, 4.0, 5.0, 6.0, 7.0), pos=37.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$PP.biomass+10)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$PP.biomass+10)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$PP.biomass+10)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modPb0<-lme(log(PP.biomass)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modPb1<-lme(log(PP.biomass)~1 + invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modPb2<-lme(log(PP.biomass)~1 + invT + week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modPb3<-lme(log(PP.biomass)~1 + invT + week + trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modPb4<-lme(log(PP.biomass)~1 + invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modPb0, modPb1)
anova(modPb1, modPb2)
anova(modPb3, modPb4)
anova(modPb4, modPb1)
anova(modPb4, modPb0)
anova(modPb4, modPb2)

AIC(modPb0, modPb1, modPb2, modPb3, modPb4)
AICs <- as.data.frame(cbind(AICc(modPb0),AICc(modPb1), AICc(modPb2), AICc(modPb3), AICc(modPb4)))
akaike.weights(AICs)

logLik(modPb0)
logLik(modPb1)
logLik(modPb2)
logLik(modPb3)
logLik(modPb4)

model2Pb0 <- update(modPb0, correlation = corAR1())
model2Pb1 <- update(modPb1, correlation = corAR1())
model2Pb2 <- update(modPb2, correlation = corAR1())
model2Pb3 <- update(modPb3, correlation = corAR1())
model2Pb4 <- update(modPb4, correlation = corAR1())

model2Pb0
model2Pb1
model2Pb2
model2Pb3
model2Pb4
summary(modPb4)
confint(modPb4)

## add lines to plot
abline(-28.604018, 0.892368, lty = 1, lwd = 3, col = 'seagreen')
abline(-19.364001, 0.629293, lty = 2, lwd = 3, col = 'brown')
abline(-20.586785, 0.681782, lty = 3, lwd = 3, col = 'blue')
legend(40.0, 4.5, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


## Does zooplankton carbon vary with temperature?  
## figures 
dataz <- data[(which(data$trophic.level != 'P')),] # because we don't have observations for the P treatments (we're just assuming they are 0) I don't htink we should analyze those tanks here.
hist(data$zoo.ug.carbon.liter)
hist(log(data$zoo.ug.carbon.liter))
hist(log(data$zoo.ug.carbon.liter+1))
plot(log(data$zoo.ug.carbon.liter)~data$Tank, pch = 19, col = data$trophic.level)

plot(log(dataz$zoo.ug.carbon.liter+1)~dataz$Tank, pch = 19, col = dataz$trophic.level)
plot(log(dataz$zoo.ug.carbon.liter+1)~dataz$week, pch = 19, col = dataz$trophic.level)
plot(log(dataz$zoo.ug.carbon.liter+1)~dataz$invT, pch = 19, col = dataz$trophic.level)

## analysis
modzpc0<-lme(log(zoo.ug.carbon.liter+1)~1, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpc1<-lme(log(zoo.ug.carbon.liter+1)~1+invT, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpc2<-lme(log(zoo.ug.carbon.liter+1)~1+invT+trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpc3<-lme(log(zoo.ug.carbon.liter+1)~1+invT+week+trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpc4<-lme(log(zoo.ug.carbon.liter+1)~1+week+invT*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpc5<-lme(log(zoo.ug.carbon.liter+1)~1+invT*week+invT*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)
modzpc6<-lme(log(zoo.ug.carbon.liter+1)~1+invT*week*trophic.level, random=~1|Tank, data=dataz, method="ML", na.action=na.omit)

model.sel(modzpc0, modzpc1, modzpc2, modzpc3, modzpc4, modzpc5, modzpc6)

anova(modzpc4, modzpc5)
anova(modzpc4, modzpc1)
anova(modzpc1, modzpc6)
anova(modzpc1, modzpc2)
anova(modzpc1, modzpc3)
anova(modzpc1, modzpc0)

summary(modzpc4)

model2zpc0 <- update(modzpc0, correlation = corAR1())
model2zpc1 <- update(modzpc1, correlation = corAR1())
model2zpc2 <- update(modzpc2, correlation = corAR1())
model2zpc3 <- update(modzpc3, correlation = corAR1())
model2zpc4 <- update(modzpc4, correlation = corAR1())

model2zpc0
model2zpc1
model2zpc2
model2zpc3
model2zpc4

summary(modzpc4)
confint(modzpc3)

## add lines to plot
abline(7.04955,-0.14377, lty = 2, lwd = 3, col = 'brown')
abline(-56.5769, 1.47017, lty = 3, lwd = 3, col = 'blue')
legend(40.1, 6.5, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')


## Does zooplankton density vary with temperature?  
## figures 
hist(dataz$total.zoo.abundance.liter, breaks = 40)
hist(log(dataz$total.zoo.abundance.liter), breaks = 40)
hist(log(dataz$total.zoo.abundance.liter+1))

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



## Does total biomass vary with temperature and trophic structure?  
## figures 
hist(data$total.carbon)
hist(log(data$total.carbon))
hist(log(data$total.carbon+1))

plot(log(data$total.carbon+1)~data$invT)
plot(log(data$total.carbon+1)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.0,41), ylim=c(1,7), xlab='inv Temperature (C)', ylab='Biomass ln(ug C/L)') 
axis(1, at=c(38.0, 38.5,39, 39.5, 40, 40.5, 41), pos=1, lwd=2, cex.lab=1.5)
axis(2, at=c(1, 2, 3, 4, 5, 6, 7), pos=38.0, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$total.carbon+1)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$total.carbon+1)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$total.carbon+1)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modTCm0<-lme(log(total.carbon+1)~1, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modTCm1<-lme(log(total.carbon+1)~1+invT, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modTCm2<-lme(log(total.carbon+1)~1+invT+week, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modTCm3<-lme(log(total.carbon+1)~1+invT+week+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)
modTCm4<-lme(log(total.carbon+1)~1+invT*week+invT*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit)

anova(modTCm0, modTCm1)
anova(modTCm1, modTCm2)
anova(modTCm2, modTCm3)
anova(modTCm3, modTCm4)

AIC(modTCm0, modTCm1, modTCm2, modTCm3, modTCm4)
AICs <- as.data.frame(cbind(AICc(modTCm0),AICc(modTCm1), AICc(modTCm2), AICc(modTCm3), AICc(modTCm4)))
akaike.weights(AICs)

logLik(modTCm0)
logLik(modTCm1)
logLik(modTCm2)
logLik(modTCm3)
logLik(modTCm4)

model2TCm0 <- update(modTCm0, correlation = corAR1(form=1|data$week))
model2TCm1 <- update(modTCm1, correlation = corAR1())
model2TCm2 <- update(modTCm2, correlation = corAR1())
model2TCm3 <- update(modTCm3, correlation = corAR1())
model2TCm4 <- update(modTCm4, correlation = corAR1())

model2TCm0
model2TCm1
model2TCm2
model2TCm3
model2TCm4

summary(modTCm1)
confint(modTCm1)

## add lines to plot
abline(coef(modTCm1)[1], coef(modTCm1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, 3, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))


## plotting ZP biomass x PP biomass
plot(data$PP.biomass ~ data$total.carbon)
abline(0, 1)

plot(data$PP.biomass ~ data$zoo.ug.carbon.liter, pch = 19, col = c(data$average.temp))
abline(0, 1)

plot(data$zoo.ug.carbon.liter ~ data$PP.biomass, pch = 19, col = c(data$average.temp))
abline(0, 1)
