### Mary's analysis of Jessica's experimet
### November 9, 2014
### for WSN talk

### set working directory and load data
setwd("/Users/maryo/Documents/temporary files/Jessicas experment")
data <- read.csv("week8.csv")
dim(data)
data <- data[-31,] # not sure what that last row was; ask Jessica
head(data)


### load libraries, define variables and add columns
library(qpcR)
library(nlme)
k <- 8.62*10^-5
data$invT <-  1/((data$average.temp + 273)*k)
data$invT <- as.numeric(as.character(data$invT))

## Does NPP vary with temperature?  
## figures
hist(data$NPP)
plot(log(data$NPP)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-3,2), xlim=c(38.5,41), xlab='inv Temperature (C)', ylab='NPP ln(mg/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$NPP)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$NPP)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$NPP)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modNPP0<-lm(log(data$NPP)~1)
modNPP1<-lm(log(data$NPP)~1+data$invT)
modNPP2<-lm(log(data$NPP)~1+data$invT+data$trophic.level)
modNPP3<-lm(log(data$NPP)~1+data$invT*data$trophic.level)
anova(modNPP0, modNPP1)
anova(modNPP0, modNPP2)
anova(modNPP2, modNPP3)
AIC(modNPP0, modNPP1, modNPP2)

summary(modNPP2)
coef(modNPP2)
confint(modNPP2)

## add lines to plot
abline(coef(modNPP2)[1], coef(modNPP2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPP2)[1]+coef(modNPP2)[3]), coef(modNPP2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPP2)[1]+coef(modNPP2)[4]), coef(modNPP2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))



## Does mass-specific NPP vary with temperature?  
## figures 
hist(data$NPP.mass)
plot(log(data$NPP.mass)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-8,2), xlim=c(38.5,41), xlab='inv Temperature (C)', ylab='NPP ln(mg O/gC/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-8, lwd=2, cex.lab=1.5)
axis(2, at=c(-8,-6,-4,-2,0,2), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$NPP.mass)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$NPP.mass)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$NPP.mass)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modNPPm0<-lm(log(data$NPP.mass)~1)
modNPPm1<-lm(log(data$NPP.mass)~1+data$invT)
modNPPm2<-lm(log(data$NPP.mass)~1+data$invT+data$trophic.level)
modNPPm3<-lm(log(data$NPP.mass)~1+data$invT*data$trophic.level)
anova(modNPPm0, modNPPm1)
anova(modNPPm1, modNPPm2)
anova(modNPPm2, modNPPm3)
summary(modNPPm2)
coef(modNPPm2)
confint(modNPPm2)

## add lines to plot
abline(coef(modNPPm2)[1], coef(modNPPm2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPPm2)[1]+coef(modNPPm2)[3]), coef(modNPPm2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPPm2)[1]+coef(modNPPm2)[4]), coef(modNPPm2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')



## Does ER vary with temperature?  
## figures 
hist(data$ER)
plot(log(data$ER)~data$invT, cex=1.5, pch='',  axes=FALSE,ylim=c(-3,2), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ER ln(mg O/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$ER)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$ER)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$ER)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')
## analysis
modER0<-lm(log(data$ER)~1)
modER1<-lm(log(data$ER)~1+data$invT)
modER2<-lm(log(data$ER)~1+data$invT+data$trophic.level)
modER3<-lm(log(data$ER)~1+data$invT*data$trophic.level)
anova(modER0, modER1)
anova(modER1, modER2)
anova(modER2, modER3)
summary(modER2)
coef(modER2)
confint(modER2)

## add lines to plot
abline(coef(modER2)[1], coef(modER2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modER2)[1]+coef(modER2)[3]), coef(modER2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modER2)[1]+coef(modER2)[4]), coef(modER2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.7, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))


## Does mass specific ER vary with temperature?  
## figures 
data$ER.mass <- data$ER/data$total.carbon
hist(data$ER.mass)
plot(log(data$ER.mass)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(-8,2), xlab='inv Temperature (C)', ylab='ER ln(mg O/g C/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-8, lwd=2, cex.lab=1.5)
axis(2, at=c(-8,-6,-4,-2,0,2), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$ER.mass)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$ER.mass)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$ER.mass)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modERm0<-lm(log(data$ER.mass)~1)
modERm1<-lm(log(data$ER.mass)~1+data$invT)
modERm2<-lm(log(data$ER.mass)~1+data$invT+data$trophic.level)
modERm3<-lm(log(data$ER.mass)~1+data$invT*data$trophic.level)
anova(modERm0, modERm1)
anova(modERm1, modERm2)
anova(modERm2, modERm3)
summary(modERm2)
coef(modERm2)
confint(modERm2)

## add lines to plot
abline(coef(modERm2)[1], coef(modERm2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modERm2)[1]+coef(modERm2)[3]), coef(modERm2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modERm2)[1]+coef(modERm2)[4]), coef(modERm2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

## Does chla vary with temperature?  
## figures 
hist(data$chla)
plot(log(data$chla)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(-3,3), xlab='inv Temperature (C)', ylab='Chl a ln(ug Chla / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$chla)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$chla)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$chla)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modchl0<-lm(log(data$chla)~1)
modchl1<-lm(log(data$chla)~1+data$invT)
modchl2<-lm(log(data$chla)~1+data$invT+data$trophic.level)
modchl3<-lm(log(data$chla)~1+data$invT*data$trophic.level)
anova(modchl0, modchl1)
anova(modchl1, modchl2)
anova(modchl1, modchl3)
summary(modchl1)
confint(modchl1)

## add lines to plot
abline(coef(modchl1)[1], coef(modchl1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, -1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

## Does total PP biomass vary with temperature and FCL?  
## figures 
hist(data$PP.biomass)
plot(log(data$PP.biomass)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(1,6), xlab='inv Temperature (C)', ylab='PP biomass ln(ug C / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=1, lwd=2, cex.lab=1.5)
axis(2, at=c(1,2,3,4,5,6), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$PP.biomass)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$PP.biomass)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$PP.biomass)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modPb0<-lm(log(data$PP.biomass)~1)
modPb1<-lm(log(data$PP.biomass)~1+data$invT)
modPb2<-lm(log(data$PP.biomass)~1+data$invT+data$trophic.level)
modPb3<-lm(log(data$PP.biomass)~1+data$invT*data$trophic.level)
anova(modPb0, modPb1)
anova(modPb1, modPb2)
anova(modPb1, modPb3)
summary(modPb1)
confint(modPb1)

## add lines to plot
abline(coef(modPb1)[1], coef(modPb1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, -1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')



## Does zooplankton carbon vary with temperature?  
## figures 
hist(data$zooplankton.carbon.per.L)
plot(log(data$zooplankton.carbon.per.L)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-1,3), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ZP biomass ln(g C / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-1, lwd=2, cex.lab=1.5)
axis(2, at=c(-1,0,1,2,3), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='PZ'),]$zooplankton.carbon.per.L)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$zooplankton.carbon.per.L)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modzpc0<-lm(log(data$zooplankton.carbon.per.L)~1)
modzpc1<-lm(log(data$zooplankton.carbon.per.L)~1+data$invT)
modzpc2<-lm(log(data$zooplankton.carbon.per.L)~1+data$invT+data$trophic.level)
modzpc3<-lm(log(data$zooplankton.carbon.per.L)~1+data$invT*data$trophic.level)
anova(modzpc0, modzpc1)
anova(modzpc1, modzpc2)
anova(modzpc1, modzpc3)
summary(modzpc3)
confint(modzpc3)

## add lines to plot
abline(coef(modzpc3)[1], coef(modzpc3)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modzpc3)[1]+coef(modzpc3)[3]), (coef(modzpc3)[2]+coef(modzpc3)[4]), lty = 3, lwd = 3, col = 'blue')
legend(38.5, 3, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')

## Does zooplankton density vary with temperature?  
## figures 
hist(data$ZP.per.L)
plot(log(data$ZP.per.L)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-2,4), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ZP density ln(ind / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-2, lwd=2, cex.lab=1.5)
axis(2, at=c(-2,-1,0,1,2,3,4), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='PZ'),]$ZP.per.L)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$ZP.per.L)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modzp0<-lm(log(data$ZP.per.L)~1)
modzp1<-lm(log(data$ZP.per.L)~1+data$invT)
modzp2<-lm(log(data$ZP.per.L)~1+data$invT+data$trophic.level)
modzp3<-lm(log(data$ZP.per.L)~1+data$invT*data$trophic.level)
anova(modzp0, modzp1)
anova(modzp1, modzp2)
anova(modzp1, modzp3)
summary(modzp3)

## add lines to plot
abline(coef(modzp3)[1], coef(modzp3)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modzp3)[1]+coef(modzp3)[3]), (coef(modzp3)[2]+coef(modzp3)[4]), lty = 3, lwd = 3, col = 'blue')
legend(38.5, 4, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')

## Does total biomass vary with temperature and trohpic structure?  
## figures 
hist(data$total.carbon)
plot(log(data$total.carbon)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(0,6), xlab='inv Temperature (C)', ylab='Biomass ln(g C/L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,2,4,6), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$total.carbon)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$total.carbon)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$total.carbon)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modTCm0<-lm(log(data$total.carbon)~1)
modTCm1<-lm(log(data$total.carbon)~1+data$invT)
modTCm2<-lm(log(data$total.carbon)~1+data$invT+data$trophic.level)
modTCm3<-lm(log(data$total.carbon)~1+data$invT*data$trophic.level)
anova(modTCm0, modTCm1)
anova(modTCm1, modTCm2)
anova(modTCm1, modTCm3)
summary(modTCm1)
confint(modTCm1)

## add lines to plot
abline(coef(modTCm1)[1], coef(modTCm1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, 3, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))

