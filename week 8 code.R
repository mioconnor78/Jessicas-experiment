### Mary's analysis of Jessica's experimet
### November 9, 2014
### for WSN talk

### set working directory and load data
setwd("/Users/maryo/Documents/temporary files/Jessicas experment")
data <- read.csv("week8.csv")
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
modNPP2<-lm(log(data$NPP)~1+data$invT*data$trophic.level)
anova(modNPP0, modNPP1)
anova(modNPP1, modNPP2)
anova(modNPP0, modNPP2)
AIC(modNPP0, modNPP1, modNPP2)

summary(modNPP2)
coef(modNPP2)

## add lines to plot
abline(coef(modNPP2)[1], coef(modNPP2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPP2)[1]+coef(modNPP2)[3]), (coef(modNPP2)[2]+coef(modNPP2)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPP2)[1]+coef(modNPP2)[4]), (coef(modNPP2)[2]+coef(modNPP2)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')



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
modNPPm2<-lm(log(data$NPP.mass)~1+data$invT*data$trophic.level)
anova(modNPPm0, modNPPm1)
anova(modNPPm1, modNPPm2)
summary(modNPPm1)
coef(modNPPm1)

## add lines to plot
abline(coef(modNPPm1)[1], coef(modNPPm1)[2], lty = 1, lwd = 3, col = 'black')
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
modER2<-lm(log(data$ER)~1+data$invT*data$trophic.level)
anova(modER0, modER1)
anova(modER1, modER2)
summary(modER2)
coef(modER2)

## add lines to plot
abline(coef(modER2)[1], coef(modER2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modER2)[1]+coef(modER2)[3]), (coef(modER2)[2]+coef(modER2)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modER2)[1]+coef(modER2)[4]), (coef(modER2)[2]+coef(modER2)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


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
modERm2<-lm(log(data$ER.mass)~1+data$invT*data$trophic.level)
anova(modERm0, modERm1)
anova(modERm1, modERm2)
summary(modERm2)
coef(modERm2)

## add lines to plot
abline(coef(modERm2)[1], coef(modERm2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modERm2)[1]+coef(modERm2)[3]), (coef(modERm2)[2]+coef(modERm2)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modERm2)[1]+coef(modERm2)[4]), (coef(modERm2)[2]+coef(modERm2)[6]), lty = 3, lwd = 3, col = 'blue')
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
modchl2<-lm(log(data$chla)~1+data$invT*data$trophic.level)
anova(modchl0, modchl1)
anova(modchl1, modchl2)
summary(modchl1)
coef(modERm2)

## add lines to plot
abline(coef(modchl1)[1], coef(modchl1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, -1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')
