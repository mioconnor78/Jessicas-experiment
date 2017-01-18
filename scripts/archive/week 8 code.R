### Mary's analysis of Jessica's experimet
### Feb 9, 2015
### building on WSN talk code for manuscript with complete and revised datafile

### load libraries
library(qpcR)
library(nlme)

### set working directory and load data
data <- read.csv("./temporal_data.csv")
dim(data)
head(data)
tail(data)

### load libraries, define variables and add columns
k <- 8.617342*10^-5  # eV/K
data$invT <-  1/((data$average.temp + 273)*k)
data$invT <- as.numeric(as.character(data$invT))
#data$kT <-  ((data$average.temp + 273)*k)
#data$kT <- as.numeric(as.character(data$kT))*100

## estimating NPP, ER and GPP
# NPP = sum of changes in O2 over the day (but not night; so our dusk-dawn)
data$NPPd <- data$dusk - data$dawn1
# ER = Rday + sum of change in O2 over the night (our dawn2 - dusk); Rday is estimated by extrapolating the mean night time hourly R over the hours of daylight.
data$ERd <- (-(data$dawn2 - data$dusk)/data$hours2)*24
# GPP = NPP + Rday, where Rday is daytime respiration (estimated below)
data$GPP <- data$NPPd + data$ERd
# NEM = oxygen gained / oxygen lost, so NPP.daily / ER.daily
data$NEM <- data$NPPd / data$ERd  # values > 1 = net oxygen producers = net autotrophic 
data$NEM2 <- data$dawn2 - data$dawn1
data$MB <- data$GPP / data$ERd

data$PP.biomass <- (data$chla*55) #chla (ug/L)* 55 C in PP / 1 chla = ugPPC/L
#data$total.carbon <- data$PP.biomass + data$zoo.carbon.liter #I'm pretty sure zp was in ugC/L
data$NPP.mass <- data$NPPd / (data$PP.biomass)
data$totalC <- ifelse(data$trophic.level == 'P', data$PP.biomass, (data$PP.biomass + data$zoo.ug.carbon.liter))
data$ER.mass <- data$ERd/(data$totalC)

## just week 8
week8 <- data[which(data$week == '8'),]

## Does NPP vary with temperature?  
## figures on invT
hist(week8$NPPd)
plot(log(week8$NPPd)~week8$Tank, pch = 19, col = week8$trophic.level)
plot(log(week8$NPPd)~week8$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-3,2), xlim=c(38.5,41), xlab='inv(Temperature) 1/eV', ylab='NPP ln(mg O/L/d)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
abline(12, -0.32, lwd = 3, col = 2)
points(log(week8[(week8$trophic.level=='P'),]$NPPd)~week8[(week8$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$NPPd)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$NPPd)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')


## analysis
modNPP0<-lm(log(week8$calc.NPP)~1)
modNPP1<-lm(log(week8$calc.NPP)~1+week8$invT)
modNPP2<-lm(log(week8$calc.NPP)~1+week8$invT+week8$trophic.level)
modNPP3<-lm(log(week8$calc.NPP)~1+week8$invT*week8$trophic.level)
anova(modNPP0, modNPP1)
anova(modNPP1, modNPP2)
anova(modNPP2, modNPP3)
AIC(modNPP0, modNPP1, modNPP2, modNPP3)

summary(modNPP2)
coef(modNPP2)
confint(modNPP2)

## add lines to plot
abline(coef(modNPP2)[1], coef(modNPP2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPP2)[1]+coef(modNPP2)[3]), coef(modNPP2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPP2)[1]+coef(modNPP2)[4]), coef(modNPP2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))



# net ecosystem metabolism
hist(log(week8$NEM))
plot(log(week8$NEM)~week8$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-3,2), xlim=c(38.5,41), xlab='inv(Temperature) 1/eV', ylab='NEM ln(mg O/L/d)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
abline(0, 0, lwd = 3, col = 1, lty = 2)
points(log(week8[(week8$trophic.level=='P'),]$NEM)~week8[(week8$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$NEM)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$NEM)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')


## analysis
modNEM0<-lm(log(week8[(week8$trophic.level!='PZN'),]$NEM)~1)
modNEM1<-lm(log(week8[(week8$trophic.level!='PZN'),]$NEM)~1+week8[(week8$trophic.level!='PZN'),]$invT)
modNEM2<-lm(log(week8[(week8$trophic.level!='PZN'),]$NEM)~1+week8[(week8$trophic.level!='PZN'),]$invT+week8[(week8$trophic.level!='PZN'),]$trophic.level)
modNEM3<-lm(log(week8[(week8$trophic.level!='PZN'),]$NEM)~1+week8[(week8$trophic.level!='PZN'),]$invT*week8[(week8$trophic.level!='PZN'),]$trophic.level)
anova(modNEM0, modNEM1)
anova(modNEM0, modNEM2)
anova(modNEM2, modNEM3)
AIC(modNEM0, modNEM1, modNEM2, modNEM3)

summary(modNEM3)
coef(modNEM3)
confint(modNEM2)

## add lines to plot
abline(coef(modNEM3)[1], coef(modNEM3)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNEM3)[1]+coef(modNEM3)[3]), (coef(modNEM3)[2]+coef(modNEM3)[4]), lty = 2, lwd = 3, col = 'brown')
#abline((coef(modNEM3)[1]+coef(modNEM3)[4]), (coef(modNEM3)[2]+coef(modNEM3)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))


## Does mass-specific NPP vary with temperature?  
## figures  
hist(week8$NPP.mass)
plot(log(week8$NPP.mass)~week8$Tank, pch = 19, col = week8$trophic.level)
plot(log(week8$chla)~week8$Tank, col = week8$trophic.level)
#data1 <- data[-which(data$Tank=='30'),]

plot(log(week8$NPP.mass)~week8$invT, cex=1.5, pch='', ylim=c(-8,2),  axes=FALSE, xlim=c(38.5,41), xlab='inv Temperature (C)', ylab='NPP ln(mg O/gC/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-8, lwd=2, cex.lab=1.5)
axis(2, at=c(-8,-6,-4,-2,0,2), pos=38.5, lwd=2, cex.lab=1.5)
abline(10, -0.32, lwd = 3, col = 2)
points(log(week8[(week8$trophic.level=='P'),]$NPP.mass)~week8[(week8$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$NPP.mass)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$NPP.mass)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modNPPm0<-lm(log(week8$NPP.mass)~1)
modNPPm1<-lm(log(week8$NPP.mass)~1+week8$invT)
modNPPm2<-lm(log(week8$NPP.mass)~1+week8$invT+week8$trophic.level)
modNPPm3<-lm(log(week8$NPP.mass)~1+week8$invT*week8$trophic.level)
anova(modNPPm0, modNPPm1)
anova(modNPPm1, modNPPm2)
anova(modNPPm2, modNPPm3)
AIC(modNPPm0, modNPPm1, modNPPm2, modNPPm3)
summary(modNPPm3)
coef(modNPPm2)
confint(modNPPm2)

## add lines to plot for modNPPm3
abline(coef(modNPPm3)[1], coef(modNPPm3)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPPm3)[1]+coef(modNPPm3)[3]), (coef(modNPPm3)[2]+coef(modNPPm3)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPPm3)[1]+coef(modNPPm3)[4]), (coef(modNPPm3)[2]+coef(modNPPm3)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, -2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))

## add lines to plot for modNPPm2
abline(coef(modNPPm2)[1], coef(modNPPm2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPPm2)[1]+coef(modNPPm2)[3]), coef(modNPPm2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPPm2)[1]+coef(modNPPm2)[4]), coef(modNPPm2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.5, -2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))



## Does ERd vary with temperature?  
## figures 
hist(week8$ERd)  # have some negative values here, which means that some tanks produced oxygen over night.
plot(log(week8$ERd)~week8$invT, cex=1.5, pch='',  axes=FALSE,ylim=c(-3,2), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ER ln(mg O/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
abline(20, -0.65, lwd = 3, col = 2)
points(log(week8[(week8$trophic.level=='P'),]$ERd)~week8[(week8$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$ERd)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$ERd)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modER0<-lm(log(week8$ERd)~1)
modER1<-lm(log(week8$ERd)~1+week8$invT)
modER2<-lm(log(week8$ERd)~1+week8$invT+week8$trophic.level)
modER3<-lm(log(week8$ERd)~1+week8$invT*week8$trophic.level)
anova(modER0, modER1)
anova(modER1, modER2)
anova(modER2, modER3)
AIC(modER0, modER1, modER2, modER3)
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
hist(log(week8$ER.mass))
plot(log(week8$ER.mass)~week8$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(-8,2), xlab='inv Temperature (C)', ylab='ER ln(mg O/g C/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-8, lwd=2, cex.lab=1.5)
axis(2, at=c(-8,-6,-4,-2,0,2), pos=38.5, lwd=2, cex.lab=1.5)
abline(21, -0.65, lwd = 3, col = 2)
points(log(week8[(week8$trophic.level=='P'),]$ER.mass)~week8[(week8$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$ER.mass)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$ER.mass)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modERm0<-lm(log(week8$ER.mass)~1)
modERm1<-lm(log(week8$ER.mass)~1+week8$invT)
modERm2<-lm(log(week8$ER.mass)~1+week8$invT+week8$trophic.level)
modERm3<-lm(log(week8$ER.mass)~1+week8$invT*week8$trophic.level)
anova(modERm0, modERm1)
anova(modERm1, modERm2)
anova(modERm2, modERm3)
AIC(modERm0, modERm1, modERm2, modERm3)

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
hist(week8$chla)
plot(log(week8$chla)~week8$Tank, pch = 19, col = week8$trophic.level)
data1 <- week8
  # week8[-which(week8$Tank=='30'),]
plot(log(data1$chla)~data1$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(-4,3), xlab='inv Temperature (C)', ylab='Chl a ln(ug Chla / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-4, lwd=2, cex.lab=1.5)
axis(2, at=c(-4,-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
abline(-25.5, 0.65, lwd = 3, col = 2)
points(log(data1[(data1$trophic.level=='P'),]$chla)~data1[(data1$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data1[(data1$trophic.level=='PZ'),]$chla)~data1[(data1$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data1[(data1$trophic.level=='PZN'),]$chla)~data1[(data1$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modchl0<-lm(log(data1$chla)~1)
modchl1<-lm(log(data1$chla)~1+data1$invT)
modchl2<-lm(log(data1$chla)~1+data1$invT+data1$trophic.level)
modchl3<-lm(log(data1$chla)~1+data1$invT*data1$trophic.level)
anova(modchl0, modchl1)
anova(modchl1, modchl2)
anova(modchl2, modchl3)
AIC(modchl0, modchl1, modchl2, modchl3)
summary(modchl3)
confint(modchl3)

## add lines to plot
abline(coef(modchl3)[1], coef(modchl3)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modchl3)[1]+coef(modchl3)[3]), (coef(modchl3)[2]+coef(modchl3)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modchl3)[1]+coef(modchl3)[4]), (coef(modchl3)[2]+coef(modchl3)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, -1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')



## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
plot(log(data$PP.biomass)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38.5,41), ylim=c(0,6), xlab='inv Temperature (C)', ylab='PP biomass ln(ug C / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,2,4,6), pos=38.5, lwd=2, cex.lab=1.5)
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
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

plot((data$PO4)~data$invT, col = data$trophic.level)
plot((data$PO4)~data$Tank, col = data$trophic.level)
plot((data$NO3.NO2)~data$invT, col = data$trophic.level)
plot((data$NO3.NO2)~data$Tank, col = data$trophic.level)

## Does zooplankton carbon vary with temperature?  
## figures 
hist(data$zoo.ug.carbon.liter)
plot(log(week8$zoo.ug.carbon.liter)~week8$Tank, pch = 19, col = week8$trophic.level)
plot(log(week8$zoo.ug.carbon.liter)~week8$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-1,6), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ZP biomass ln(ug C / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-1, lwd=2, cex.lab=1.5)
axis(2, at=c(-1,0,1,2,3,4,5,6), pos=38.5, lwd=2, cex.lab=1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$zoo.ug.carbon.liter)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$zoo.ug.carbon.liter)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
#week8 <- week8[-which(week8$Tank=='2'),]
modzpc0<-lm(log(week8$zoo.ug.carbon.liter)~1)
modzpc1<-lm(log(week8$zoo.ug.carbon.liter)~1+week8$invT)
modzpc2<-lm(log(week8$zoo.ug.carbon.liter)~1+week8$invT+week8$trophic.level)
modzpc3<-lm(log(week8$zoo.ug.carbon.liter)~1+week8$invT*week8$trophic.level)
anova(modzpc0, modzpc1)
anova(modzpc1, modzpc2)
anova(modzpc2, modzpc3)
AIC(modzpc0, modzpc1, modzpc2, modzpc3)
summary(modzpc2)
confint(modzpc2)

## add lines to plot
abline(coef(modzpc2)[1], coef(modzpc2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modzpc2)[1]+coef(modzpc2)[3]), (coef(modzpc2)[2]), lty = 3, lwd = 3, col = 'blue')
legend(38.5, 3, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')


## Does zooplankton density vary with temperature?  
## figures 
hist(log(week8$total.zoo.abundance.liter))
plot(log(week8$total.zoo.abundance.liter)~week8$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-2,4), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ZP density ln(ind / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-2, lwd=2, cex.lab=1.5)
axis(2, at=c(-2,-1,0,1,2,3,4), pos=38.5, lwd=2, cex.lab=1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$total.zoo.abundance.liter)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$total.zoo.abundance.liter)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modzp0<-lm(log(week8[(week8$trophic.level!='P'),]$total.zoo.abundance.liter)~1)
modzp1<-lm(log(week8[(week8$trophic.level!='P'),]$total.zoo.abundance.liter)~1+week8[(week8$trophic.level!='P'),]$invT)
modzp2<-lm(log(week8[(week8$trophic.level!='P'),]$total.zoo.abundance.liter)~1+week8[(week8$trophic.level!='P'),]$invT+week8[(week8$trophic.level!='P'),]$trophic.level)
modzp3<-lm(log(week8[(week8[(week8$trophic.level!='P'),]$trophic.level!='P'),]$total.zoo.abundance.liter)~1+week8[(week8$trophic.level!='P'),]$invT*week8[(week8$trophic.level!='P'),]$trophic.level)
anova(modzp0, modzp1)
anova(modzp1, modzp2)
anova(modzp1, modzp3)


## add lines to plot
abline(coef(modzp1)[1], coef(modzp1)[2], lty = 2, lwd = 3, col = 1)
#abline((coef(modzp3)[1]+coef(modzp3)[3]), (coef(modzp3)[2]+coef(modzp3)[4]), lty = 3, lwd = 3, col = 'blue')
legend(38.5, 4, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')



## Does daphnia:copepod  vary with temperature?  
## figures 
hist(log(week8$Daphnia.Copepod.Ratio))
plot(log(week8$Daphnia.Copepod.Ratio)~week8$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-3,3), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='Daphnia: Copepod') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2,3), pos=38.5, lwd=2, cex.lab=1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$Daphnia.Copepod.Ratio)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$Daphnia.Copepod.Ratio)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modzp0<-lm(log(week8$total.zoo.abundance.liter)~1)
modzp1<-lm(log(week8$total.zoo.abundance.liter)~1+week8$invT)
modzp2<-lm(log(week8$total.zoo.abundance.liter)~1+week8$invT+week8$trophic.level)
modzp3<-lm(log(week8$total.zoo.abundance.liter)~1+week8$invT*week8$trophic.level)
anova(modzp0, modzp1)
anova(modzp1, modzp2)
anova(modzp1, modzp3)



## Does zooplankton body size vary with temperature?  
## figures 
hist(log(week8$community.size))
plot(log(week8$community.size)~week8$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-1,0), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ZP mean size') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-1, lwd=2, cex.lab=1.5)
axis(2, at=c(-1,-0.5, 0), pos=38.5, lwd=2, cex.lab=1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$community.size)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$community.size)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modzp0<-lm(log(week8$community.size)~1)
modzp1<-lm(log(week8$community.size)~1+week8$invT)
modzp2<-lm(log(week8$community.size)~1+week8$invT+week8$trophic.level)
modzp3<-lm(log(week8$community.size)~1+week8$invT*week8$trophic.level)
anova(modzp0, modzp1)
anova(modzp1, modzp2)
anova(modzp2, modzp3)

abline(coef(modzp2)[1], coef(modzp2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modzp2)[1]+coef(modzp2)[3]), (coef(modzp2)[2]), lty = 3, lwd = 3, col = 'blue')
legend(38.5, 8, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')

## Does adult zooplankton density vary with temperature?  
## figures 
hist(week8$total.adults)
plot((week8$total.adults)~week8$invT, cex=1.5, pch='',  axes=FALSE,  ylim=c(0,8), xlim=c(38.5,41),  xlab='inv Temperature (C)', ylab='ZP adult density (ind / L)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,2,4,6,8), pos=38.5, lwd=2, cex.lab=1.5)
points((week8[(week8$trophic.level=='PZ'),]$total.adults)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points((week8[(week8$trophic.level=='PZN'),]$total.adults)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modzp0<-lm((week8$total.adults)~1)
modzp1<-lm((week8$total.adults)~1+week8$invT)
modzp2<-lm((week8$total.adults)~1+week8$invT+week8$trophic.level)
modzp3<-lm((week8$total.adults)~1+week8$invT*week8$trophic.level)
anova(modzp0, modzp1)
anova(modzp1, modzp2)
anova(modzp1, modzp3)


## add lines to plot
abline(coef(modzp3)[1], coef(modzp3)[2], lty = 2, lwd = 3, col = 'brown')
#abline((coef(modzp3)[1]+coef(modzp3)[3]), (coef(modzp3)[2]+coef(modzp3)[4]), lty = 3, lwd = 3, col = 'blue')
legend(38.5, 8, c('2 TL','3 TL'), pch = c(15, 17), col = c('brown', 'blue'), bty = 'n')



## Does total biomass vary with temperature and trophic structure?  
## figures 
hist(log(week8$totalC))
plot(log(week8$totalC)~week8$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38,41), ylim=c(0,8), xlab='inv Temperature (C)', ylab='Biomass ln(ug C/L)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,2,4,6,8), pos=38, lwd=2, cex.lab=1.5)
points(log(week8[(week8$trophic.level=='P'),]$totalC)~week8[(week8$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$totalC)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$totalC)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modTCm0<-lm(log(data$totalC)~1)
modTCm1<-lm(log(data$totalC)~1+data$invT)
modTCm2<-lm(log(data$totalC)~1+data$invT+data$trophic.level)
modTCm3<-lm(log(data$totalC~1+data$invT*data$trophic.level)
anova(modTCm0, modTCm1)
anova(modTCm1, modTCm2)
anova(modTCm1, modTCm3)
summary(modTCm1)
confint(modTCm1)

## add lines to plot
abline(coef(modTCm1)[1], coef(modTCm1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, 3, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))



## Does H:A vary with temperature and trophic structure?  
## figures 
week8$HA <- week8$zoo.ug.carbon.liter / week8$Pcarbon
hist(log(week8$HA))
plot(log(week8$HA)~week8$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38,41), ylim=c(-4, 0), xlab='inv Temperature (C)', ylab='Biomass ln(ug C/L)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=-4, lwd=2, cex.lab=1.5)
axis(2, at=c(-4,-2,0), pos=38, lwd=2, cex.lab=1.5)
points(log(week8[(week8$trophic.level=='PZ'),]$HA)~week8[(week8$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(week8[(week8$trophic.level=='PZN'),]$HA)~week8[(week8$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

modHA0<-lm(log(week8[(week8$trophic.level!='P'),]$HA)~1)
modHA1<-lm(log(week8[(week8$trophic.level!='P'),]$HA)~1+week8[(week8$trophic.level!='P'),]$invT)
modHA2<-lm(log(week8[(week8$trophic.level!='P'),]$HA)~1+week8[(week8$trophic.level!='P'),]$invT+week8[(week8$trophic.level!='P'),]$trophic.level)
modHA3<-lm(log(week8[(week8$trophic.level!='P'),]$HA)~1+week8[(week8$trophic.level!='P'),]$invT*week8[(week8$trophic.level!='P'),]$trophic.level)
            anova(modHA0, modHA1)
            anova(modHA1, modHA2)
            anova(modHA1, modHA3)
            summary(modTCm1)
            confint(modTCm1)

## plotting ZP biomass x PP biomass
plot(data$PP.biomass ~ data$total.carbon)
abline(0, 1)

plot(data$PP.biomass ~ data$zooplankton.carbon.per.L, pch = 19, col = c(data$average.temp))
abline(0, 1)
