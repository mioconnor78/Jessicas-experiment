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
modNPP1<-lm(log(data$NPP)~data$invT)
modNPP2<-lm(log(data$NPP)~data$invT*data$trophic.level)
anova(modNPP1, modNPP2)
summary(modNPP2)
coef(modNPP2)

## add lines to plot
abline(coef(modNPP2)[1], coef(modNPP2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPP2)[1]+coef(modNPP2)[3]), (coef(modNPP2)[2]+coef(modNPP2)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPP2)[1]+coef(modNPP2)[4]), (coef(modNPP2)[2]+coef(modNPP2)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')



## Does mass-specific NPP vary with temperature?  
## figures
hist(data$NPP)
plot(log(data$NPP)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-3,2), xlim=c(38.5,41), xlab='inv Temperature (C)', ylab='NPP ln(mg/L/hr)') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2), pos=38.5, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='P'),]$NPP)~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(data[(data$trophic.level=='PZ'),]$NPP)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$NPP)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modNPP1<-lm(log(data$NPP)~data$invT)
modNPP2<-lm(log(data$NPP)~data$invT*data$trophic.level)
anova(modNPP1, modNPP2)
summary(modNPP2)
coef(modNPP2)

## add lines to plot
abline(coef(modNPP2)[1], coef(modNPP2)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((coef(modNPP2)[1]+coef(modNPP2)[3]), (coef(modNPP2)[2]+coef(modNPP2)[5]), lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPP2)[1]+coef(modNPP2)[4]), (coef(modNPP2)[2]+coef(modNPP2)[6]), lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')
