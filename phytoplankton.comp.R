setwd("/Users/maryo/Documents/temporary files/Jessicas experment")
pp <- read.csv("PPweek7.csv")
dim(pp)
k <- 8.62*10^-5
pp$invT <-  1/((pp$average.temp + 273)*k)
pp$invT <- as.numeric(as.character(pp$invT))
names(pp)
hist(pp$proportion.cyano, breaks=20, xlim=c(0,1), col='gray')
hist(pp$TOTAL, breaks=20, col='gray')

plot(pp$Diatom.cyano.ratio~pp$average.temp, cex=1.5, pch='',  axes=FALSE, ylim=c(0,0.5), xlim=c(14,28), xlab='inv Temperature (C)', ylab='Diatom/cyano density') 
axis(1, at=c(14,18,22,24,28), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,.25,.5), pos=14, lwd=2, cex.lab=1.5)
points(pp[(pp$trophic.level=='P'),]$Diatom.cyano.ratio~pp[(pp$trophic.level=='P'),]$average.temp, pch=19, col = 'seagreen', cex = 1.5)
points(pp[(pp$trophic.level=='PZ'),]$Diatom.cyano.ratio~pp[(pp$trophic.level=='PZ'),]$average.temp, pch=15, col = 'brown', cex = 1.5)
points(pp[(pp$trophic.level=='PZN'),]$Diatom.cyano.ratio~pp[(pp$trophic.level=='PZN'),]$average.temp, pch=17, col = 'blue')
legend(26, 0.5, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))

## does the proportion of cells that are cyanos change with temperature?
## C axes
plot(pp$proportion.cyano~pp$average.temp, cex=1.5, pch='',  axes=FALSE, ylim=c(0,1), xlim=c(14,28), xlab='Temperature (C)', ylab='Cyanos/total cells') 
axis(1, at=c(14,16,18,20,22,24,26,28), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,.25,.5,0.75,1), pos=14, lwd=2, cex.lab=1.5)
points(pp[(pp$trophic.level=='P'),]$proportion.cyano~pp[(pp$trophic.level=='P'),]$average.temp, pch=19, col = 'seagreen', cex = 1.5)
points(pp[(pp$trophic.level=='PZ'),]$proportion.cyano~pp[(pp$trophic.level=='PZ'),]$average.temp, pch=15, col = 'brown', cex = 1.5)
points(pp[(pp$trophic.level=='PZN'),]$proportion.cyano~pp[(pp$trophic.level=='PZN'),]$average.temp, pch=17, col = 'blue')
legend(15, 1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))

## analysis
hist(pp$proportion.cyano, breaks=20, xlim=c(0,1), col='gray')
modPc0<-lm((pp$proportion.cyano)~1)
modPc1<-lm((pp$proportion.cyano)~1+pp$average.temp)
modPc2<-lm((pp$proportion.cyano)~1+pp$average.temp+pp$trophic.level)
modPc3<-lm((pp$proportion.cyano)~1+pp$average.temp*pp$trophic.level)
anova(modPc0, modPc1)
anova(modPc1, modPc2)
anova(modPc1, modPc3)
AIC(modPPc0, modPPc1, modPPc2, modPPc3)

abline(coef(modPc1)[1], coef(modPc1)[2], lty = 1, lwd = 3, col = 'black')

## does the proportion of cells that are cyanos change with temperature?
## invT axes
plot(pp$proportion.cyano~pp$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(0,1), xlim=c(38.5,41), xlab='Temperature (C)', ylab='Cyanos/total cells') 
axis(1, at=c(38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,.25,.5,0.75,1), pos=38.5, lwd=2, cex.lab=1.5)
points(pp[(pp$trophic.level=='P'),]$proportion.cyano~pp[(pp$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(pp[(pp$trophic.level=='PZ'),]$proportion.cyano~pp[(pp$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(pp[(pp$trophic.level=='PZN'),]$proportion.cyano~pp[(pp$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')
legend(40.5, 1, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))

## analysis
hist(pp$proportion.cyano, breaks=20, xlim=c(0,1), col='gray')
modPc0<-lm((pp$proportion.cyano)~1)
modPc1<-lm((pp$proportion.cyano)~1+pp$invT)
modPc2<-lm((pp$proportion.cyano)~1+pp$invT+pp$trophic.level)
modPc3<-lm((pp$proportion.cyano)~1+pp$invT*pp$trophic.level)
anova(modPc0, modPc1)
anova(modPc1, modPc2)
anova(modPc1, modPc3)
AIC(modPPc0, modPPc1, modPPc2, modPPc3)

abline(coef(modPc1)[1], coef(modPc1)[2], lty = 1, lwd = 3, col = 'black')

plot(pp$TOTAL~pp$average.temp, cex=1.5, pch='',ylim=c(0,400),  axes=FALSE,  xlim=c(14,28), xlab='Temperature (C)', ylab='total cells') 
axis(1, at=c(14,16,18,20,22,24,26,28), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,100,200,300,400), pos=14, lwd=2, cex.lab=1.5)
points(pp[(pp$trophic.level=='P'),]$TOTAL~pp[(pp$trophic.level=='P'),]$average.temp, pch=19, col = 'seagreen', cex = 1.5)
points(pp[(pp$trophic.level=='PZ'),]$TOTAL~pp[(pp$trophic.level=='PZ'),]$average.temp, pch=15, col = 'brown', cex = 1.5)
points(pp[(pp$trophic.level=='PZN'),]$TOTAL~pp[(pp$trophic.level=='PZN'),]$average.temp, pch=17, col = 'blue')
legend(26, 400, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))

## analysis
modPPc0<-lm((pp$TOTAL)~1)
modPPc1<-lm((pp$TOTAL)~1+pp$average.temp)
modPPc2<-lm((pp$TOTAL)~1+pp$average.temp+pp$trophic.level)
modPPc3<-lm((pp$TOTAL)~1+pp$average.temp*pp$trophic.level)
anova(modPPc0, modPPc1)
anova(modPPc0, modPPc2)
anova(modPPc0, modPPc3)
AIC(modPPc0, modPPc1, modPPc2, modPPc3)

summary(modPPc3)
confint(modPPc3)
coef(modPPc3)

abline(coef(modPPc3)[1], coef(modPPc3)[2], lty = 2, lwd = 3, col = 'seagreen')
abline((coef(modzpc3)[1]+coef(modzpc3)[3]), (coef(modzpc3)[2]+coef(modzpc3)[4]), lty = 3, lwd = 3, col = 'blue')
