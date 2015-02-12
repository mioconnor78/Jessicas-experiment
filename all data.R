### Mary's analysis of Jessica's experimet
### Feb 9, 2015
### building on WSN talk code for manuscript with complete and revised datafile
### figs and analyses for all data

## temp was missing for week 9, tank 11, so I just used the temp for week 8.

## use data file and libraries loaded in week 8.

## Does NPP vary with temperature?  
## figures on invT
hist(data$NPPd)
plot(log(data$NPPd)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPPd)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-4,2), xlim=c(37.5,41), xlab='inv(Temperature) 1/eV', ylab='NPP ln(mg O/L/d)') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40,40.5, 41), pos=-4, lwd=2, cex.lab=1.5)
axis(2, at=c(-4, -3,-2,-1,0,1,2), pos=37.5, lwd=2, cex.lab=1.5)
abline(12, -0.32, lwd = 3, col = 2)
points(log(data[(data$trophic.level=='P'&data$week=='2'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='2'),]$invT, pch=19, col = 'green', cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='3'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='3'),]$invT, pch=19, col = 'seagreen', cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='4'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='4'),]$invT, pch=19, col = 'seagreen1', cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='5'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='5'),]$invT, pch=19, col = 493, cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='6'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='6'),]$invT, pch=19, col = 91, cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='7'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='7'),]$invT, pch=19, col = 393, cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='8'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='8'),]$invT, pch=19, col = 394, cex = 1)
points(log(data[(data$trophic.level=='P'&data$week=='9'),]$NPPd)~data[(data$trophic.level=='P'&data$week=='9'),]$invT, pch=19, col = 396, cex = 1)

points(log(data[(data$trophic.level=='PZ'&data$week=='2'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='2'),]$invT, pch=15, col = 'green', cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='3'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='3'),]$invT, pch=15, col = 'seagreen', cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='4'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='4'),]$invT, pch=15, col = 'seagreen1', cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='5'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='5'),]$invT, pch=15, col = 493, cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='6'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='6'),]$invT, pch=15, col = 91, cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='7'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='7'),]$invT, pch=15, col = 393, cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='8'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='8'),]$invT, pch=15, col = 394, cex = 1)
points(log(data[(data$trophic.level=='PZ'&data$week=='9'),]$NPPd)~data[(data$trophic.level=='PZ'&data$week=='9'),]$invT, pch=15, col = 396, cex = 1)

points(log(data[(data$trophic.level=='PZN'&data$week=='2'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='2'),]$invT, pch=17, col = 'green', cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='3'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='3'),]$invT, pch=17, col = 'seagreen', cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='4'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='4'),]$invT, pch=17, col = 'seagreen1', cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='5'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='5'),]$invT, pch=17, col = 493, cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='6'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='6'),]$invT, pch=17, col = 91, cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='7'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='7'),]$invT, pch=17, col = 393, cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='8'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='8'),]$invT, pch=17, col = 394, cex = 1)
points(log(data[(data$trophic.level=='PZN'&data$week=='9'),]$NPPd)~data[(data$trophic.level=='PZN'&data$week=='9'),]$invT, pch=17, col = 396, cex = 1)


points(log(data[(data$trophic.level=='PZ'),]$NPPd)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = rainbow(data$week), cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$NPPd)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = rainbow(data$week))

#find rich.colors()

## analysis
modNPP0<-lm(log(data$NPPd)~1)
modNPP1<-lm(log(data$NPPd)~1+data$invT)
modNPP2<-lm(log(data$NPPd)~1+data$invT+data$week)
modNPP3<-lm(log(data$NPPd)~1+data$invT*data$week)
modNPP4<-lm(log(data$NPPd)~1+data$invT*data$week+data$trophic.level)
modNPP5<-lm(log(data$NPPd)~1+data$invT*data$week+data$invT*data$trophic.level)
modNPP6<-lm(log(data$NPPd)~1+data$invT*data$week*data$trophic.level)
anova(modNPP0, modNPP1)
anova(modNPP1, modNPP2)
anova(modNPP2, modNPP3)
anova(modNPP3, modNPP4)
anova(modNPP4, modNPP5)
anova(modNPP5, modNPP6)
AIC(modNPP0, modNPP1, modNPP2, modNPP3,modNPP4, modNPP5, modNPP6)

summary(modNPP6)
coef(modNPP6)
confint(modNPP2)

## add lines to plot
abline((coef(modNPP6)[1]+coef(modNPP6)[3]), (coef(modNPP6)[2]+coef(modNPP6)[6]), lty = 1, lwd = 3, col = 'green')
abline((coef(modNPP6)[1]+2*coef(modNPP6)[3]+coef(modNPP6)[4]+2*coef(modNPP6)[9]), (coef(modNPP6)[2]+coef(modNPP6)[7]+2*coef(modNPP6)[11]), lty = 2, lwd = 3, col = 'green')
abline((coef(modNPP6)[1]+coef(modNPP6)[5]), (coef(modNPP6)[2]+coef(modNPP6)[8]), lty = 3, lwd = 3, col = 'green')

abline((coef(modNPP2)[1]+coef(modNPP2)[3]), coef(modNPP2)[2], lty = 2, lwd = 3, col = 'brown')
abline((coef(modNPP2)[1]+coef(modNPP2)[4]), coef(modNPP2)[2], lty = 3, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))
