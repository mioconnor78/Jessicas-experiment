### Mary's analysis of Jessica's experimet
### Feb 9, 2015
### building on WSN talk code for manuscript with complete and revised datafile
### figs and analyses for all data

## temp was missing for week 9, tank 11, so I just used the temp for week 8.

## use data file and libraries loaded in week 8.

## create datafile of weeks post bloom
data.pb <- data[which(data$week >= 4),]

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
modNPP0<-lm(log(dat$NPPd)~1)
modNPP1<-lm(log(dat$NPPd)~1+dat$invT)
modNPP2<-lm(log(dat$NPPd)~1+dat$invT+dat$week)
modNPP3<-lm(log(dat$NPPd)~1+dat$invT*dat$week)
modNPP4<-lm(log(dat$NPPd)~1+dat$invT*dat$week+dat$trophic.level)
modNPP5<-lm(log(dat$NPPd)~1+dat$invT*dat$week+dat$invT*dat$trophic.level)
modNPP6<-lm(log(dat$NPPd)~1+dat$invT*dat$week*dat$trophic.level)
modNPP7<-lm(log(dat$NPPd)~1+dat$invT*dat$trophic.level*dat$week)
anova(modNPP0, modNPP1)
anova(modNPP1, modNPP2)
anova(modNPP2, modNPP3)
anova(modNPP3, modNPP4)
anova(modNPP3, modNPP5)
anova(modNPP4, modNPP6)
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




### NPP simpler plot #####
#############
dat <- data.pb
hist(log(dat$NPPd))
plot(log(dat$NPPd)~dat$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-4,3), xlim=c(37.5,41), xlab='inv(Temperature) 1/eV', ylab='NPPd ln(mg O/L/d)') 
axis(1, at=c(37.5, 38,38.5,39, 39.5, 40,40.5, 41), pos=-4, lwd=2, cex.lab=1.5)
axis(2, at=c(-4,-3,-2,-1,0,1,2,3), pos=37.5, lwd=2, cex.lab=1.5)
abline(0, 0, lwd = 3, col = 1, lty = 2)
points(log(dat[(dat$trophic.level=='P'),]$NPPd)~dat[(dat$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZ'),]$NPPd)~dat[(dat$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZN'),]$NPPd)~dat[(dat$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')



### ERd #####
#############
dat <- data
hist(log(dat$ERd))
plot(log(dat$ERd)~dat$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-3,3), xlim=c(38,41), xlab='inv(Temperature) 1/eV', ylab='ERd ln(mg O/L/d)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=-3, lwd=2, cex.lab=1.5)
axis(2, at=c(-3,-2,-1,0,1,2,3), pos=38, lwd=2, cex.lab=1.5)
abline(0, 0, lwd = 3, col = 1, lty = 2)
points(log(dat[(dat$trophic.level=='P'),]$ERd)~dat[(dat$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZ'),]$ERd)~dat[(dat$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZN'),]$ERd)~dat[(dat$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modERd0<-lme(log(ERd)~1, random = ~ 1|Tank, data = dat, method = "ML")
modERd1<-lme(log(ERd)~1+invT, random = ~ 1|Tank, data = dat, method = "ML")
modERd2<-lme(log(ERd)~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML")
modERd3<-lme(log(ERd)~1+invT+week+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modERd4<-lme(log(ERd)~1+invT*week+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modERd5<-lme(log(ERd)~1+invT*trophic.level+week, random = ~ 1|Tank, data = dat, method = "ML")
modERd6<-lme(log(ERd)~1+invT*trophic.level*week , random = ~ 1|Tank, data = dat, method = "ML")
anova(modERd0, modERd1)
anova(modERd1, modERd2)
anova(modERd2, modERd3)
anova(modERd3, modERd4)
anova(modERd3, modERd5)
anova(modERd5, modERd6)
summary(modERd5)
fixef(modERd5)
coef(modERd5)

abline(fixef(modERd5)[1], fixef(modERd5)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((fixef(modERd5)[1]+fixef(modERd5)[3]), (fixef(modERd5)[2]+fixef(modERd5)[6]), lty = 1, lwd = 3, col = 'brown')
abline((fixef(modERd5)[1]+fixef(modERd5)[4]), (fixef(modERd5)[2]+fixef(modERd5)[7]), lty = 1, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')



### NEM #####
#############
dat <- data.pb[(dat$trophic.level!='PZN'),]
hist(log(dat$NEM))
plot(log(dat$NEM)~dat$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-5,3), xlim=c(38,41), xlab='inv(Temperature) 1/eV', ylab='NEM ln(mg O/L/d)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=-5, lwd=2, cex.lab=1.5)
axis(2, at=c(-5,-4,-3,-2,-1,0,1,2,3), pos=38, lwd=2, cex.lab=1.5)
abline(0, 0, lwd = 3, col = 1, lty = 2)
points(log(dat[(dat$trophic.level=='P'),]$NEM)~dat[(dat$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZ'),]$NEM)~dat[(dat$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZN'),]$NEM)~dat[(dat$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modNEM0<-lme(log(NEM)~1, random = ~ 1|Tank, data = dat, method = "ML")
modNEM1<-lme(log(NEM)~1+invT, random = ~ 1|Tank, data = dat, method = "ML")
modNEM1.1<-lme(log(NEM)~1+trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modNEM2.1<-lme(log(NEM)~1+trophic.level+week, random = ~ 1|Tank, data = dat, method = "ML")
modNEM2<-lme(log(NEM)~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML")
modNEM3<-lme(log(NEM)~1+invT+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
#modNEM4<-lme(log(NEM)~1+invT*week, random = ~ 1|Tank, data = dat, method = "ML")
modNEM5<-lme(log(NEM)~1+invT*trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modNEM5.1 <- lme(log(NEM)~1+invT*week + trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modNEM6<-lme(log(NEM)~1+invT+trophic.level*week , random = ~ 1|Tank, data = dat, method = "ML")
modNEM7<-lme(log(NEM)~1+invT*trophic.level*week , random = ~ 1|Tank, data = dat, method = "ML")
anova(modNEM0, modNEM1)
anova(modNEM1, modNEM2)
anova(modNEM1, modNEM3)
anova(modNEM3, modNEM5)
anova(modNEM3, modNEM6)
anova(modNEM6, modNEM7)
summary(modNEM6)
fixef(modNEM6)

abline(fixef(modNEM6)[1], fixef(modNEM6)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((fixef(modNEM6)[1]+fixef(modNEM6)[3]), fixef(modNEM6)[2], lty = 1, lwd = 3, col = 'brown')
abline((fixef(modNEM6)[1]+fixef(modNEM6)[4]), fixef(modNEM6)[2], lty = 1, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


abline(fixef(modNEM3)[1], fixef(modNEM3)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((fixef(modNEM3)[1]+fixef(modNEM3)[3]), fixef(modNEM3)[2], lty = 1, lwd = 3, col = 'brown')
abline((fixef(modNEM3)[1]+fixef(modNEM3)[4]), fixef(modNEM3)[2], lty = 1, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

abline(fixef(modNEM5)[1], fixef(modNEM5)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((fixef(modNEM5)[1]+fixef(modNEM5)[3]), (fixef(modNEM5)[2]+fixef(modNEM5)[4]), lty = 1, lwd = 3, col = 'brown')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

### GPP #####
#############
dat <- data
hist(log(dat$GPP))
plot(log(dat$GPP)~dat$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(-1,3), xlim=c(38,41), xlab='inv(Temperature) 1/eV', ylab='GPP ln(mg O/L/d)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=-5, lwd=2, cex.lab=1.5)
axis(2, at=c(-1,0,1,2,3), pos=38, lwd=2, cex.lab=1.5)
abline(0, 0, lwd = 3, col = 1, lty = 2)
points(log(dat[(dat$trophic.level=='P'),]$GPP)~dat[(dat$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZ'),]$GPP)~dat[(dat$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZN'),]$GPP)~dat[(dat$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modGPP0<-lme(log(GPP)~1, random = ~ 1|Tank, data = dat, method = "ML")
modGPP1<-lme(log(GPP)~1+invT, random = ~ 1|Tank, data = dat, method = "ML")
modGPP2<-lme(log(GPP)~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML")
modGPP3<-lme(log(GPP)~1+invT+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
#modNEM4<-lme(log(NEM)~1+invT*week, random = ~ 1|Tank, data = dat, method = "ML")
modGPP5<-lme(log(GPP)~1+invT*trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modNEM6<-lme(log(NEM)~1+invT*trophic.level*week , random = ~ 1|Tank, data = dat, method = "ML")
anova(modGPP0, modGPP1)
anova(modGPP1, modGPP2)
anova(modGPP2, modGPP3)
anova(modGPP2, modGPP5)
anova(modNEM3, modNEM6)
summary(modNEM2)
fixef(modNEM3)

abline(fixef(modNEM3)[1], fixef(modNEM3)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((fixef(modNEM3)[1]+fixef(modNEM3)[3]), fixef(modNEM3)[2], lty = 1, lwd = 3, col = 'brown')
abline((fixef(modNEM3)[1]+fixef(modNEM3)[4]), fixef(modNEM3)[2], lty = 1, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


### Metabolic Balance (MB) #####
#############
dat <- data
hist(log(dat$MB))
hist(sqrt(dat$MB))
hist((dat$MB))
plot(log(dat$MB)~dat$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(0,3), xlim=c(38,41), xlab='inv(Temperature) 1/eV', ylab='MB ln(mg O/L/d)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(-1,0,1,2,3), pos=38, lwd=2, cex.lab=1.5)
abline(0, 0, lwd = 3, col = 1, lty = 2)
points(log(dat[(dat$trophic.level=='P'),]$MB)~dat[(dat$trophic.level=='P'),]$invT, pch=1, col = 'seagreen', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZ'),]$MB)~dat[(dat$trophic.level=='PZ'),]$invT, pch=22, col = 'brown', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZN'),]$MB)~dat[(dat$trophic.level=='PZN'),]$invT, pch=2, col = 'blue')

## analysis
modMB0<-lme(log(MB)~1, random = ~ 1|Tank, data = dat, method = "ML")
modMB1<-lme(log(MB)~1+invT, random = ~ 1|Tank, data = dat, method = "ML")
#modMB2<-lme(log(MB)~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML")
modMB3<-lme(log(MB)~1+invT+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
#modNEM4<-lme(log(NEM)~1+invT*week, random = ~ 1|Tank, data = dat, method = "ML")
modMB5<-lme(log(MB)~1+invT*trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modMB6<-lme(log(MB)~1+invT*trophic.level*week , random = ~ 1|Tank, data = dat, method = "ML")
anova(modMB0, modMB1)
anova(modMB1, modMB3)
anova(modMB3, modMB5)
anova(modMB3, modMB6)
anova(modNEM3, modNEM6)
summary(modMB3)
fixef(modNEM3)

abline(fixef(modMB3)[1], fixef(modMB3)[2], lty = 1, lwd = 3, col = 'seagreen')
abline((fixef(modMB3)[1]+fixef(modMB3)[3]), fixef(modMB3)[2], lty = 1, lwd = 3, col = 'brown')
abline((fixef(modMB3)[1]+fixef(modMB3)[4]), fixef(modMB3)[2], lty = 1, lwd = 3, col = 'blue')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')


## Does total PP biomass vary with temperature and FCL?   
## figures 
dat <- data.pb
hist(log(dat$PP.biomass))
plot(log(dat$PP.biomass)~dat$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38,41), ylim=c(0,8), xlab='inv Temperature (C)', ylab='PP biomass ln(ug C / L)') 
axis(1, at=c(38, 38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,2,4,6,8), pos=38, lwd=2, cex.lab=1.5)
points(log(dat[(dat$trophic.level=='P'),]$PP.biomass)~dat[(dat$trophic.level=='P'),]$invT, pch=19, col = 'seagreen', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZ'),]$PP.biomass)~dat[(dat$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(dat[(dat$trophic.level=='PZN'),]$PP.biomass)~dat[(dat$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modPb0<-lme(log(PP.biomass)~1, random = ~ 1|Tank, data = dat, method = "ML")
modPb1<-lme(log(PP.biomass)~1+invT, random = ~ 1|Tank, data = dat, method = "ML")
modPb2<-lme(log(PP.biomass)~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML")
modPb3<-lme(log(PP.biomass)~1+invT+week+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modPb4<-lme(log(PP.biomass)~1+invT*week, random = ~ 1|Tank, data = dat, method = "ML")
modPb5<-lme(log(PP.biomass)~1+invT+trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modPb5<-lme(log(PP.biomass)~1+invT*week , random = ~ 1|Tank, data = dat, method = "ML")
anova(modPb0, modPb1)
anova(modPb1, modPb2)
anova(modPb2, modPb3)
anova(modPb2, modPb4)
summary(modPb3)
confint(modPb1)

## add lines to plot
abline(coef(modPb1)[1], coef(modPb1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

plot((data$PO4)~data$invT, col = data$trophic.level)
plot((data$PO4)~data$Tank, col = data$trophic.level)
plot((data$NO3.NO2)~data$invT, col = data$trophic.level)
plot((data$NO3.NO2)~data$Tank, col = data$trophic.level)


## Does PP Diversity?  
## figures on invT
data$PPdiv <- data$Phyto.Shannon.Wiener.Diversity
hist(data$PPdiv)
plot((data$PPdiv)~data$Tank, pch = 19, col = data$trophic.level)
plot((data$PPdiv)~data$invT, cex=1.5, pch='',  axes=FALSE, ylim=c(0,1.5), xlim=c(37.5,41), xlab='inv(Temperature) 1/eV', ylab='PP Shannon Div') 
axis(1, at=c(37.5, 38.0, 38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,.5, 1, 1.5), pos=37.5, lwd=2, cex.lab=1.5)
abline(12, -0.32, lwd = 3, col = 2)
points(data[(data$trophic.level=='P'),]$PPdiv~data[(data$trophic.level=='P'),]$invT, pch=19, col = 'green', cex = 1)
points(data[(data$trophic.level=='PZ'),]$PPdiv~data[(data$trophic.level=='PZ'),]$invT, pch=19, col = 'brown', cex = 1)
points(data[(data$trophic.level=='PZN'),]$PPdiv~data[(data$trophic.level=='PZN'),]$invT, pch=19, col = 'blue', cex = 1)


## analysis
dat <- data
modPD0<-lme(PPdiv~1, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD1<-lme(PPdiv~1+invT, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD2<-lme(PPdiv~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD3<-lme(PPdiv~1+invT+week+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD4<-lme(PPdiv~1+trophic.level, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD5<-lme(PPdiv~1+week, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD6<-lme(PPdiv~1+trophic.level*week, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD7<-lme(PPdiv~1+trophic.level*invT, random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
modPD5<-lme(PPdiv~1+invT*week , random = ~ 1|Tank, data = dat, method = "ML", na.action = na.omit)
anova(modPD0, modPD1)
anova(modPD0, modPD2)
anova(modPD0, modPD4)
anova(modPD0, modPD7)



## Does total ZP biomass vary with temperature and FCL?   
## figures 
dat <- data

plot(log(data$zoo.ug.carbon.liter+0.001) ~ data$trophic.level, ylab = 'log ZP carbon + 0.01', xlab = 'trophic level')
summary(aov(log(data[(data$trophic.level!='P'),]$zoo.ug.carbon.liter+0.001) ~ data[(data$trophic.level!='P'),]$trophic.level))


hist((dat$zoo.ug.carbon.liter))
plot((dat$zoo.ug.carbon.liter)~dat$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38,41), ylim=c(0,200), xlab='inv Temperature (C)', ylab='ZP biomass (ug C / L)') 
axis(1, at=c(38, 38.5,39, 39.5, 40,40.5, 41), pos=0, lwd=2, cex.lab=1.5)
axis(2,  pos=38, lwd=2, cex.lab=1.5) #at=c(0,2,4,6,8),
points((dat[(dat$trophic.level=='PZ'),]$zoo.ug.carbon.liter)~dat[(dat$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points((dat[(dat$trophic.level=='PZN'),]$zoo.ug.carbon.liter)~dat[(dat$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

## analysis
modZb0<-lme(log(zoo.ug.carbon.liter)~1, random = ~ 1|Tank, data = dat, method = "ML")
modZb1<-lme(log(zoo.ug.carbon.liter)~1+invT, random = ~ 1|Tank, data = dat, method = "ML")
modZb2<-lme(log(zoo.ug.carbon.liter)~1+invT+week, random = ~ 1|Tank, data = dat, method = "ML")
modZb3<-lme(log(zoo.ug.carbon.liter)~1+invT+week+ trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modZb4<-lme(log(zoo.ug.carbon.liter)~1+invT*week, random = ~ 1|Tank, data = dat, method = "ML")
modZb5<-lme(log(zoo.ug.carbon.liter)~1+invT+trophic.level, random = ~ 1|Tank, data = dat, method = "ML")
modZb5<-lme(log(zoo.ug.carbon.liter)~1+invT*week , random = ~ 1|Tank, data = dat, method = "ML")
anova(modPb0, modPb1)
anova(modPb1, modPb2)
anova(modPb2, modPb3)
anova(modPb2, modPb4)
summary(modPb3)
confint(modPb1)

## add lines to plot
abline(coef(modPb1)[1], coef(modPb1)[2], lty = 1, lwd = 3, col = 'black')
legend(40.5, 2, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'), bty = 'n')

plot((data$PO4)~data$invT, col = data$trophic.level)
plot((data$PO4)~data$Tank, col = data$trophic.level)
plot((data$NO3.NO2)~data$invT, col = data$trophic.level)
plot((data$NO3.NO2)~data$Tank, col = data$trophic.level)



## Does H:A vary with temperature and trophic structure?  
## figures 

plot(data$HA ~ data$trophic.level)
summary(aov(data$HA ~ data$trophic.level))


data$HA <- data$zoo.ug.carbon.liter / data$Pcarbon
hist(log(data$HA+0.001))
plot(log(data$HA+0.001)~data$invT, cex=1.5, pch='',  axes=FALSE, xlim=c(38,41), ylim=c(-8, 1), xlab='inv Temperature (C)', ylab='H/A biomass ratio ln(ug C/L)') 
axis(1, at=c(38,38.5,39, 39.5, 40,40.5, 41), pos=-8, lwd=2, cex.lab=1.5)
axis(2, at=c(-8,-4,-2,0,2), pos=38, lwd=2, cex.lab=1.5)
points(log(data[(data$trophic.level=='PZ'),]$HA+0.001)~data[(data$trophic.level=='PZ'),]$invT, pch=15, col = 'brown', cex = 1.5)
points(log(data[(data$trophic.level=='PZN'),]$HA+0.001)~data[(data$trophic.level=='PZN'),]$invT, pch=17, col = 'blue')

modHA0<-lm(log(data[(data$trophic.level!='P'),]$HA+0.001)~1, na.action = na.omit)
modHA1<-lm(log(data[(data$trophic.level!='P'),]$HA+0.001)~1+data[(data$trophic.level!='P'),]$invT, na.action = na.omit)
modHA2<-lm(log(week8[(week8$trophic.level!='P'),]$HA+0.001)~1+week8[(week8$trophic.level!='P'),]$invT+week8[(week8$trophic.level!='P'),]$trophic.level)
modHA3<-lm(log(week8[(week8$trophic.level!='P'),]$HA+0.001)~1+week8[(week8$trophic.level!='P'),]$invT*week8[(week8$trophic.level!='P'),]$trophic.level)
anova(modHA0, modHA1)
anova(modHA1, modHA2)
anova(modHA1, modHA3)
summary(modTCm1)
confint(modTCm1)

## plotting ZP biomass x PP biomass
plot(log(data[(data$trophic.level!='PZN'),]$zoo.ug.carbon.liter+0.001) ~ log(data[(data$trophic.level!='PZN'),]$Pcarbon), pch = 19, ylab='ZP carbon ln(ug C / L)', xlab = 'Phyto carbon ln(ug C / L)', col = 'brown')
points(log(data[(data$trophic.level!='PZ'),]$zoo.ug.carbon.liter+0.001) ~ log(data[(data$trophic.level!='PZ'),]$Pcarbon),pch = 19, col = 'blue')


abline(0, 1)

plot(data$PP.biomass ~ data$zooplankton.carbon.per.L, pch = 19, col = c(data$average.temp))
abline(0, 1)
