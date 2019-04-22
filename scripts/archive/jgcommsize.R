### Jessica's community size code
### Jan 31 2017

### load libraries
library(MuMIn)
library(nlme)
library(plyr)
library(tidyverse) 
library(broom)
library(reshape2)
library(lubridate)
library(hms)
library(zoo)

### set working directory and load 
CommSize <- read.csv("./data/CommunitySizes.csv")
iew(CommSize)
dim(CommSize)
head(CommSize)
tail(CommSize)
#data <- data[-(241:255),]
CommSize$Tank <- as.character(CommSize$Tank)
#View(data)

# process temperature data ------------------------------------------------

### Define temperature as inverse temperature
CSdata <- CommSize
CSdata$invTi <- as.numeric(as.character(CSdata$invTi))
CSdata$invTT <- as.numeric(as.character(CSdata$invTT))

### data prep complete

# analysis begins ---------------------------------------------------------
## trying the method suggested by Van de pol and Wright 2009
# center by within-tank temperature
hist(CSdata$size)
hist(log(CSdata$size))
hist(log(CSdata$size+1))

## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invTi - invTT), with Tank as a random intercept effect.  
modCS0 <- lme(log(size) ~ 1, random = ~ 1 | Tank, data=CSdata, na.action=na.omit, method="ML")  
modCS1 <- lme(log(size) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=CSdata, na.action=na.omit, method="ML")
modCS2 <- lme(log(size) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=CSdata, method="ML", na.action=na.omit)  
modCS4 <- lme(log(size) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=CSdata, method="ML", na.action=na.omit) 
modCS5 <- lme(log(size) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=CSdata, method="ML", na.action=na.omit) 
modCS6 <- lme(log(size) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=CSdata, method="ML", na.action=na.omit) 
modCS7 <- lme(log(size) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=CSdata, method="ML", na.action=na.omit) 

model.sel(modCS0, modCS2, modCS4,modCS1, modCS5, modCS6, modCS7)

anova(modCS2, modCS6)
anova(modCS6, modCS4)
anova(modCS4, modCS5)
anova(modCS5, modCS1)
anova(modCS1, modCS7)
anova(modCS7, modCS0)

## Best model: create fitted values to use later for plotting
modCS2r <- lme(log(size+1) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=CSdata, method="REML", na.action=na.omit)
summary(modCS2r)
intervals(modCS2r, which = "fixed")
mod.coefs <- augment(modCS2r, effect = "random")

## best averaged model: ## next challenge: figure out how to get coefs from the averaged model, given random effects. i'm not sure this is possible. 
## don't think it's possible. So, plot the fixed effect coefs for averaged models (need to pencil this out; leave off within group lines?)
m.avg <- model.avg(modCS2, modCS6, modCS4, modCS5)
summary(m.avg)
confint(m.avg, level=.95)
pred.data <- as.data.frame(predict(m.avg, se.fit = FALSE)) # this is promising to get predicted values for plotting cheat below.
data.CS <- cbind(data, pred.data)

### exploring ways to extract coefficients from an averaging object
coeffs(m.avg)
coefTable(m.avg)

### i changed these functions, they were not right. Because the data here are just for the two trophic levels (PZ and PZN), we can't just copy the expression from other code used on data with 3 trophic levels. the numbers for the coefficients, and even their formulas, didn't align because the 'base' trophic level here is the PZ level, not the P. 
CS.funcPZ <- function(x) { (fixef(modCS2r)[1] - fixef(modCS2r)[3]*mean(mod.coefs$invTT)) + (fixef(modCS2r)[3])*x } # for trophic level 2
yvalsPZ <- CS.funcPZ(mod.coefs$invTT)

CS.funcPZN <- function(x) { (fixef(modCS2r)[1] + fixef(modCS2r)[4] - fixef(modCS2r)[3]*mean(mod.coefs$invTT)) + (fixef(modCS2r)[3])*x } # for trophic level 3
yvalsPZN <- CS.funcPZN(mod.coefs$invTT)

######## open symbols greyscale within-tank black among-tank

CS.plot <- ggplot(data = CommSize, aes(x = invTi, y = log(size+1), colour = trophic.level)) + 
  theme_bw(20) +
  theme(legend.position = "none") +
  geom_point(aes(group = Tank, shape = trophic.level), alpha = 1/2, size = 2) + #color = trophic.level
  scale_shape_manual(values = c(1, 2)) +
  scale_colour_manual(limits=c("PZ", "PZN"), breaks=c("PZ", "PZN"), values= c("grey50", "grey20")) +
  xlab("Temperature 1/kTi") +
  ylab("Average community size(mm)")
CS.plot

A <-  CS.plot + geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23, size = .8)
A
A1.1 <- A + geom_smooth(data=mod.coefs, aes(x = invTT, y = yvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 2, size = 1.5)
A1.1
A2 <- A1.1 + geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 4, size = 1.5)
A2

## i just rewrote the plot code this way, in the way i'm used to, but yours is fine.
CS.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23, size = .8) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 2, size = 1.5)

ggsave("CSplot.png", device = "png", width = 5, height = 3)
