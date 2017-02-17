### Garzke et al temperature experiment
### this is the current working code file
### MO made a new file on Aug29 from June file when I decided to take week out as a fixed effect, and instead model autocorrelation. 

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

### set working directory and load data
data <- read.csv("./data/temporal_dataFEB12.csv")
temps <- read.csv("./data/dailytemps.csv")
dim(data)
head(data)
tail(data)
data <- data[-(241:255),]
data$Tank <- as.character(data$Tank)
#View(data)

# process temperature data ------------------------------------------------

### extract temps from datalogger data, and only use these temps.
temps2 <- melt(temps, id = c("Hours", "Date", "Week"))
names(temps2) <- c('time', 'date','week','Tank', 'temp')  
temps3 <- tidyr::separate(temps2, Tank, c("X", "Tank"), sep = 1)
temps3 <- temps3[,-4]
temps3 <- tidyr::separate(temps3, date, c("Day", "Month", "Year"), sep = "/")

## average temp over each week for each tank
temps.wk <- 
  temps3 %>% 
  group_by(Tank, week) %>%
  summarize(mean(temp))

names(temps.wk) <- c("Tank", "week", "temp.wk")
data.t <- join(data, temps.wk, by = c("week", "Tank")) # add weekly temps to data file

## average temp over each tank over all weeks. 
temps.Tmn <- 
  temps3 %>% 
  group_by(Tank) %>%
  summarize(mean(temp, na.rm = TRUE))

names(temps.Tmn) <- c("Tank", "temp.Tmn")
data.t <- join(data.t, temps.Tmn, by = c("Tank")) 

## temperature at the time of measurement of oxygen for oxygen exchange corrections. These temps are only needed for the abiotic corrections on oxygen flux. 
## get dates and times so that I can pick the time I want from the temps3 file
data.t$date <- (dmy(data.t$date))
data.t <- tidyr::separate(data.t, dawn1time, c("d1Hour", "d1Min", "d1Sec"), sep = ":")
data.t <- tidyr::separate(data.t, dusktime, c("dkHour", "dkMin", "dkSec"), sep = ":")
data.t <- tidyr::separate(data.t, dawn2time, c("d2Hour", "d2Min", "d2Sec"), sep = ":")
data.t <- tidyr::separate(data.t, date, c("Year", "Month", "Date"), sep = "-")
data.t <- select(data.t, -contains("Min"))
data.t <- select(data.t, -contains("Sec"))
data.t <- select(data.t, -contains("calc"))


temps3$Year <- as.numeric(temps3$Year)
temps3$Month <- as.numeric(temps3$Month)
temps3$Day <- as.numeric(temps3$Day)
temps3$time <- as.numeric(temps3$time)
data.t$Month <- as.numeric(data.t$Month)
data.t$Date <- as.numeric(data.t$Date)

data.t2 <- data.t %>%  
  mutate(d1Hour = as.numeric(data.t$d1Hour)) %>% #time.d1
  mutate(d2Hour = as.numeric(data.t$d2Hour)) %>% #time.d2
  mutate(dkHour = as.numeric(data.t$dkHour)) #dkHour

data.t2 <- data.t2 %>% 
  unite(date_complete, Year, Month, Date, sep = "-") %>%
  mutate(date_formatted = ymd(date_complete)) 

## add date column to temps3
temps3$Year <- rep(2012, length(temps3[,1]))
temps4 <- temps3 %>% 
  unite(date_complete, Year, Month, Day, sep = "-") %>%
  mutate(date_formatted = ymd(date_complete)) %>% 
  filter(!is.na(date_formatted))
head(temps4)

temps4 <- temps4 %>% 
  mutate(T4hrs = rollmean(temp, 4, align = "right", fill = "NA"))

## join temps4 and data by the date, time and tank for each oxygen sampling time (hour)
data.t3 <- left_join(data.t2, temps4, by = c("date_formatted", "week", "Tank", "d1Hour" = "time")) #, suffix = c(".x", ".d1")
data.t3 <- dplyr::rename(data.t3, temp.d1 = T4hrs)

data.t3 <- left_join(data.t3, temps4, by = c("date_formatted", "week", "Tank", "dkHour" = "time")) #, suffix = c(".x", ".dk")
data.t3 <- dplyr::rename(data.t3, temp.dk = T4hrs)

data.t3 <- data.t3 %>%
  mutate(d2.date = date_formatted + 1)

data.t4 <- left_join(data.t3, temps4, by = c("d2.date" = "date_formatted", "Tank", "d2Hour" = "time")) #, suffix = c(".x", ".d2")
data.t4 <- dplyr::rename(data.t4, temp.d2 = T4hrs)


### Define temperature as inverse temperature
data <- data.t4
data <- dplyr::rename(data, week = week.x)
k <- 8.617342*10^-5  # eV/K
data$invTi <-  1/((data$temp.wk + 273)*k) # average temp of the tank each week
data$invTT <-  1/((data$temp.Tmn + 273)*k) # average temp of the tank over all weeks
data$invTi <- as.numeric(as.character(data$invTi))
data$invTT <- as.numeric(as.character(data$invTT))

## some plots to visualize for the two temperature terms:
plot((data$invTi - data$invTT) ~ data$week)
plot((data$invTi - data$invTT) ~ data$invTT)
plot((data$invTi) ~ data$invTT)
plot((data$invTi - data$invTT) ~ data$Tank)

### estimate biomass from chlorophyll concentration
data$PP.biomass <- (data$chla*55) #chla (ug/L)* 55 C in PP / 1 chla = ugPPC/L
data$ZP.carbon1 <- ifelse(data$trophic.level=='P',  0, data$zoo.ug.carbon.liter) # for adding
data$total.carbon <- data$PP.biomass + data$ZP.carbon1 #I'm assuming 0 for ZP in P treatments here.

### Estimating NPP and ER from the raw data: 

## correct of water-air oxygen flux
### adjusting for effects of temperature on physical exchange
C.star <- function(T) exp(7.7117 - 1.31403*log(T+45.93)) - 0.035
# the -0.035 is an approximate adjustment for elevation
C.star(T) # yields the oxygen concentration expected at a given temperature (T in C).


#calculate NPP and ER (hourly), in terms of umol O2 / l / hr, following yvon durochers 2010.
# O2 has molar mass of 32g/mol. so 1 umol = 32 ug. so take ug/32
data$NPP2 <- (((data$dusk - data$dawn1) - (C.star(data$temp.dk) - C.star(data$temp.d1)))*1000)/(32)  # oxygen produced umol / L /day, net all respiration. raw data is mg/L. Subtract o2 water/atm flux due to change in temperature 
data$ER2 <- -(24/data$hours2)*(((data$dawn2 - data$dusk) - (C.star(data$temp.d2) - C.star(data$temp.dk)))*1000)/(32)  # amount of oxygen consumed per day via respiaration. negative to get the change in oxygen umol / L /day; oxygen used in the dark and daylight. MeanER can be greater than meanNPP, because NPP reflects ER already.
data$GPP <- data$NPP2+(data$ER2/24)*data$hours1 # daily oxygen production (NPP2) + estimated daytime community respiration (daily R / 24 * hours daylight)
data$NEM <- data$ER2/data$GPP  # following Yvon Durochers 2010. NEM > 1 means the system is respiring more than it's fixing per day. This does not need to be logged.
data$NPP.mass <- data$NPP2 / (data$PP.biomass)  # NPP on ummol 02/L/day/ugCPP
data$ER.mass <- data$ER2/(data$total.carbon) # ER on ummol 02/L/day/ugTPP

data <- data[data$week >= '4',]
### data prep complete




# analysis begins ---------------------------------------------------------

### ANALYSIS
## trying the method suggested by Van de pol and Wright 2009
# center by within-tank temperature

data1 <- data
data1 <- data[(data$NPP2 >= 0.01),]
## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invTi - invTT), with Tank as a random intercept effect.  
modNPP0 <- lme(log(NPP2) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modNPP1 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")
modNPP2 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)  
modNPP4 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP5 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP6 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPP7 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modNPP0, modNPP2, modNPP4,modNPP1, modNPP5, modNPP6, modNPP7)

## Best model: create fitted values to use later for plotting
modNPP5r <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit)  
mod.coefs <- augment(modNPP5r, effect = "random")
## come back to get confints, might need qpCR

## Patrick's figure idea: dplyr
NPPmod = mod.coefs %>% 
  gather(key = "EF", value = "Rate", NPP2, ER2, PP.biomass) %>%
  group_by(Tank, trophic.level, EF) %>% 
  do(fitTank = lm(log(Rate) ~ invTi, data = .)) %>%
  tidy(., fitTank) %>%
  filter(term == 'invTi') %>%
  select(trophic.level, estimate) %>%
  mutate(level = "Tank")

xtanks <- data.frame(trophic.level = c("P", "PZ", "PZN"), estimate = c(-1.31564, -0.8241136, -0.4654028, fixef(modER2r)[3], fixef(modER2r)[3], fixef(modER2r)[3], 1.681596, 4.070138, 2.135455), Tank = "all", level = "treatment", EF = rep(c("NPP2", "ER2", "PP.biomass"), each = 3))

bind_rows(NPPmod, xtanks) -> data3

plot1 <- ggplot(data3, aes(x = trophic.level, y = estimate, color = trophic.level, size = level, pch = level)) +
  geom_point() +
  facet_grid(~EF) +
  theme_bw()

ggsave("plot1.png", device = "png", width = 5, height = 3)

#scales = "free" #in facet

plot(NPPmod.coef$trophic.level, NPPmod.coef$estimate)

## best averaged model: ## next challenge: figure out how to get coefs from the averaged model, given random effects. i'm not sure this is possible. 
## don't think it's possible. So, plot the fixed effect coefs for averaged models (need to pencil this out; leave off within group lines?)
## approach on 1/15/17 is to use best model
m.avg <- model.avg(modNPP5, modNPP2, modNPP1)
summary(m.avg)
pred.data <- as.data.frame(predict(m.avg, se.fit = FALSE)) # this is promising to get predicted values for plotting cheat below.
data.NPP <- cbind(data1, pred.data)

### exploring ways to extract coefficients from an averaging object
coeffs(m.avg)
coefTable(m.avg)

### SOME BASIC PLOTS
#data1 <- data[(data$NPP2 >= 0.5),] # three negative values and one very small value now, not sure what to do about them.
hist(data[(data$NPP2 >= 0.05),]$NPP2)
hist(log(data$NPP2))
hist((data1$NPP2))

plot(log(data1$NPP2)~data1$Tank, pch = 19, col = data1$trophic.level)
plot(log(data1$NPP2)~data1$week, pch = 19, col = data1$trophic.level)
plot(log(data1$NPP2)~I(data1$invTi-data1$invTT), pch = 19, col = data1$Tank, ylim = c(0,6))
plot(log(data1$NPP2)~I(data1$invTi-data1$invTT), pch = 19, col = data1$Tank, ylim = c(0,6))
abline(3.5435395, -0.6159088, lwd = 2, col = 1)
plot(log(data1$NPP2)~data1$invTi, pch = data1$Tank, col = data1$Tank)

#### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

NPP.plot <- ggplot(data = data1, aes(x = invTi, y = log(NPP2), ymin = -2)) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(aes(group = Tank, shape = trophic.level), alpha = 1/2, size = 2) + #color = trophic.level
  xlab("Temperature 1/kTi") +
  ylab("ln(NPPi)")


# NPP.PLOT ----------------------------------------------------------------

## NPP.PLOT 1: Individual regression lines fitted within and among groups
NPP.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = trophic.level), color = "gray3") 
 # geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black')
ggsave("NPPplot.png", device = "png", width = 5, height = 3)

## PLOT 2: Raw data and fitted lines from the model. Added the predictions of the model to the original dataset, then fit lines to those using linear regressions
## ok now that there are trophic group effects, we need to make different slopes: 
NPP.funcP <- function(x) { (fixef(modNPP5r)[1] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT)) + fixef(modNPP5r)[5]*x}
yvalsP <- NPP.funcP(mod.coefs$invTT)

NPP.funcPZ <- function(x) { (fixef(modNPP5r)[1] + fixef(modNPP5r)[3] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT) - fixef(modNPP5r)[8]*mean(mod.coefs$invTT)) + (fixef(modNPP5r)[5] + fixef(modNPP5r)[8])*x}
yvalsPZ <- NPP.funcPZ(mod.coefs$invTT)

NPP.funcPZN <- function(x) { (fixef(modNPP5r)[1] + fixef(modNPP5r)[4] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT) - fixef(modNPP5r)[9]*mean(mod.coefs$invTT)) + (fixef(modNPP5r)[5] + fixef(modNPP5r)[9])*x}
yvalsPZN <- NPP.funcPZN(mod.coefs$invTT)


NPP.plot +
  #geom_smooth(data = mod.coefs, aes(x = invTi, y = yvalsP, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23, size = .8) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 2, size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 4, size = 1.5) +
  geom_text(label = "Y.P = -1.32x + 54.96", x = 38.9, y = -1) +
  geom_text(label = "Y.PZ = -0.82x + 36.12", x = 38.9, y = -1.5) +
  geom_text(label = "Y.PZN = -0.46x + 22.03", x = 38.9, y = -2)
  #geom_point(data = mod.coefs, aes(x = invTi, y = pred.data))
  
ggsave("NPPplot.png", device = "png", width = 5, height = 3)


################################################
## Does mass-specific NPP vary with temperature?  
################################################
### NPP in umol/L/day per ugC
data1 <- data[(data$NPP2 >= 0.05),]

## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invTi - invTT), with Tank as a random intercept effect.  
modNPPm0 <- lme(log(NPP.mass) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modNPPm1 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")
modNPPm2 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)  
modNPPm4 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm5 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm6 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modNPPm7 <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modNPPm0, modNPPm2, modNPPm4, modNPPm1, modNPPm5, modNPPm6, modNPPm7)

## Best model: create fitted values to use later for plotting
modNPPm5r <- lme(log(NPP.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
intervals(modNPPm5r, which = "fixed") 
mod.coefs <- augment(modNPPm5r, effect = "random") # puts fitted values back in the original dataset.

## figures  
### NPP in umol/L/day per ugC
hist(data1$NPP.mass)
hist(log(data1$NPP.mass))
# hist(log(data$NPP.mass+.001)) # no reason to add a number, there are no -inf values.

plot(log(data$NPP.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$NPP.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$NPP.mass)~I(data1$invTT-mean(data1$invTT)), pch = 19, col = data1$trophic.level)
plot(log(data1$NPP.mass)~I(data1$invT-(data1$invTT)), pch = 19, col = data1$trophic.level)
abline(-2.62, -1.83, lwd = 2, col = 1)
abline((-2.62+1.53), (-1.83-1.40), lwd = 2, col = 2)
abline((-2.62+0.13), (-1.83-0.05), lwd = 2, col = 3)


### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

NPPm.plot <- ggplot(data = data1, aes(x = invTi, y = log(NPP.mass))) + 
  theme_bw() +
  geom_point(aes(shape = trophic.level), colour = 312) +
  xlab("Temperature 1/kTi") +
  ylab("ln(NPPmi)")

## PLOT 1: Individual regression lines fitted within and among groups
NPPm.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black')

## PLOT 2: Fit lines to fitted data from the model [close here - go get the coefficients as estimated for NPP, and write them out here with a z value. should plot up well.]
NPPm.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black')

## PLOT 3: Plot lines from model
z <- 0.5 #invTi - invTT for each tank; 

# PP coefs
NPPmPP.func <- function(x) { (fixef(modNPPm6r)[1] - fixef(modNPPm6r)[3]*mean(data1$invTT) - fixef(modNPPm6r)[6]*z*mean(mod.coefs$invTT)) + (fixef(modNPPm6r)[3] + fixef(modNPPm6r)[6]*z)*x } #x = invTT # slope = -2.47
# use this function to compute yvals for plotting.
yvalsPP <- NPPmPP.func(mod.coefs$invTT)

# ZP coefs
NPPmZP.func <- function(x) { (fixef(modNPPm6r)[1] - fixef(modNPPm6r)[3]*mean(data1$invTT) - fixef(modNPPm6r)[6]*z*mean(mod.coefs$invTT) + fixef(modNPPm6r)[4] - fixef(modNPPm6r)[7]*mean(mod.coefs$invTT)) + (fixef(modNPPm6r)[3] + fixef(modNPPm6r)[6]*z + fixef(modNPPm6r)[7])*x } #x = invTT #slope = -3.98
# use this function to compute yvals for plotting.
yvalsZP <- NPPmZP.func(mod.coefs$invTT)

# PZN coefs
NPPmPZN.func <- function(x) { (fixef(modNPPm6r)[1] - fixef(modNPPm6r)[3]*mean(data1$invTT) - fixef(modNPPm6r)[6]*z*mean(mod.coefs$invTT) + fixef(modNPPm6r)[5] - fixef(modNPPm6r)[8]*mean(mod.coefs$invTT)) + (fixef(modNPPm6r)[3] + fixef(modNPPm6r)[6]*z + fixef(modNPPm6r)[8])*x } #x = invTT, #slope = -2.38
yvalsPZN <- NPPmPZN.func(mod.coefs$invTT)

NPPm.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = Tank, colour = trophic.level), size = 0.5, method = "lm", se = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(values = c("#999999","#666666","#000000")) +
  #geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsPP, ymin = yvalsPP - 0.3, ymax = yvalsPP + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsPP), lwd = 2, color = "#999999") +
  #geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsZP, ymin = yvalsZP - 0.3, ymax = yvalsZP + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsZP), lwd = 2, color = "#666666") +
  #geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsPZN, ymin = yvalsPZN - 0.3, ymax = yvalsPZN + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsPZN), lwd = 2, color = "#000000") +

ggsave("NPPmplot.png", device = "png")


######################################
## Does ER vary with temperature?  
#####################################
hist(data$ER2)
data1 <- data[(data$ER2 >= 0),]
hist(log(data1$ER2))

#analysis
modER0<-lme(log(ER2) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modER1<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modER2<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER4<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER5 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)
modER6 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modER7 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modER0, modER1, modER2, modER4, modER5, modER6, modER7)

## Best model: create fitted values to use later for plotting 
modER2r<-lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
summary(modER2r)
intervals(modER2r, which = "fixed")
mod.coefsER <- augment(modER2r, effect = "random") 

m.avg <- model.avg(modER2, modER4)
summary(m.avg)

### SOME BASIC PLOTS
## figures 
plot(log(data$ER2)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER2)~data$week, pch = 19, col = data$trophic.level)
plot(log(data1$ER2)~I(data1$invTT-mean(data1$invTT)), pch = 19, col = data1$trophic.level, ylim = c(0,6))
abline(4.27, -0.43, lwd = 2, col = 1)
abline((4.27+0.35), (-0.43), lwd = 2, col = 2)
abline((3.88-0.04), (-0.43), lwd = 2, col = 3)

### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs
ER.plot <- ggplot(data = mod.coefs, aes(x = invTi, y = log(ER2))) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(aes(group = Tank, shape = trophic.level), alpha = 1/2, size = 2) +
  xlab("Temperature 1/kTi") +
  ylab("ln(ERi)")

ER.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23) +
  geom_smooth(method = "lm", se = FALSE, aes(group = trophic.level), color = "gray3")

ggsave("ERplot.png", device = "png", width = 5, height = 4)

ER.funcP <- function(x) { (fixef(modER2r)[1] - fixef(modER2r)[3]*mean(mod.coefs$invTT)) + (fixef(modER2r)[3])*x } # for trophic level 1
yvalsP <- ER.funcP(mod.coefs$invTT)

ER.funcPZ <- function(x) { (fixef(modER2r)[1] + fixef(modER2r)[4] - fixef(modER2r)[3]*mean(mod.coefs$invTT)) + (fixef(modER2r)[3])*x } # for trophic level 2
yvalsPZ <- ER.funcPZ(mod.coefs$invTT)

ER.funcPZN <- function(x) { (fixef(modER2r)[1] + fixef(modER2r)[5] - fixef(modER2r)[3]*mean(mod.coefs$invTT)) + (fixef(modER2r)[3])*x } # for trophic level 3
yvalsPZN <- ER.funcPZN(mod.coefs$invTT)

ER.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23, size = .8) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 2, size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 4, size = 1.5) +
  geom_text(label = "Y.P = -0.69x + 31.30", x = 38.87, y = 3.25) +
  geom_text(label = "Y.PZ = -0.69x + 31.79", x = 38.87, y = 3) +
  geom_text(label = "Y.PZN = -0.69x + 31.49", x = 38.87, y = 2.75)
#geom_point(data = mod.coefs, aes(x = invTi, y = pred.data))
  
ggsave("ERplot.png", device = "png", width = 5, height = 3)

  



####
##################################################
## Does mass specific ER vary with temperature? 
##################################################
## figures 
data1 <- data[(data$ER2 >= 0),]
hist(data1$ER.mass)
hist(log(data1$ER.mass))

# ER mass -----------------------------------------------------------------


plot(log(data$ER.mass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$ER.mass)~I(data$invTT-mean(data$invTT)), pch = 19, col = data$trophic.level)
abline(-1.89, -1.49, lwd = 2, col = 1)
abline((-1.89+1.54), (-1.49-0.74), lwd = 2, col = 2)
abline((-1.89+0.21), (-1.49+0.43), lwd = 2, col = 3)

#analysis
modERm0 <- lme(log(ER.mass) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modERm1 <- lme(log(ER.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML") 
modERm2 <- lme(log(ER.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modERm4<-lme(log(ER.mass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modERm5 <- lme(log(ER.mass) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)
modERm6 <- lme(log(ER.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)
modERm7 <- lme(log(ER.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modERm0, modERm1, modERm2, modERm4, modERm5, modERm6, modERm7)

## Best model: create fitted values to use later for plotting 
modERm6r <- lme(log(ER.mass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit) 
summary(modERm6r)
intervals(modERm6r, which = "fixed")
mod.coefs <- augment(modERm6r, effect = "random") 

### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

ERm.plot <- ggplot(data = data1, aes(x = invTi, y = log(ER.mass))) + 
  theme_bw() +
  geom_point(aes(group = Tank, color = trophic.level)) +
  xlab("Temperature 1/kTi") +
  ylab("ln(ERmi)")

## PLOT 1: Individual regression lines fitted within and among groups
ERm.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = 'black')

## PLOT 2: Fit lines to fitted data from the model [close here - go get the coefficients as estimated for NPP, and write them out here with a z value. should plot up well.]
ERm.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black')

## PLOT 3: Plot lines from model
# PP coefs
ERmPP.func <- function(x) { (fixef(modERm6r)[1] - fixef(modERm6r)[3]*mean(data1$invTT) - fixef(modERm6r)[6]*z*mean(mod.coefs$invTT)) + (fixef(modERm6r)[3] + fixef(modERm6r)[6]*z)*x } #x = invTT # slope = -2.47
# use this function to compute yvals for plotting.
yvalsPP <- ERmPP.func(mod.coefs$invTT)

# ZP coefs
ERmZP.func <- function(x) { (fixef(modERm6r)[1] - fixef(modERm6r)[3]*mean(data1$invTT) - fixef(modERm6r)[6]*z*mean(mod.coefs$invTT) + fixef(modERm6r)[4] - fixef(modERm6r)[7]*mean(mod.coefs$invTT)) + (fixef(modERm6r)[3] + fixef(modERm6r)[6]*z + fixef(modERm6r)[7])*x } #x = invTT #slope = -3.98
# use this function to compute yvals for plotting.
yvalsZP <- ERmZP.func(mod.coefs$invTT)

# PZN coefs
ERmPZN.func <- function(x) { (fixef(modERm6r)[1] - fixef(modERm6r)[3]*mean(data1$invTT) - fixef(modERm6r)[6]*z*mean(mod.coefs$invTT) + fixef(modERm6r)[5] - fixef(modERm6r)[8]*mean(mod.coefs$invTT)) + (fixef(modERm6r)[3] + fixef(modERm6r)[6]*z + fixef(modERm6r)[8])*x } #x = invTT, #slope = -2.38
yvalsPZN <- ERmPZN.func(mod.coefs$invTT)

z <- 0.5 #invTi - invTT for each tank; 

ERm.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsPP, ymin = yvalsPP - 0.3, ymax = yvalsPP + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsPP), lwd = 2) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsZP, ymin = yvalsZP - 0.3, ymax = yvalsZP + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsZP), lwd = 2, color = "seagreen") +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsPZN, ymin = yvalsPZN - 0.3, ymax = yvalsPZN + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsPZN), lwd = 2, color = "blue")
  
  ggsave("ERmplot.png", device = "png")


  
##### BIOMASS RESULTS #####

# PP biomass --------------------------------------------------------------
#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))

plot(log(data$PP.biomass)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$PP.biomass)~I(data$invTT-mean(data$invTT)), pch = 19, col = data$trophic.level)

abline(6.11, 0.87, col = 1, lwd = 2)
abline((6.11 -1.23), (.87 + 1.47), col = 2, lwd = 2)
abline((6.11 - 0.01), (0.87 + 0.25), col = 3, lwd = 2)

#analysis
modPP0 <- lme(log(PP.biomass) ~ 1, random = ~ 1 | Tank, data=data, na.action=na.omit, method="ML")  
modPP1 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, na.action=na.omit, method="ML") 
modPP2 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit) 
modPP4<-lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit) 
modPP5 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPP6 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPP7 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit) 

model.sel(modPP0, modPP1, modPP2, modPP4, modPP5, modPP6, modPP7) 

## Best model: create fitted values to use later for plotting  
modPP6r <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data, method="REML", na.action=na.omit) 
summary(modPP6r)
intervals(modPP6r, which = "fixed")
mod.coefs <- augment(modPP6r, effect = "random") 


### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

PP.plot <- ggplot(data = data, aes(x = invTi, y = log(PP.biomass))) + 
  theme_bw() +
  theme(legend.position = "FALSE") +
  geom_point(aes(group = Tank, shape = trophic.level), alpha = 1/2, size = 2) +
  xlab("Temperature 1/kTi") +
  ylab("ln(PPi)")

## PLOT 1: Individual regression lines fitted within and among groups
PP.plot +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = trophic.level), color = 'black')

ggsave("PPplot.png", device = "png", height = 3, width = 5)

## PLOT 2: Use fitted lines from the model. 
PP.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = .fitted, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = .fitted), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black')

## PLOT 3: Plot lines from model
# PP coefs
PP.PP.func <- function(x) { (fixef(modPP6r)[1] - fixef(modPP6r)[3]*mean(data$invTT) - fixef(modPP6r)[6]*z*mean(mod.coefs$invTT)) + (fixef(modPP6r)[3] + fixef(modPP6r)[6]*z)*x } #x = invTT # slope = -2.47
# use this function to compute yvals for plotting.
yvalsPP <- PP.PP.func(mod.coefs$invTT)

# ZP coefs
PP.ZP.func <- function(x) { (fixef(modPP6r)[1] - fixef(modPP6r)[3]*mean(data$invTT) - fixef(modPP6r)[6]*z*mean(mod.coefs$invTT) + fixef(modPP6r)[4] - fixef(modPP6r)[7]*mean(mod.coefs$invTT)) + (fixef(modPP6r)[3] + fixef(modPP6r)[6]*z + fixef(modPP6r)[7])*x } #x = invTT #slope = -3.98
# use this function to compute yvals for plotting.
yvalsZP <- PP.ZP.func(mod.coefs$invTT)

# PZN coefs
PP.PZN.func <- function(x) { (fixef(modPP6r)[1] - fixef(modPP6r)[3]*mean(data$invTT) - fixef(modPP6r)[6]*z*mean(mod.coefs$invTT) + fixef(modPP6r)[5] - fixef(modPP6r)[8]*mean(mod.coefs$invTT)) + (fixef(modPP6r)[3] + fixef(modPP6r)[6]*z + fixef(modPP6r)[8])*x } #x = invTT, #slope = -2.38
yvalsPZN <- PP.PZN.func(mod.coefs$invTT)

z <- 0.5 #invTi - invTT for each tank; 

## this is good, it has an attempt at confidence intervals. but they're not quite the right intervals. so the next plot below has none.
PP.plot +
  geom_smooth(data = mod.coefs, aes(x = invTi, y = (.fitted), group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsPP, ymin = yvalsPP - 0.3, ymax = yvalsPP + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsPP), lwd = 2) +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsZP, ymin = yvalsZP - 0.3, ymax = yvalsZP + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsZP), lwd = 2, color = "seagreen") +
  geom_ribbon(aes(x = (mod.coefs$invTT), y = yvalsPZN, ymin = yvalsPZN - 0.3, ymax = yvalsPZN + 0.3), fill = "grey70", alpha = 0.6) +
  geom_line(aes(x = (mod.coefs$invTT), y = yvalsPZN), lwd = 2, color = "blue")

ggsave("PPplot.png", device = "png")


PP.plot +
  #geom_smooth(data = mod.coefs, aes(x = invTi, y = yvalsP, group = Tank, color = trophic.level), method = "lm", se = FALSE, inherit.aes = FALSE) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Tank, color = trophic.level), alpha = 0.23, size = .8) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsZP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 2, size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 4, size = 1.5) +
  geom_text(label = "Y.P = 1.68x - 59.75", x = 39.9, y = 3) +
  geom_text(label = "Y.PZ = 4.07x - 155.09", x = 39.9, y = 2.5) +
  geom_text(label = "Y.PZN = 2.13x -78.10", x = 39.9, y = 2)
#geom_point(data = mod.coefs, aes(x = invTi, y = pred.data))

ggsave("PPplot.png", device = "png", width = 5, height = 3)



#### community size plots

CS.plot <- ggplot(data = data, aes(x = invTi, y = log(size+1), colour = factor(trophic.level))) + 
  theme_bw() +
  geom_point()









## Does zooplankton carbon vary with temperature?  ### leaving this for now. 
## figures 
dataz <- data[(which(data$trophic.level != 'P')),] # because we don't have observations for the P treatments (we're just assuming they are 0) I don't htink we should analyze those tanks here.
hist(dataz$zoo.ug.carbon.liter)
hist(log(dataz$zoo.ug.carbon.liter))
hist(log(dataz$zoo.ug.carbon.liter+1))
plot(log(dataz$zoo.ug.carbon.liter)~dataz$Tank, pch = 19, col = data$trophic.level)

plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$Tank, pch = 19, col = dataz$trophic.level)
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$week, pch = 19, col = dataz$trophic.level)
plot(log(dataz[dataz$week >= '4',]$zoo.ug.carbon.liter+1)~dataz[dataz$week >= '4',]$invT, pch = 19, col = dataz$trophic.level)

## analysis
## determine need for random effects in the full model: 
modZa<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modZb<-lme(log(zoo.ug.carbon.liter+1)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 
anova(modZa, modZb)

modTCc<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~1|Tank, data=data, method="ML", na.action=na.omit)
anova(modTCc, modTCb)

modTC0<-lme(log(total.carbon) ~ 1, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modTC1<-lme(log(total.carbon)~1+I(invT-mean(invT)), random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modTC2<-lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))
modTC4<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random=~1|Tank, data=data, method="ML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

model.sel(modTC0, modTC1, modTC2, modTC4)

# for model fitting: 
modTC2r <- lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random=~1|Tank, data=data, method="REML", na.action=na.omit, correlation = corAR1(form = ~week|Tank))

summary(modTC2r)




### TOTAL CARBON #### 
## Does total biomass vary with temperature and trophic structure?  
## have to analyze this for weeks 4 and later, because we have no ZP data from weeks 2 and 3. 
## figures 
hist(data$total.carbon)
hist(log(data$total.carbon))

plot(log(data$total.carbon)~data$Tank, pch = 19, col = data$trophic.level)
plot(log(data$total.carbon)~data$week, pch = 19, col = data$trophic.level)
plot(log(data$total.carbon)~I(data$invT-mean(data$invT)), pch = 19, col = data$trophic.level)

## analysis
## determine need for random effects in the full model: 
modTCa<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTCb<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~1|week, data=data, method="ML", na.action=na.omit) 
anova(modTCa, modTCb)

modTC0<-lme(log(total.carbon) ~ 1, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTC1<-lme(log(total.carbon)~1+I(invT-mean(invT)), random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTC2<-lme(log(total.carbon)~1+I(invT-mean(invT))+trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  
modTC4<-lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="ML", na.action=na.omit)  

model.sel(modTC0, modTC1, modTC2, modTC4)
anova(modTC4, modTC2)
anova(modTC1, modTC2)
anova(modTC1, modTC0)

# for model fitting: 
modTC4r <- lme(log(total.carbon)~1+I(invT-mean(invT))*trophic.level, random = ~I(invT-mean(invT))|week, data=data, method="REML", na.action=na.omit)  

summary(modTC4r)


