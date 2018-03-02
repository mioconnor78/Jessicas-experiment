### Garzke et al temperature experiment
### scripts for final manuscript
### Mary O'Connor 


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

## first set of code (to line 152) is processing

### set working directory and load data
data <- read.csv("./data/temporal_dataFEB12.csv")
temps <- read.csv("./data/dailytemps.csv")

data <- data[-(241:255),]
data$Tank <- as.character(data$Tank)

# process temperature data ------------------------------------------------

### extract temps from datalogger data, and only use these temps.
temps2 <- melt(temps, id = c("Hours", "Date", "Week"))
names(temps2) <- c('time', 'date','week','Tank', 'temp')  
temps3 <- tidyr::separate(temps2, Tank, c("X", "Tank"), sep = 1)
temps3 <- temps3[,-4]
temps3 <- tidyr::separate(temps3, date, c("Day", "Month", "Year"), sep = "/")

## Format date
temps3$Year <- rep(2012, length(temps3[,1]))
temps4 <- temps3 %>% 
  unite(date_complete, Year, Month, Day, sep = "-") %>%
  mutate(date_formatted = ymd(date_complete)) %>% 
  filter(!is.na(date_formatted))
View(temps4)

## there are too many readings for 7/3/2012, so remove this date:
temps4 <- temps4 %>% 
  filter(., date_formatted != "2012-07-03")

## estimate moving avg temp over four hours
temps4 <- temps4 %>% 
  mutate(T4hrs = rollmean(temp, 4, align = "right", fill = "NA"))

# calculate mean weekly temps
temps.wk <- temps4 %>%
  group_by(week, Tank) %>%
  summarise(., avg = mean(temp, na.rm = TRUE))

names(temps.wk) <- c("week","Tank", "temp.wk")
temps.wk <- temps.wk[!is.na(temps.wk$week),]
temps.wk$week <- as.integer(temps.wk$week)
temps.wk$Tank <- as.character(temps.wk$Tank)

# add weekly temps to data file
data.t <- left_join(data, temps.wk, by = c("week", "Tank")) 

## average temp over each tank over all weeks. 
temps.Tmn <- 
  temps4 %>% 
  group_by(Tank) %>%
  summarise(., avg = mean(temp, na.rm = TRUE), sd = sd(temp, na.rm = TRUE)) %>%
  arrange(as.numeric(Tank))

names(temps.Tmn) <- c("Tank", "temp.Tmn", "temp.Tsd")
data.t <- left_join(data.t, temps.Tmn, by = c("Tank")) 


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

data.t$Month <- as.numeric(data.t$Month)
data.t$Date <- as.numeric(data.t$Date)

data.t2 <- data.t %>%  
  mutate(d1Hour = as.numeric(data.t$d1Hour)) %>% #time.d1
  mutate(d2Hour = as.numeric(data.t$d2Hour)) %>% #time.d2
  mutate(dkHour = as.numeric(data.t$dkHour)) #dkHour

#format dates for merging with temp file
data.t2 <- data.t2 %>% 
  unite(date_complete, Year, Month, Date, sep = "-") %>%
  mutate(date_formatted = ymd(date_complete)) 

## join temps4 and data by the date, time and tank for each oxygen sampling time (hour)
data.t3 <- left_join(data.t2, temps4, by = c("date_formatted", "week", "Tank", "d1Hour" = "time")) #, suffix = c(".x", ".d1")
data.t3 <- dplyr::rename(data.t3, temp.d1 = T4hrs)

data.t3 <- left_join(data.t3, temps4, by = c("date_formatted", "week", "Tank", "dkHour" = "time")) #, suffix = c(".x", ".dk")
data.t3 <- dplyr::rename(data.t3, temp.dk = T4hrs)

data.t3 <- data.t3 %>%
  mutate(d2.date = date_formatted + 1)

data.t4 <- left_join(data.t3, temps4, by = c("d2.date" = "date_formatted", "Tank", "d2Hour" = "time")) #, suffix = c(".x", ".d2")
data.t4 <- dplyr::rename(data.t4, temp.d2 = T4hrs)



## create a column for tank rank within trophic trt
names(data.t4)
data.t5 <- as.tibble(data.t4) %>%
  group_by(Tank, trophic.level) %>%
  summarise(., mean.temp = mean(temp.wk)) %>%
  arrange(trophic.level, mean.temp) %>%
  select(-mean.temp) 
  
data.t5$Tankn <- rep(c(1:10), 3)

head(data.t5) 
data.t6 <- left_join(data.t4, data.t5, by = c("Tank", "trophic.level")) 
  
### Define temperature as inverse temperature
data <- data.t6
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

#calculate NPP and ER (hourly), in terms of umol O2 / l / hr, following yvon durochers 2010.
z <- 31.25 # for converting mg/L to umol O2

data$NPP2 <- (data$dusk - data$dawn1)/z  # oxygen produced umol / L /day, net all respiration. raw data is mg/L. 
data$ER2 <- -(24/data$hours2)*((data$dawn2 - data$dusk))/z  # amount of oxygen consumed per day via respiaration. negative to get the change in oxygen umol / L /day; oxygen used in the dark and daylight. MeanER can be greater than meanNPP, because NPP reflects ER already.
data$GPP <- data$NPP2+(data$ER2/24)*data$hours1 # daily oxygen production (NPP2) + estimated daytime community respiration (daily R / 24 * hours daylight)
data$NEM <- data$ER2/data$GPP  # following Yvon Durochers 2010. NEM > 1 means the system is respiring more than it's fixing per day. This does not need to be logged.
data$NPP.mass <- data$NPP2 / (data$PP.biomass)  # NPP on ummol 02/L/day/ugCPP
data$ER.mass <- data$ER2/(data$total.carbon) # ER on ummol 02/L/day/ugTPP

data <- data[data$week >= '4',]

### data prep complete


# analysis begins ---------------------------------------------------------

### ANALYSES
## following Van de pol and Wright 2009
## center by within-tank temperature
## model with mean tank temperature (invTT) and the weekly deviation from that long-term average (invTi - invTT), with Tank as a random intercept effect.  

# ANALYSES FOR FIGURE 2 ---------------------------------------------------
data1 <- data
data1 <- data[(data$NPP2 >= 0.001),] #remove 18 negative values 

# NPP candidate model set -------------------------------------------------
modNPPF <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP8 <- lme(log(NPP2) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP7 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP6 <- lme(log(NPP2) ~ 1 + trophic.level*I(invTi - invTT), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP5 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP4 <- lme(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP3 <- lme(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP2 <- lme(log(NPP2) ~ 1 + I(invTi - invTT), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP1 <- lme(log(NPP2) ~ 1 + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

modNPP0 <- lme(log(NPP2) ~ 1, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)

model.sel(modNPP0, modNPP1, modNPP2, modNPP3, modNPP4, modNPP5, modNPP6, modNPP7, modNPP8, modNPPF)

## Average best models: 
m.avgN <- model.avg(modNPP8, modNPPF)
confint(m.avgN)

## calculating confidence intervals for activation energies.
vcov(m.avgN) # diagonals of vcov are variances
vcov(m.avgN)[1,1] # intercept variance is:

df <- length(data1$NPP2)-(length(coefficients(m.avgN))+1)-1
t.stat <- qt(0.975, df = df) #calculates critical t-value for the threshold (first value) and df (= n - p - 1)

# slope for NPP.PP: 
Sl1 <- coefficients(m.avgN)[5]
Sl1.l <- confint(m.avgN)[5,1]
Sl1.u <- confint(m.avgN)[5,2]
Sl1.l2 <- Sl1 - t.stat * sqrt(vcov(m.avgN)[5,5])
Sl1.u2 <- Sl1 + t.stat * sqrt(vcov(m.avgN)[5,5])

# slope for NPP.ZP: 
Sl2 <- coefficients(m.avgN)[5] + coefficients(m.avgN)[8]
Sl2.l <- Sl2 - t.stat * sqrt(vcov(m.avgN)[5,5] + vcov(m.avgN)[8,8] + 2*vcov(m.avgN)[8,5])
Sl2.u <- Sl2 + t.stat * sqrt(vcov(m.avgN)[5,5] + vcov(m.avgN)[8,8] + 2*vcov(m.avgN)[8,5])
  
# slope for NPP.PZN: 
Sl3 <- coefficients(m.avgN)[5] + coefficients(m.avgN)[9]
Sl3.l <- Sl3 - t.stat * sqrt(vcov(m.avgN)[5,5] + vcov(m.avgN)[9,9] + 2*vcov(m.avgN)[9,5])
Sl3.u <- Sl3 + t.stat * sqrt(vcov(m.avgN)[5,5] + vcov(m.avgN)[9,9] + 2*vcov(m.avgN)[9,5])

slopesNPP <- (cbind(c(Sl1, Sl2, Sl3), c(Sl1.l2, Sl2.l, Sl3.l), c(Sl1.u2, Sl2.u, Sl3.u)))
rownames(slopesNPP) <- c("P", "PZ", "PZN")
colnames(slopesNPP) <- c("S", "l", "u")

# ints for NPP.PP:
I1 <- coefficients(m.avgN)[1] - coefficients(m.avgN)[5]*mean(data1$invTT)

# ints for NPP.ZP: m.avgN includes all int terms, but we leave out the invTi term for among group lines; ints here are at mean(invTT) 
I2 <- coefficients(m.avgN)[1] + coefficients(m.avgN)[2] - coefficients(m.avgN)[5]*mean(data1$invTT) - coefficients(m.avgN)[8]*mean(data1$invTT)

# ints for NPP.PZN: 
I3 <- coefficients(m.avgN)[1] + coefficients(m.avgN)[3] - coefficients(m.avgN)[5]*mean(data1$invTT) - coefficients(m.avgN)[9]*mean(data1$invTT)

#IntsNPP <- (cbind(c(I1, I2, I3), c(I1.l2, I2.l, I3.l), c(I1.u2, I2.u, I3.u)))
#rownames(IntsNPP) <- c("P", "PZ", "PZN")
#colnames(IntsNPP) <- c("I", "l", "u")

# FIGURE 2: 
### plotting within- and among-group regressions and model outputs
labels <- c(P = "Phytoplankton", PZ = "Phytoplankton + Grazers", PZN = "Phyto. + Grazers + Predators")

NPP.plot <- ggplot(data = data1, aes(x = -invTi, y = log(NPP2), ymin = -2, ymax = 6)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(.~trophic.level, labeller=labeller(trophic.level = labels)) + ## this sets it up as facets
  geom_point(aes(group = as.character(Tankn), color = as.character(Tankn), shape = as.factor(week), alpha = Tankn), size = 2) + 
  scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none") +
  scale_alpha("Tankn", guide = "none") +
  xlab("Temperature 1/kTi") +
  ylab("Oxygen Production (NPP) \n ln(umol O2 / L / hr)")

NPP.plot
ggsave("NPPplot.png", device = "png", width = 7, height = 3) # save for appendix

## PLOT 2A: Raw data and fitted lines from the averaged model. Added the predictions of the model to the original dataset (mod.coefsN), then fit lines to those using linear regressions
NPP.funcP <- function(x) { I1 + Sl1*x}
NvalsP <- NPP.funcP(data1[(data1$trophic.level=="P"),]$invTT)

NPP.funcPZ <- function(x) { I2 + Sl2*x }
NvalsPZ <- NPP.funcPZ(data1[(data1$trophic.level=="PZ"),]$invTT)

NPP.funcPZN <- function(x) { I3 + Sl3*x }
NvalsPZN <- NPP.funcPZN(data1[(data1$trophic.level=="PZN"),]$invTT)


# Figure 2A ---------------------------------------------------------------
# the within-group lines here are lms fitted to the actual data; I think these should be the modeled data too...
# the among-group lines are model fits based on the best model (so this could be model averaged coefficients too)
Fig2A <- 
  NPP.plot + 
  geom_smooth(data = subset(data1), method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTi, y = log(NPP2), group = Tank),  size = .8, color = alpha("steelblue", 0.5)) + 
  geom_smooth(data = subset(data1, trophic.level == "P"), aes(x = -invTT, y = NvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
 geom_smooth(data = subset(data1, trophic.level == "PZ"), aes(x = -invTT, y = NvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data1, trophic.level == "PZN"), aes(x = -invTT, y = NvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) 

ggsave("Fig2A-C.png", device = "png", width = 7, height = 3)




# ER candidate model set --------------------------------------------------

hist(data$ER2)
data2 <- data[(data$ER2 >= 0),] #no values removed
data2 <- data[!is.na(data$ER2),]
hist(log(data2$ER2))

#analysis
### might just ax the random int...or test for it: 
modERF <- lme(log(ER2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data2, method="ML", na.action=na.omit) 
modERa <- lm(log(ER2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data2, na.action=na.omit) #don't need random effect
## proceed without random effect
modERF <- lme(log(ER2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER8 <- lme(log(ER2) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER7 <- lme(log(ER2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER6 <- lme(log(ER2) ~ 1 + trophic.level*I(invTi - invTT), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER5 <- lme(log(ER2) ~ 1 + I(invTi - invTT) + trophic.level, random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER4 <- lme(log(ER2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER3 <- lme(log(ER2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER2 <- lme(log(ER2) ~ 1 + I(invTi - invTT), random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER1 <- lme(log(ER2) ~ 1 + trophic.level, random = ~ 1 | Tank, data=data2, na.action=na.omit)
modER0 <- lme(log(ER2) ~ 1, random = ~ 1 | Tank, data=data2, na.action=na.omit)

model.sel(modER0, modER1, modER2, modER3, modER4, modER5, modER6, modER7, modER8, modERF)

## calculating confidence intervals for activation energies.
modER <- modER7
vcov(modER) # diagonals of vcov are variances
vcov(modER)[1,1] # intercept variance is:

df <- length(data2$ER2)-(length(coef(modER))+1)-1
t.stat <- qt(0.975, df = df) #calculates critical t-value for the threshold (first value) and df (= n - p - 1)

# slope for ER.PP: 
SlER1 <- coefficients(modER)[5]
SlER1.l <- confint(modER)[5,1]
SlER1.u <- confint(modER)[5,2]
SlER1.l2 <- SlER1 - t.stat * sqrt(vcov(modER)[5,5])
SlER1.u2 <- SlER1 + t.stat * sqrt(vcov(modER)[5,5])

# slope for ER.ZP: 
SlER2 <- coefficients(modER)[5] + coefficients(modER)[6]
SlER2.l <- SlER2 - t.stat * sqrt(vcov(modER)[5,5] + vcov(modER)[6,6] + 2*vcov(modER)[6,5])
SlER2.u <- SlER2 + t.stat * sqrt(vcov(modER)[5,5] + vcov(modER)[6,6] + 2*vcov(modER)[6,5])

# slope for ER.PZN: 
SlER3 <- coefficients(modER)[5] + coefficients(modER)[7]
SlER3.l <- SlER3 - t.stat * sqrt(vcov(modER)[5,5] + vcov(modER)[7,7] + 2*vcov(modER)[7,5])
SlER3.u <- SlER3 + t.stat * sqrt(vcov(modER)[5,5] + vcov(modER)[7,7] + 2*vcov(modER)[7,5])

slopesER <- (cbind(c(SlER1, SlER2, SlER3), c(SlER1.l2, SlER2.l, SlER3.l), c(SlER1.u2, SlER2.u, SlER3.u)))
rownames(slopesER) <- c("P", "PZ", "PZN")
colnames(slopesER) <- c("S", "l", "u")

# ints for ER.PP: 
IER1 <- coefficients(modER)[1] - coefficients(modER)[5]*mean(data2$invTT)
IER1.l <- confint(modER)[1,1]
IER1.u <- confint(modER)[1,2]
IER1.l2 <- IER1 - t.stat * sqrt(vcov(modER)[1,1])

# ints for ER.ZP: modER includes all int terms, but we leave out the invTi term for among group lines; ints here are at mean(invTT) 
IER2 <- coefficients(modER)[1] + coefficients(modER)[3] - coefficients(modER)[5]*mean(data2$invTT) - coefficients(modER)[6]*mean(data2$invTT)
#IER2.l <- IER2 - t.stat * sqrt(vcov(modER)[1,1] + vcov(modER)[3,3] + 2*vcov(modER)[3,1] + vcov(modER)[5,5]*(mean(data2$invTT)^2) + vcov(modER)[6,6]*(mean(data2$invTT)^2)) 
#IER2.u <- IER2 + t.stat * sqrt(vcov(modER)[1,1] + vcov(modER)[3,3] + 2*vcov(modER)[3,1])

# ints for ER.PZN: 
IER3 <- coefficients(modER)[1] + coefficients(modER)[4] - coefficients(modER)[5]*mean(data2$invTT) - coefficients(modER)[7]*mean(data2$invTT)
#IER3.l <- IER3 - t.stat * sqrt(vcov(modER)[1,1] + vcov(modER)[4,4] + 2*vcov(modER)[4,1]) 
#IER3.u <- IER3 + t.stat * sqrt(vcov(modER)[1,1] + vcov(modER)[4,4] + 2*vcov(modER)[4,1])

slopesER <- (cbind(c(SlER1, SlER2, SlER3), c(SlER1.l, SlER2.l, SlER3.l), c(SlER1.u, SlER2.u, SlER3.u)))
rownames(slopesER) <- c("P", "PZ", "PZN")
colnames(slopesER) <- c("S", "l", "u")


# Fig 2 D-E: ER ------------------------------------------------------------

### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs
xlab <- expression(paste('Temperature (',~degree,'C)',sep=''))
ER.plot <- ggplot(data = data2, aes(x = -invTi, y = log(ER2), min = 0)) + 
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  facet_grid(.~trophic.level, labeller = labeller(trophic.level = labels)) + ## this sets it up as facets
  geom_point(aes(group = as.character(Tankn), color = as.character(Tankn), shape = as.factor(week), alpha = Tankn), size = 2) + 
  scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none") +
  scale_alpha("Tankn", guide = "none") +
  scale_shape(name = "Week", guide = guide_legend(ncol = 2)) +
  scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~((1/(k*-.))-273), name = xlab)) +
  xlab("") + #xlab("Temperature 1/kTi") +
  ylab("Oxygen Consumption (ER) \n ln(umol O2 / L / hr)")

ER.plot
ggsave("ERplot.png", device = "png", width = 7, height = 3) # save for appendix

## PLOT 2: Raw data and fitted lines from the model. Added the predictions of the model to the original dataset, then fit lines to those using linear regressions
ER.funcP <- function(x) {IER1 + SlER1*x} # for trophic level 1
x <- data2[(data2$trophic.level=="P"),]$invTT
RvalsP <- ER.funcP(x)

ER.funcPZ <- function(x) { IER2 + SlER2*x } # for trophic level 2
x <- data2[(data2$trophic.level=="PZ"),]$invTT
RvalsPZ <- ER.funcPZ(x)

ER.funcPZN <- function(x) { IER3 + SlER3*x } # for trophic level 3
x <- data2[(data2$trophic.level=="PZN"),]$invTT
RvalsPZN <- ER.funcPZN(x)


Fig2D <- 
ER.plot +
  geom_smooth(method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTi, y = log(ER2), group = Tank),  size = .8, color = alpha("steelblue", 0.5)) +
  #geom_smooth(method = "lm", se = FALSE, aes(group = Tank), color = "gray40", alpha = 0.23, size = .8) +
  geom_smooth(data = data2[(data2$trophic.level=="P"),], aes(x = -invTT, y = RvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data2, trophic.level == "PZ"), aes(x = -invTT, y = RvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 1, size = 1.5) +
  geom_smooth(data = subset(data2, trophic.level == "PZN"), aes(x = -invTT, y = RvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 1, size = 1.5) 

  #geom_text(label = "Y.P = -0.69x + 31.30", x = 38.87, y = 1.75) +
  #geom_text(label = "Y.PZ = -0.69x + 31.79", x = 38.87, y = 1) +
  #geom_text(label = "Y.PZN = -0.69x + 31.49", x = 38.87, y = .25)

ggsave("Fig2D-F.png", device = "png", width = 7, height = 3)


##### Chlorophyll A RESULTS #####

# Chla --------------------------------------------------------------
#### Chla ###
## Does total Chla vary with temperature and FCL?   
## figures 
hist(data$chla)
hist(log(data$chla))

#analysis
### might just ax the random int...or test for it: 
modPBF <- lme(log(chla) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit) 
modPBa <- lm(log(chla) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data, na.action=na.omit) 
## proceed without random effect...

modPB8 <- lme(log(chla) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB7 <- lme(log(chla) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB6 <- lme(log(chla) ~ 1 + trophic.level*I(invTi - invTT), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB5 <- lme(log(chla) ~ 1 + I(invTi - invTT) + trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB4 <- lme(log(chla) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB3 <- lme(log(chla) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB2 <- lme(log(chla) ~ 1 + I(invTi - invTT), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB1 <- lme(log(chla) ~ 1 + trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB0 <- lme(log(chla) ~ 1, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)

model.sel(modPB0, modPB1, modPB2, modPB3, modPB4, modPB5, modPB6, modPB7, modPB8, modPBF)

## calculating confidence intPBvals for activation enPBgies.
modPB <- modPB7
vcov(modPB) # diagonals of vcov are variances
vcov(modPB)[1,1] # intPBcept variance is:

df <- length(data$PP.biomass)-(length(coef(modPB))+1)-1
t.stat <- qt(0.975, df = df) #calculates critical t-value for the threshold (first value) and df (= n - p - 1)

# slope for PB.PP: 
SlPB1 <- fixef(modPB)[5]
SlPB1.l <- intervals(modPB)[1]
SlPB1.u <- intervals(modPB)[5,2]
SlPB1.l2 <- SlPB1 - t.stat * sqrt(vcov(modPB)[5,5])
SlPB1.u2 <- SlPB1 + t.stat * sqrt(vcov(modPB)[5,5])

# slope for PB.ZP: 
SlPB2 <- fixef(modPB)[5] + fixef(modPB)[6]
SlPB2.l <- SlPB2 - t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[6,6] + 2*vcov(modPB)[6,5])
SlPB2.u <- SlPB2 + t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[6,6] + 2*vcov(modPB)[6,5])

# slope for PB.PZN: 
SlPB3 <- fixef(modPB)[5] + fixef(modPB)[7]
SlPB3.l <- SlPB3 - t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[7,7] + 2*vcov(modPB)[7,5])
SlPB3.u <- SlPB3 + t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[7,7] + 2*vcov(modPB)[7,5])

slopesPB <- (cbind(c(SlPB1, SlPB2, SlPB3), c(SlPB1.l2, SlPB2.l, SlPB3.l), c(SlPB1.u2, SlPB2.u, SlPB3.u)))
rownames(slopesPB) <- c("P", "PZ", "PZN")
colnames(slopesPB) <- c("S", "l", "u")

# ints for PB.PP: 
IPB1 <- fixef(modPB)[1] - fixef(modPB)[5]*mean(data$invTT)
IPB1.l <- confint(modPB)[1,1]
IPB1.u <- confint(modPB)[1,2]
IPB1.l2 <- IPB1 - t.stat * sqrt(vcov(modPB)[1,1])

# ints for PB.ZP: modPB includes all int tPBms, but we leave out the invTi tPBm for among group lines; ints hPBe are at mean(invTT) 
IPB2 <- fixef(modPB)[1] + fixef(modPB)[3] - fixef(modPB)[5]*mean(data$invTT) - fixef(modPB)[6]*mean(data$invTT)
#IPB2.l <- IPB2 - t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[3,3] + 2*vcov(modPB)[3,1] + vcov(modPB)[5,5]*(mean(data$invTT)^2) + vcov(modPB)[6,6]*(mean(data$invTT)^2)) 
#IPB2.u <- IPB2 + t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[3,3] + 2*vcov(modPB)[3,1])

# ints for PB.PZN: 
IPB3 <- fixef(modPB)[1] + fixef(modPB)[4] - fixef(modPB)[5]*mean(data$invTT) - fixef(modPB)[7]*mean(data$invTT)
IPB3.l <- IPB3 - t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[4,4] + 2*vcov(modPB)[4,1]) 
IPB3.u <- IPB3 + t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[4,4] + 2*vcov(modPB)[4,1])

slopesPB <- (cbind(c(SlPB1, SlPB2, SlPB3), c(SlPB1.l, SlPB2.l, SlPB3.l), c(SlPB1.u, SlPB2.u, SlPB3.u)))
rownames(slopesPB) <- c("P", "PZ", "PZN")
colnames(slopesPB) <- c("S", "l", "u")

### WITHIN AND AMONG GROUP PLOTS
### plotting within- and among-group regressions and model outputs

PP.plot <- ggplot(data = data, aes(x = -invTi, y = log(chla), min = 0)) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        legend.position = c(0.94, 0.30), 
        legend.text=element_text(size=6), 
        legend.title = element_text(size = 7), 
        legend.key = element_rect(fill = NA), 
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  facet_grid(.~trophic.level) + 
  geom_point(aes(group = as.character(Tankn), 
                 color = as.character(Tankn), 
                 shape = as.factor(week), 
                 alpha = Tankn), 
             size = 2) + 
  scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none") +
  scale_alpha("Tankn", guide = "none") +
  scale_x_continuous("Temperature (1/kTi)", 
                     sec.axis = sec_axis(~((1/(k*-.))-273), 
                                         name = xlab)) +
  scale_shape(name = "Week", 
              guide = guide_legend(ncol = 2, 
                                   keywidth = .8, 
                                   keyheight = .8)) +
  xlab("Temperature 1/kTi") +
  ylab("Phytoplankton ug Chl a / L")

PP.plot
ggsave("PBplot.png", device = "png", width = 7, height = 3) 

PB.funcP <- function(x) {IPB1 + SlPB1*x} # for trophic level 1
x <- data[(data$trophic.level=="P"),]$invTT
PBvalsP <- PB.funcP(x)

PB.funcPZ <- function(x) { IPB2 + SlPB2*x } # for trophic level 2
x <- data2[(data$trophic.level=="PZ"),]$invTT
PBvalsPZ <- PB.funcPZ(x)

PB.funcPZN <- function(x) { IPB3 + SlPB3*x } # for trophic level 3
x <- data[(data$trophic.level=="PZN"),]$invTT
PBvalsPZN <- PB.funcPZN(x)

Fig2G <-
  PP.plot +
  geom_smooth(method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTi, y = log(chla), group = Tank),  size = .8, color = alpha("steelblue", 0.5)) +
  geom_smooth(data = subset(data, trophic.level == "P"), aes(x = -invTT, y = PBvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data, trophic.level == "PZ"), aes(x = -invTT, y = PBvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 1, size = 1.5) +
  geom_smooth(data = subset(data, trophic.level == "PZN"), aes(x = -invTT, y = PBvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 1, size = 1.5) 
  
ggsave("PPplot.png", device = "png")



## maybe put NPP back where was, and just use this:
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

Figure2 <- multiplot(Fig2A, Fig2D, Fig2G, cols = 1)
Fig2A
Fig2D
Fig2G
ggsave("Figure2.png", plot = Figure2, width = 7, height = 7)

png('Figure2.png', width = 7, height = 7, units = 'in', res = 300)
multiplot(Fig2A, Fig2D, Fig2G, cols = 1)
dev.off()


#### community size plots

CS.plot <- ggplot(data = data, aes(x = invTi, y = log(size+1), colour = factor(trophic.level))) + 
  theme_bw() +
  geom_point()




# PP biomass --------------------------------------------------------------
#### PP BIOMASS ###
## Does total PP biomass vary with temperature and FCL?   
## figures 
hist(data$PP.biomass)
hist(log(data$PP.biomass))

#analysis
### might just ax the random int...or test for it: 
modPBF <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit) 
modPBa <- lm(log(PP.biomass) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data, na.action=na.omit) 
## proceed without random effect...

modPB8 <- lme(log(PP.biomass) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB7 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB6 <- lme(log(PP.biomass) ~ 1 + trophic.level*I(invTi - invTT), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB5 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB4 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB3 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB2 <- lme(log(PP.biomass) ~ 1 + I(invTi - invTT), random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB1 <- lme(log(PP.biomass) ~ 1 + trophic.level, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)
modPB0 <- lme(log(PP.biomass) ~ 1, random = ~ 1 | Tank, data=data, method="ML", na.action=na.omit)

model.sel(modPB0, modPB1, modPB2, modPB3, modPB4, modPB5, modPB6, modPB7, modPB8, modPBF)

## calculating confidence intPBvals for activation enPBgies.
modPB <- modPB7
vcov(modPB) # diagonals of vcov are variances
vcov(modPB)[1,1] # intPBcept variance is:

df <- length(data$PP.biomass)-(length(coef(modPB))+1)-1
t.stat <- qt(0.975, df = df) #calculates critical t-value for the threshold (first value) and df (= n - p - 1)

# slope for PB.PP: 
SlPB1 <- fixef(modPB)[5]
SlPB1.l <- intervals(modPB)[1]
SlPB1.u <- intervals(modPB)[5,2]
SlPB1.l2 <- SlPB1 - t.stat * sqrt(vcov(modPB)[5,5])
SlPB1.u2 <- SlPB1 + t.stat * sqrt(vcov(modPB)[5,5])

# slope for PB.ZP: 
SlPB2 <- fixef(modPB)[5] + fixef(modPB)[6]
SlPB2.l <- SlPB2 - t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[6,6] + 2*vcov(modPB)[6,5])
SlPB2.u <- SlPB2 + t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[6,6] + 2*vcov(modPB)[6,5])

# slope for PB.PZN: 
SlPB3 <- fixef(modPB)[5] + fixef(modPB)[7]
SlPB3.l <- SlPB3 - t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[7,7] + 2*vcov(modPB)[7,5])
SlPB3.u <- SlPB3 + t.stat * sqrt(vcov(modPB)[5,5] + vcov(modPB)[7,7] + 2*vcov(modPB)[7,5])

slopesPB <- (cbind(c(SlPB1, SlPB2, SlPB3), c(SlPB1.l2, SlPB2.l, SlPB3.l), c(SlPB1.u2, SlPB2.u, SlPB3.u)))
rownames(slopesPB) <- c("P", "PZ", "PZN")
colnames(slopesPB) <- c("S", "l", "u")

# ints for PB.PP: 
IPB1 <- fixef(modPB)[1] - fixef(modPB)[5]*mean(data$invTT)
IPB1.l <- confint(modPB)[1,1]
IPB1.u <- confint(modPB)[1,2]
IPB1.l2 <- IPB1 - t.stat * sqrt(vcov(modPB)[1,1])

# ints for PB.ZP: modPB includes all int tPBms, but we leave out the invTi tPBm for among group lines; ints hPBe are at mean(invTT) 
IPB2 <- fixef(modPB)[1] + fixef(modPB)[3] - fixef(modPB)[5]*mean(data$invTT) - fixef(modPB)[6]*mean(data$invTT)
#IPB2.l <- IPB2 - t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[3,3] + 2*vcov(modPB)[3,1] + vcov(modPB)[5,5]*(mean(data$invTT)^2) + vcov(modPB)[6,6]*(mean(data$invTT)^2)) 
#IPB2.u <- IPB2 + t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[3,3] + 2*vcov(modPB)[3,1])

# ints for PB.PZN: 
IPB3 <- fixef(modPB)[1] + fixef(modPB)[4] - fixef(modPB)[5]*mean(data$invTT) - fixef(modPB)[7]*mean(data$invTT)
IPB3.l <- IPB3 - t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[4,4] + 2*vcov(modPB)[4,1]) 
IPB3.u <- IPB3 + t.stat * sqrt(vcov(modPB)[1,1] + vcov(modPB)[4,4] + 2*vcov(modPB)[4,1])

slopesPB <- (cbind(c(SlPB1, SlPB2, SlPB3), c(SlPB1.l, SlPB2.l, SlPB3.l), c(SlPB1.u, SlPB2.u, SlPB3.u)))
rownames(slopesPB) <- c("P", "PZ", "PZN")
colnames(slopesPB) <- c("S", "l", "u")




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





# BACKTRANSFORMED FIGURES -------------------------------------------------
NPP.BT <- ggplot(data = mod.coefs, aes(x = invTT, y = NPP2, ymin = -2)) + 
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(aes(group = Tank, shape = trophic.level), alpha = 1/2, size = 2) + #color = trophic.level
  xlab("Temperature") +
  ylab("NPPi")

NPP.BT

#y = m*x + b
#y = exp(b)*x^m

intervals(modNPP5r, which = "fixed")[1] -> ints
ints[[1]][1,1]
b <- mean(mod.coefs$invTT)

NPP.funcP <- function(x) { (fixef(modNPP5r)[1] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT)) + fixef(modNPP5r)[5]*x}

#NPP.btP <- function(x) { exp(fixef(modNPP5r)[1] - fixef(modNPP5r)[5]*mean(x)) * exp(x)^(fixef(modNPP5r)[5])}
yvalsP <- NPP.funcP(mod.coefs$invTT)

#NPP.btPlow <- function(x) { exp(ints[[1]][1,1] - ints[[1]][5,1]*mean(x)) * exp(x)^(ints[[1]][5,1])}
#yvalsPl <- NPP.btPlow(mod.coefs$invTT)
#NPP.btPhi <- function(x) { exp(ints[[1]][1,3] - ints[[1]][5,3]*mean(x)) * exp(x)^(ints[[1]][5,3])}
#yvalsPh <- NPP.btPhi(mod.coefs$invTT)

NPP.funcPZ <- function(x) { (fixef(modNPP5r)[1] + fixef(modNPP5r)[3] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT) - fixef(modNPP5r)[8]*mean(mod.coefs$invTT)) + (fixef(modNPP5r)[5] + fixef(modNPP5r)[8])*x}
yvalsPZ <- NPP.funcPZ(mod.coefs$invTT)

NPP.btPZ <- function(x) { exp(fixef(modNPP5r)[1] + fixef(modNPP5r)[3] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT) - fixef(modNPP5r)[8]*mean(mod.coefs$invTT)) * exp(x)^(fixef(modNPP5r)[5] + fixef(modNPP5r)[8]) }
yvalsPZ <- NPP.btPZ(mod.coefs$invTT)

NPP.funcPZN <- function(x) { (fixef(modNPP5r)[1] + fixef(modNPP5r)[4] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT) - fixef(modNPP5r)[9]*mean(mod.coefs$invTT)) + (fixef(modNPP5r)[5] + fixef(modNPP5r)[9])*x}
yvalsPZN <- NPP.funcPZN(mod.coefs$invTT)

NPP.btPZN <- function(x) { exp(fixef(modNPP5r)[1] + fixef(modNPP5r)[4] - fixef(modNPP5r)[5]*mean(mod.coefs$invTT) - fixef(modNPP5r)[9]*mean(mod.coefs$invTT)) * exp(x)^(fixef(modNPP5r)[5] + fixef(modNPP5r)[9]) }
yvalsPZN <- NPP.btPZN(mod.coefs$invTT)

NPP.BT +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsP), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) + 
  # geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPl), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey60', size = 1.5) +
  # geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPh), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey60', size = 1.5) # not sure if this is right.
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZ), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey80', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = invTT, y = yvalsPZN), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey60', size = 1.5) +
  
  #  geom_segment(aes(x = min(invTT), xend = max(invTT), y = exp(intervals(modNPP5r, which = "fixed")[[1]][1,1]), yend = exp(intervals(modNPP5r, which = "fixed")[[1]][1,3])), linetype=1, color = "gray40")
  # trying to plot the curves with CIs. For the fixefs we can get these CIs from intervals (below). then, I suppose to plot them I'd make a different line - is it just the line defined by the parameters at the lower edge of their CI, and another line above?
  #thinking through this: we log transoformed the data, then fit a model to estimate temperature effects. I want to describe this temperature effect (slope parameter) and our confidence in it. Confidence will come from the confidence intervals on the fixed effect (for P, this is just the centered temp effect), but for the PZ level, it's that term + the base term... with which uncertainty?
  
  intervals(modNPP5r, which = "fixed")[1] -> ints
ints[[1]][1,1]

ggsave("figure2Aalt2.png", device = "png", width = 5, height = 3)



ER.BT <- ggplot(data = mod.coefsER, aes(x = invTT, y = ER2)) + 
  theme_bw() +
  theme(legend.position = "none") +
  #geom_point(aes(group = Tank, shape = trophic.level), alpha = 1/2, size = 2) + #color = trophic.level
  xlab("Temperature") +
  ylab("ERi")

ER.BT

ER.btP <- function(x) { exp(fixef(modER2r)[1] - fixef(modER2r)[3]*mean(x)) * exp(x)^(fixef(modER2r)[3]) } # for trophic level 1
yvalsP <- ER.btP(mod.coefsER$invTT)

ER.btPZ <- function(x) { exp(fixef(modER2r)[1] + fixef(modER2r)[4] - fixef(modER2r)[3]*mean(x)) * exp(x)^(fixef(modER2r)[3]) } # for trophic level 2
yvalsPZ <- ER.btPZ(mod.coefsER$invTT)

ER.btPZN <- function(x) { exp(fixef(modER2r)[1] + fixef(modER2r)[5] - fixef(modER2r)[3]*mean(x)) * exp(x)^(fixef(modER2r)[3]) } # for trophic level 3
yvalsPZN <- ER.btPZN(mod.coefsER$invTT)


ER.BT +
  geom_smooth(data = mod.coefsER, aes(x = invTT, y = yvalsP), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefsER, aes(x = invTT, y = yvalsPZ), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey80', size = 1.5) +
  geom_smooth(data = mod.coefsER, aes(x = invTT, y = yvalsPZN), method = loess, se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey60', size = 1.5)

ggsave("figure2Balt2.png", device = "png", width = 5, height = 3)


# FIGURE 2D ---------------------------------------------------------------

## 
## Patrick's figure idea: dplyr
NPPmod = mod.coefs %>% 
  gather(key = "EF", value = "Rate", NPP2, ER2, PP.biomass) %>%
  group_by(Tank, trophic.level, EF) %>% 
  do(fitTank = lm(log(Rate) ~ invTi, data = .)) %>%
  tidy(., fitTank) %>%
  filter(term == 'invTi') %>%
  select(trophic.level, estimate) %>%
  mutate(level = "Tank")

xtanks <- data.frame(trophic.level = c("P", "PZ", "PZN"), estimate = c(-1.31564, -0.8241136, -0.4654028, fixef(modER2r)[3], fixef(modER2r)[3], fixef(modER2r)[3], 1.281337, 3.669879, 1.735196), Tank = "all", level = "treatment", EF = rep(c("NPP2", "ER2", "PP.biomass"), each = 3))

bind_rows(NPPmod, xtanks) -> data3

plot1 <- ggplot(data3, aes(x = trophic.level, y = estimate, shape = level)) +
  scale_y_continuous(name ="Activation Energy (Ea)", labels=c("-4","-2","0", "2", "4"), limits = c(-4, 4.2)) +
  xlab("Food Chain Length") +
  geom_point(size = 3, color = "gray50") +
  geom_point(data = subset(data3, level == "treatment"), size = 4) +
  guides(shape=guide_legend(title=NULL)) +
  scale_shape_discrete(breaks = c("Tank", "treatment"), labels = c("Within Ecosystem", "Across Ecosystems")) +
  facet_grid(~EF) +
  theme_bw() +
  theme(legend.position="top")
plot1
ggsave("figure1D.png", device = "png", width = 5, height = 3)

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


