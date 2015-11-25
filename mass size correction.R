#################################################
#### Garzke, O'Connor and Sommer temperature experiment
#### File for creating size-corrected biomass estimates. Necessary for mixed effects models.
#### Code by Mary O'Connor, revised November 2015
#################################################



### size-corrected biomass analysis

library(plyr)

### following Barneche et al 2014

### for each tank, we will want: 

### size corrected biomass

sdata <- read.csv("./zooplankton size for correction no eggs.csv")
sdata$drywt <- ifelse(sdata$species == 'Daphnia', exp(log(1.468 + 2.829*log(sdata$size))), exp(log(1.821 + 0.654*log(sdata$size))))
sdata$C.ug <- sdata$drywt*0.5

#add in notonectids
sdata$noto.wt.corr <- ifelse(sdata$treatment == 'PZN', 2*((10.8/4)^.75), 0) #adding in a mass-corrected biomass estimate for notonectids
mass.data <- ddply(sdata, .(week, tank, treatment), summarise, sum((C.ug^0.75))) # estimating a mass corrected biomass for adult zooplankton in all tanks.
mass.data2 <- ddply(sdata, .(week, tank, treatment), summarise, sum((C.ug^0.75)+noto.wt.corr)) # estimating a mass corrected biomass for adult zooplankton incl Noto in all tanks.
mass.data3 <- ddply(sdata, .(week, tank, treatment), summarise, sum((C.ug))) 
mass.data$M.corrz <- mass.data2$..1
mass.data$zpC.uncorr <- mass.data3$..1
names(mass.data) <- c('week', 'tank', 'treatment', 'M.corrT', 'M.corrz', 'M.uncorrz')

