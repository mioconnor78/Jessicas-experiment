#################################################
#### Garzke, O'Connor and Sommer temperature experiment
#### File for creating size-corrected biomass estimates. Necessary for mixed effects models.
#### Code by Mary O'Connor, revised November 2015
#################################################



### size-corrected biomass analysis

library(plyr)
#library(broom)
#library(dplyr)

### following Barneche et al 2014

### for each tank, we will want: 

### size corrected biomass

sdata <- read.csv("./zooplankton size for correction no eggs.csv")
names(sdata) <- c('treatment', 'tank', 'week', 'species', 'stage', 'size')
sdata$drywt <- ifelse(sdata$species == 'Daphnia', exp(1.468 + 2.829*log(sdata$size)), exp(1.821 + 0.654*log(sdata$size)))
sdata$C.ug <- sdata$drywt*0.5

#add in notonectids
sdata$noto.wt.corr <- ifelse(sdata$treatment == 'PZN', 2*((10.8/4)^.75), 0) #adding in a mass-corrected biomass estimate for notonectids
mass.data <- ddply(sdata, .(week, tank, treatment), summarise, sum((C.ug^0.75))) # estimating a mass corrected biomass for adult zooplankton in all tanks.
mass.data2 <- ddply(sdata, .(week, tank, treatment), summarise, sum((C.ug^0.75)+noto.wt.corr)) # estimating a mass corrected biomass for adult zooplankton incl Noto in all tanks.
mass.data3 <- ddply(sdata, .(week, tank, treatment), summarise, sum((C.ug))) 
mass.data$M.corrz <- mass.data2$..1 # corrected total zooplankton
mass.data$zpC.uncorr <- mass.data3$..1  # uncorrected total (no notonectids)
names(mass.data) <- c('week', 'tank', 'treatment', 'M.corrT', 'M.corrz', 'M.uncorrz')


## estimating a ratio of metabolic biomass to total biomass. T is total biomass.
msizefunc <- function(n, T) n*((T/n)^0.75)
plot(seq(1,60,1), msizefunc(seq(1,60,1), 12), xlim=c(0,60), ylim = c(0,60), xlab = 'number of individuals', ylab = 'metabolic biomass')
abline(h = 12)
abline(v = 12)
