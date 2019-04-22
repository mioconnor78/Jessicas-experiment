## size data estimation
## Garzke et al
## code written by M. O'Connor, Nov 3 2018

library(tidyverse) 
library(dplyr)

pptaxa <- read.csv("./data/pptaxa.csv")
ppinfo <- read.csv("./data/phytoinfo.csv")
View(pptaxa)
View(ppinfo)

## clean up names
pptaxa <- pptaxa %>%
  mutate(`taxon` = str_replace(`taxon`, "Coleosphaerium", "Coelosphaerium"))

pptaxa$taxon <- as.factor(pptaxa$taxon)

## week numbers in this file may not align with our other data. When i used the original week numbers, there was no relationship between chla (from other datafile) and estimated cell volume and size. Based on stephanie's report, I think the week numbers here need to be adjusted by one week. She reported that week 3 = July 10. July 10 was week 2 in our other datafile. So I will adjust the week numbers here. 

colnames(pptaxa)[1] <- 'weekSC'

pptaxa <- pptaxa %>%
  mutate(week = weekSC - 1)
  
## create a cell mass column for each taxon. first in ppinfo, then it can be joined to pptaxa
## first will do it by estimating volume from length information - done below, now check 
## can come back later to convert volume to biomass (joey has done this)

sphere = function(d) (4/3)*pi*((d/2)^3)
cylinder = function(d, l)  pi*((d/2)^2)*l
oval = function(d, l, w) (4/3)*pi*(d*l*w)


## cell vol is in units of um^3
  ppinfo <- ppinfo %>%
    mutate(., d = ifelse(shape == "sphere", 0.5*(min.length.um + max.length.um), min.length.um)) %>%
    mutate(., cellvol = ifelse(shape == "sphere", sphere(d), cylinder(d, max.length.um))) %>%
    mutate(., cellvol = ifelse(shape == "oval", oval(d/2, max.length.um/2, d/2), cellvol)) %>%
    mutate(., cellvol = ifelse(taxon == "Chroococcus", oval(10/2, max.length.um/2, 10/2), cellvol)) 
  
  #%>% #produces volum in um^3
    # mutate(., colvol = ifelse(colonial != "no", oval(d, max.length.um, d), vol)) #need to go back and check whether colonies are sheets, spheres or filaments
    
    
## merge cell vols with pptaxa
pptaxa1 <- pptaxa %>%
    left_join(., ppinfo, by = c("taxon", "group")) %>%
  mutate(., cellmass = 0.209*cellvol^0.991)
 
## create vector of sizes for each tank and date
pptaxa2 <- pptaxa1 %>%
  filter(., ind.L > 0)

pptaxa3 <- pptaxa2[rep(seq(nrow(pptaxa2)), pptaxa2$ind.L),]   

## ok, now if we group by week and tank, we will have the list of sizes
## using functions from below in the 'proof of concept' part
a <- 0.75
Mb = function(x) sum(x) #total biomass
mba <- function(x) sum(x^(a))/length(x) # average body size that accounts for size dependent changes in metabolic rate, based on a
mb <- function(x) sum(x)/length(x)
mb2 <- mba <- function(x) sum(x^(a))
mba1 <- function(x) sum(x^(a-1)) # average body size that accounts for size dependent changes in metabolic rate, avg weighted by biomass rather than density (as in mba); equivalent to Mb = density/vol * (mba)
mtriangle <- function(x) mb2(x)/Mb(x)

## so here we have total biomass, two estimates of average body size... do i have the total mass corrected size?

## remove ciliates; note later that we did see them. 

pptaxaMb <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., Mb = Mb(cellmass))

pptaxamba <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., mba = mba(cellmass))

pptaxamba1 <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., mba1 = mba1(cellmass))

pptaxaN <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., pp.N = length(cellmass))

pptaxaAvgSize <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., AvgS = mean(cellmass))
  

## do it for cyanos and non-cyanos separately
Cyano.mba1 <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group == "Cyanobacteria") %>%
  dplyr::summarize(., Cmba1 = mba1(cellmass))

Cyano.Mb <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group == "Cyanobacteria") %>%
  dplyr::summarize(., CMb = Mb(cellmass))

Cyano.mba <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group == "Cyanobacteria") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., Cmba = mba(cellmass))

Cyano.N <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group == "Cyanobacteria") %>%
  dplyr::summarize(., C.N = length(cellmass))

Cyano.AvgS <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group == "Cyanobacteria") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., CAvgS = mean(cellmass))

Oth.mba1 <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Cyanobacteria") %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., Omba1 = mba1(cellmass))

Oth.Mb <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Cyanobacteria") %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., OMb = Mb(cellmass))

Oth.mba <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Cyanobacteria") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., Omba = mba(cellmass))

Oth.N <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Cyanobacteria") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., O.N = length(cellmass))

Oth.AvgS <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Cyanobacteria") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., OAvgS = mean(cellmass))

CN.mba1 <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(colonial == "no") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., CNmba1 = mba1(cellmass))

CN.Mb <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(colonial == "no") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., CNMb = Mb(cellmass))

CN.mba <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(colonial == "no") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., CNmba = mba(cellmass))

CN.N <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(colonial == "no") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., CN.N = length(cellmass))

CN.AvgS <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(colonial == "no") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., CN.AvgS = mean(cellmass))

FN.mba1 <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(shape != "Filament") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., FNmba1 = mba1(cellmass))

FN.Mb <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(shape != "Filament") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., FNMb = Mb(cellmass))

FN.mba <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(shape != "Filament") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., FNmba = mba(cellmass))

FN.N <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(shape != "Filament") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., FN.N = length(cellmass))

FN.AvgS <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(shape != "Filament") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  dplyr::summarize(., FN.AvgS = mean(cellmass))

size.data <- pptaxaMb %>%
  left_join(., pptaxamba, by = c("week", "tank")) %>%
  left_join(., pptaxamba1, by = c("week", "tank")) %>%
  left_join(., pptaxaN, by = c("week", "tank")) %>%
  left_join(., pptaxaAvgSize, by = c("week", "tank")) %>%
  left_join(., Cyano.Mb, by = c("week", "tank")) %>%
  left_join(., Cyano.mba, by = c("week", "tank")) %>%
  left_join(., Cyano.mba1, by = c("week", "tank")) %>%
  left_join(., Cyano.N, by = c("week", "tank")) %>%
  left_join(., Cyano.AvgS, by = c("week", "tank")) %>%
  left_join(., Oth.Mb, by = c("week", "tank")) %>%
  left_join(., Oth.mba, by = c("week", "tank")) %>%
  left_join(., Oth.mba1, by = c("week", "tank")) %>%
  left_join(., Oth.N, by = c("week", "tank")) %>%
  left_join(., Oth.AvgS, by = c("week", "tank")) %>%
  left_join(., CN.Mb, by = c("week", "tank")) %>%
  left_join(., CN.mba, by = c("week", "tank")) %>%
  left_join(., CN.mba1, by = c("week", "tank")) %>%
  left_join(., CN.N, by = c("week", "tank")) %>%
  left_join(., CN.AvgS, by = c("week", "tank")) %>%
  left_join(., FN.Mb, by = c("week", "tank")) %>%
  left_join(., FN.mba, by = c("week", "tank")) %>%
  left_join(., FN.mba1, by = c("week", "tank")) %>%
  left_join(., FN.N, by = c("week", "tank")) %>%
  left_join(., FN.AvgS, by = c("week", "tank"))

## looking at oddball values
hist(size.data$CMb)
dplyr::filter(size.data, CMb > 100000)
week1 <- pptaxa %>%
  dplyr::filter(., week == 1) %>%
  dplyr::filter(., tank == 14)

week2 <- pptaxa %>%
  dplyr::filter(., week == 2) %>%
  dplyr::filter(., tank == c(19,25))
## no clear reason to exclude them

hist(size.data$Mb)
hist(size.data$OMb)

write.csv(size.data, file = "PPsizes.csv")


plot(size.data$mba ~ size.data$Mb)
plot(size.data$mba1 ~ size.data$mba)
plot(size.data$Mb ~ size.data$pp.N)



plot(size.data$mba ~ size.data$Mb)
plot(size.data$mba1 ~ size.data$mba)
plot(size.data$Oth.Mb ~ size.data$Oth.N)

##### analysis of mean size
pptaxa4 <- as.tibble(pptaxa3) %>%
  group_by(week) %>%
  left_join(., temps.pwr, by = c("week", "power")) %>%
  mutate(., invTavg = 1/((avgTemp + 273)*k)) %>%
  filter(., week != "0") %>%
  filter(., week != "1") %>%
  filter(group != "Ciliate") %>%
  filter(group != "Amoeboid") %>%
  filter(shape != "filament")

modS1 <- lme(cellmass ~ 1 + invTavg*trophic.level, random = ~1|tank, data = pptaxa4, na.action = na.omit)
modS2 <- lme(cellmass ~ 1 + invTavg + trophic.level, random = ~1|tank, data = pptaxa4, na.action = na.omit)
modS3 <- lme(cellmass ~ 1 + invTavg, random = ~1|tank, data = pptaxa4, na.action = na.omit)
modS4 <- lme(cellmass ~ 1 + trophic.level, random = ~1|tank, data = pptaxa4, na.action = na.omit)
modS5 <- lme(cellmass ~ 1, random = ~1|tank, data = pptaxa4, na.action = na.omit)

modSres <- data.frame(model.sel(modS1, modS2, modS3, modS4, modS5))

plot(cellmass ~ invTavg, data = pptaxa4)
lines(a = , b = 181.5)

hist(log(pptaxa4$cellmass))
hist(log(pptaxa4[(pptaxa4$trophic.level == "P"), ]$cellmass))
hist(log(pptaxa4[(pptaxa4$trophic.level == "PZ"), ]$cellmass), hold = TRUE)

## this isn't what i wanted, leaving it now
ggplot(data = pptaxa4, aes(x = invTavg)) + 
  geom_histogram(data = subset(pptaxa4, trophic.level == "P"), fill = "red", alpha = 0.2) + 
  geom_histogram(data = subset(pptaxa4, trophic.level == "PZ"), fill = "blue", alpha = 0.2) +
  geom_histogram(data = subset(pptaxa4, trophic.level == "PZN"), fill = "green", alpha = 0.2)

PPsize <- ggplot(data = pptaxa4, aes(x = invTavg, y = log(cellmass)), group = week, shape = as.factor(week)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(aes(group = as.character(trophic.level), color = as.character(trophic.level), shape = as.factor(week)), size = 2)
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),  
        legend.text = element_text(size=6), 
        legend.title = element_text(size = 7), 
        strip.background = element_blank()) +
  ylab("Phytoplankton Cell Mass \n ln(ug)")

PPsize




##### metabolic biomass
## proof of concept:
sizes <- c(2,2,3,2.5,1.2, 1, 1.1, 1.2, 1, 4)
sizes2 <- c(7, 5, 5, 2)
sizes3 <- c(.6, .5, .4, .6, .3, .2, 1, 1.2, 2, .2, .3, .4, .3, .5, .1, 2, .1, 1, 4, 3, .3)
sizes4 <- c(19)

a <- 0.75
Mb = function(x) sum(x) #total biomass
mba <- function(x) sum(x^(a))/length(x) # average body size that accounts for size dependent changes in metabolic rate, based on a
mb <- function(x) sum(x)/length(x)
mba1 <- function(x) sum(x^(a-1)) # average body size that accounts for size dependent changes in metabolic rate, avg weighted by biomass rather than density (as in mba); equivalent to Mb = density/vol * (mba)

#mass correction is:
mba(sizes)/mb(sizes)
mba(sizes4)/mb(sizes4)

#mass corrected biomass is: 
Mb(sizes)*(mba(sizes)/mb(sizes))
########


#Mb <- PPbiomass (chlorophyll*carbon...) = sum(biovolumes), 
### ok, i have gotten sizes, coloniality, etc, so now have to bring the data in here and estimate volumes, taking into account shape. 
## step 1: vector of biomasses:


Mb <- sum()
#Mba <- 
cells.tot <- pptaxa %>%
  group_by(tank, week) %>%
  dplyr::summarize(cells.L = sum(ind.L*100, na.rm = TRUE)) #%>% #data is for ind/10mL





#refs https://books.google.ca/books?hl=en&lr=&id=C2QXv266DYQC&oi=fnd&pg=PA189&dq=%22cell+volume%22+Chrysosphaerella&ots=mYsIyCqSf-&sig=96fYamPVph0txxfgGIKYXq0PHzE#v=onepage&q&f=false

