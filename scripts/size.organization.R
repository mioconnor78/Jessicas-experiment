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
  
## create a cell mass column for each taxon. first in ppinfo, then it can be joined to pptaxa
## first will do it by estimating volume from length information - done below, now check 
## can come back later to convert volume to biomass (joey has done this)

sphere = function(d) (4/3)*pi*((d/2)^3)
cylinder = function(d, l)  pi*((d/2)^2)*l
oval = function(d, l, w) (4/3)*pi*(d*l*w)


## cell vol is in units of um^3
## these are seeming ok...
  ppinfo <- ppinfo %>%
    mutate(., d = ifelse(shape == "sphere", 0.5*(min.length.um + max.length.um), min.length.um)) %>%
    mutate(., cellvol = ifelse(shape == "sphere", sphere(d), cylinder(d, max.length.um))) %>%
    mutate(., cellvol = ifelse(shape == "oval", oval(d/2, max.length.um/2, d/2), cellvol)) #%>% #produces volum in um^3
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
mba1 <- function(x) sum(x^(a-1)) # average body size that accounts for size dependent changes in metabolic rate, avg weighted by biomass rather than density (as in mba); equivalent to Mb = density/vol * (mba)

## remove ciliates; note later that we did see them. 

pptaxaMb <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., Mb = Mb(cellmass))

pptaxamba <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., mba = mba(cellmass))

pptaxamba1 <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., mba1 = mba1(cellmass))

pptaxaN <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., pp.N = length(cellmass))

pptaxaAvgSize <- pptaxa3 %>%
  group_by(tank, week) %>%
  filter(group != "Ciliate") %>%
  dplyr::summarize(., AvgS = mean(cellmass))
  
size.data <- pptaxaMb %>%
  left_join(., pptaxamba, by = c("week", "tank")) %>%
  left_join(., pptaxamba1, by = c("week", "tank")) %>%
  left_join(., pptaxaN, by = c("week", "tank")) %>%
  left_join(., pptaxaAvgSize, by = c("week", "tank"))
  
plot(size.data$mba ~ size.data$Mb)
plot(size.data$mba1 ~ size.data$mba)
plot(size.data$Mb ~ size.data$pp.N)

write.csv(size.data, file = "PPsizes.csv")

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

