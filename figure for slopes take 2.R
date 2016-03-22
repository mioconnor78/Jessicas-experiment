## making a figure for slopes only

## Figure X. Estimated activation energies for ecosystem carbon stocks and oxygen fluxes under three trophic conditions. Mean activation energies are estimated by taking the average of each weekly estimated activation energy (fixed effect + weekly random effect, so n = 6). Because the estimated activation energies are derived from the MEM, non-independence of treatments within weeks is considered in the random effect.

## So in practice: get the estimate (fixed + random) for each tank at each time; then take the means (and SEMs), and plot those in this figure.

library(broom)
library(dplyr)

se <- function(x) sd(x)/sqrt(length(x))

############
#modNPP4r
modNPP <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)  #modNPP4r


###Figure 2a, slopes
rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(NPP2))
names(rand.cat) <- c('Week', 'Tank', 'invT','TL', 'NPP2')
Entry.coefs <- data.frame(coef(modNPP)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted intercepts for trophic groups
S$ITL.term <- ifelse(S$TL == 'P', S$X.Intercept., 0)
S$ITL.term <- ifelse(S$TL == 'PZ', (S$X.Intercept. + S$trophic.levelPZ), S$ITL.term)
S$ITL.term <- ifelse(S$TL == 'PZN', (S$X.Intercept. + S$trophic.levelPZN), S$ITL.term)


#add week ranefs
ranefs <- data.frame(ranef(modNPP)$week)
ranefs$week <- as.numeric(rownames(ranefs))
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)
names(S2) <- c('Week',"Tank", "invT", "TL", "NPP2","X.Intercept..x","I.invT...mean.invT...x", "trophic.levelPZ" ,"trophic.levelPZN","I.invT...mean.invT...trophic.levelPZ" ,"I.invT...mean.invT...trophic.levelPZN", "TL.term" ,"ITL.term","Int.ranef" ,"sl.ranef")

## summary with appropriate number of observations: 
NPP.TLs <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
NPP.TLi <- ddply(S2, .(Week, TL), summarize, mean(ITL.term))
NPP.ranefs <- ddply(S2, .(Week, TL), summarize, mean(sl.ranef))
NPP.ranefsi <- ddply(S2, .(Week, TL), summarize, mean(Int.ranef))
NPP.S <- merge(NPP.TLs, NPP.ranefs, by = c('Week', 'TL'), all= FALSE)
NPP.I <- merge(NPP.TLi, NPP.ranefsi, by = c('Week', 'TL'), all= FALSE)
NPPs <- merge(NPP.S, NPP.I, by = c('Week', 'TL'), all = FALSE)
names(NPPs) <- c("Week","TL","TL.term","sl.ranef","ITL.term","Int.ranef")
NPPs$Ea <- NPPs$TL.term + NPPs$sl.ranef
NPPs$Bo <- NPPs$ITL.term + NPPs$Int.ranef

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumNPP <- ddply(NPPs, .(TL), summarize, mean(Ea))
Ea.sumNPP2 <- ddply(NPPs, .(TL), summarize, se(Ea))
Ea.sumNPPs <- merge(Ea.sumNPP, Ea.sumNPP2, by = 'TL')
names(Ea.sumNPPs) <- c('TL', 'means','se')

Bo.sumNPP <- ddply(NPPs, .(TL), summarize, mean(Bo))
Bo.sumNPP2 <- ddply(NPPs, .(TL), summarize, se(Bo))
Bo.sumNPPs <- merge(Bo.sumNPP, Bo.sumNPP2, by = 'TL')
names(Bo.sumNPPs) <- c('TL', 'means','se')


## repeat for ES
#################

#modERr
modER4r <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data[(data$ER2 >= 0),], REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(ER2))
names(rand.cat) <- c('Week', 'Tank',  'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modER4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted intercepts for trophic groups
S$ITL.term <- ifelse(S$TL == 'P', S$X.Intercept., 0)
S$ITL.term <- ifelse(S$TL == 'PZ', (S$X.Intercept. + S$trophic.levelPZ), S$ITL.term)
S$ITL.term <- ifelse(S$TL == 'PZN', (S$X.Intercept. + S$trophic.levelPZN), S$ITL.term)

#add week ranefs
ranefs <- data.frame(ranef(modER4r)$week)
ranefs$week <- as.numeric(rownames(ranefs))
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)
names(S2) <- c('Week',"Tank", "invT", "TL", "NPP2","X.Intercept..x","I.invT...mean.invT...x", "trophic.levelPZ" ,"trophic.levelPZN","I.invT...mean.invT...trophic.levelPZ" ,"I.invT...mean.invT...trophic.levelPZN", "TL.term" ,"ITL.term","Int.ranef" ,"sl.ranef")

## summary with appropriate number of observations: 
ER.TLs <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
ER.TLi <- ddply(S2, .(Week, TL), summarize, mean(ITL.term))
ER.ranefs <- ddply(S2, .(Week, TL), summarize, mean(sl.ranef))
ER.ranefsi <- ddply(S2, .(Week, TL), summarize, mean(Int.ranef))
ER.S <- merge(ER.TLs, ER.ranefs, by = c('Week', 'TL'), all= FALSE)
ER.I <- merge(ER.TLi, ER.ranefsi, by = c('Week', 'TL'), all= FALSE)
ERs <- merge(ER.S, ER.I, by = c('Week', 'TL'), all = FALSE)
names(ERs) <- c("Week","TL","TL.term","sl.ranef","ITL.term","Int.ranef")
ERs$Ea <- ERs$TL.term + ERs$sl.ranef
ERs$Bo <- ERs$ITL.term + ERs$Int.ranef

Ea.sumER <- ddply(ERs, .(TL), summarize, mean(Ea))
Ea.sumER2 <- ddply(ERs, .(TL), summarize, se(Ea))
Ea.sumERs <- merge(Ea.sumER, Ea.sumER2, by = 'TL')
names(Ea.sumERs) <- c('TL', 'means','se')

Bo.sumER <- ddply(ERs, .(TL), summarize, mean(Bo))
Bo.sumER2 <- ddply(ERs, .(TL), summarize, se(Bo))
Bo.sumERs <- merge(Bo.sumER, Bo.sumER2, by = 'TL')
names(Bo.sumERs) <- c('TL', 'means','se')

## repeat for NPPm
#################

#modNPPm4r
modNPPm4r <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(NPP.mass))
names(rand.cat) <- c('Week', 'Tank','invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modNPPm4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted intercepts for trophic groups
S$ITL.term <- ifelse(S$TL == 'P', S$X.Intercept., 0)
S$ITL.term <- ifelse(S$TL == 'PZ', (S$X.Intercept. + S$trophic.levelPZ), S$ITL.term)
S$ITL.term <- ifelse(S$TL == 'PZN', (S$X.Intercept. + S$trophic.levelPZN), S$ITL.term)


#add week ranefs
ranefs <- data.frame(ranef(modNPPm4r)$week)
ranefs$week <- as.numeric(rownames(ranefs))
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)
names(S2) <- c('Week',"Tank", "invT", "TL", "NPP2","X.Intercept..x","I.invT...mean.invT...x", "trophic.levelPZ" ,"trophic.levelPZN","I.invT...mean.invT...trophic.levelPZ" ,"I.invT...mean.invT...trophic.levelPZN", "TL.term" ,"ITL.term","Int.ranef" ,"sl.ranef")

## summary with appropriate numbNPPm of obsNPPmvations: 
NPPm.TLs <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
NPPm.TLi <- ddply(S2, .(Week, TL), summarize, mean(ITL.term))
NPPm.ranefs <- ddply(S2, .(Week, TL), summarize, mean(sl.ranef))
NPPm.ranefsi <- ddply(S2, .(Week, TL), summarize, mean(Int.ranef))
NPPm.S <- merge(NPPm.TLs, NPPm.ranefs, by = c('Week', 'TL'), all= FALSE)
NPPm.I <- merge(NPPm.TLi, NPPm.ranefsi, by = c('Week', 'TL'), all= FALSE)
NPPms <- merge(NPPm.S, NPPm.I, by = c('Week', 'TL'), all = FALSE)
names(NPPms) <- c("Week","TL","TL.term","sl.ranef","ITL.term","Int.ranef")
NPPms$Ea <- NPPms$TL.term + NPPms$sl.ranef
NPPms$Bo <- NPPms$ITL.term + NPPms$Int.ranef

Ea.sumNPPm <- ddply(NPPms, .(TL), summarize, mean(Ea))
Ea.sumNPPm2 <- ddply(NPPms, .(TL), summarize, se(Ea))
Ea.sumNPPms <- merge(Ea.sumNPPm, Ea.sumNPPm2, by = 'TL')
names(Ea.sumNPPms) <- c('TL', 'means','se')

Bo.sumNPPm <- ddply(NPPms, .(TL), summarize, mean(Bo))
Bo.sumNPPm2 <- ddply(NPPms, .(TL), summarize, se(Bo))
Bo.sumNPPms <- merge(Bo.sumNPPm, Bo.sumNPPm2, by = 'TL')
names(Bo.sumNPPms) <- c('TL', 'means','se')


## repeat for ERm
#################

#modERm4r
modERm4r <- lmer(log(ER.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data[(data$ER2 >= 0),], REML = TRUE, na.action=na.omit) 

rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(ER.mass))
names(rand.cat) <- c('Week', 'Tank', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modERm4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted intercepts for trophic groups
S$ITL.term <- ifelse(S$TL == 'P', S$X.Intercept., 0)
S$ITL.term <- ifelse(S$TL == 'PZ', (S$X.Intercept. + S$trophic.levelPZ), S$ITL.term)
S$ITL.term <- ifelse(S$TL == 'PZN', (S$X.Intercept. + S$trophic.levelPZN), S$ITL.term)


#add week ranefs
ranefs <- data.frame(ranef(modERm4r)$week)
ranefs$week <- as.numeric(rownames(ranefs))
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)
names(S2) <- c('Week',"Tank", "invT", "TL", "ERm2","X.Intercept..x","I.invT...mean.invT...x", "trophic.levelPZ" ,"trophic.levelPZN","I.invT...mean.invT...trophic.levelPZ" ,"I.invT...mean.invT...trophic.levelPZN", "TL.term" ,"ITL.term","Int.ranef" ,"sl.ranef")

## summary with appropriate number of observations: 
ERm.TLs <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
ERm.TLi <- ddply(S2, .(Week, TL), summarize, mean(ITL.term))
ERm.ranefs <- ddply(S2, .(Week, TL), summarize, mean(sl.ranef))
ERm.ranefsi <- ddply(S2, .(Week, TL), summarize, mean(Int.ranef))
ERm.S <- merge(ERm.TLs, ERm.ranefs, by = c('Week', 'TL'), all= FALSE)
ERm.I <- merge(ERm.TLi, ERm.ranefsi, by = c('Week', 'TL'), all= FALSE)
ERms <- merge(ERm.S, ERm.I, by = c('Week', 'TL'), all = FALSE)
names(ERms) <- c("Week","TL","TL.term","sl.ranef","ITL.term","Int.ranef")
ERms$Ea <- ERms$TL.term + ERms$sl.ranef
ERms$Bo <- ERms$ITL.term + ERms$Int.ranef

Ea.sumERm <- ddply(ERms, .(TL), summarize, mean(Ea))
Ea.sumERm2 <- ddply(ERms, .(TL), summarize, se(Ea))
Ea.sumERms <- merge(Ea.sumERm, Ea.sumERm2, by = 'TL')
names(Ea.sumERms) <- c('TL', 'means','se')

Bo.sumERm <- ddply(ERms, .(TL), summarize, mean(Bo))
Bo.sumERm2 <- ddply(ERms, .(TL), summarize, se(Bo))
Bo.sumERms <- merge(Bo.sumERm, Bo.sumERm2, by = 'TL')
names(Bo.sumERms) <- c('TL', 'means','se')



#modBBr
modPP4r <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(PP.biomass))
names(rand.cat) <- c('Week', 'Tank', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modPP4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted intercepts for trophic groups
S$ITL.term <- ifelse(S$TL == 'P', S$X.Intercept., 0)
S$ITL.term <- ifelse(S$TL == 'PZ', (S$X.Intercept. + S$trophic.levelPZ), S$ITL.term)
S$ITL.term <- ifelse(S$TL == 'PZN', (S$X.Intercept. + S$trophic.levelPZN), S$ITL.term)

#add week ranefs
ranefs <- data.frame(ranef(modPP4r)$week)
ranefs$week <- as.numeric(rownames(ranefs))
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)
names(S2) <- c('Week',"Tank", "invT", "TL", "PP2","X.Intercept..x","I.invT...mean.invT...x", "trophic.levelPZ" ,"trophic.levelPZN","I.invT...mean.invT...trophic.levelPZ" ,"I.invT...mean.invT...trophic.levelPZN", "TL.term" ,"ITL.term","Int.ranef" ,"sl.ranef")

## summary with appropriate number of observations: 
PP.TLs <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
PP.TLi <- ddply(S2, .(Week, TL), summarize, mean(ITL.term))
PP.ranefs <- ddply(S2, .(Week, TL), summarize, mean(sl.ranef))
PP.ranefsi <- ddply(S2, .(Week, TL), summarize, mean(Int.ranef))
PP.S <- merge(PP.TLs, PP.ranefs, by = c('Week', 'TL'), all= FALSE)
PP.I <- merge(PP.TLi, PP.ranefsi, by = c('Week', 'TL'), all= FALSE)
PPs <- merge(PP.S, PP.I, by = c('Week', 'TL'), all = FALSE)
names(PPs) <- c("Week","TL","TL.term","sl.ranef","ITL.term","Int.ranef")
PPs$Ea <- PPs$TL.term + PPs$sl.ranef
PPs$Bo <- PPs$ITL.term + PPs$Int.ranef

Ea.sumPP <- ddply(PPs, .(TL), summarize, mean(Ea))
Ea.sumPP2 <- ddply(PPs, .(TL), summarize, se(Ea))
Ea.sumPPs <- merge(Ea.sumPP, Ea.sumPP2, by = 'TL')
names(Ea.sumPPs) <- c('TL', 'means','se')

Bo.sumPP <- ddply(PPs, .(TL), summarize, mean(Bo))
Bo.sumPP2 <- ddply(PPs, .(TL), summarize, se(Bo))
Bo.sumPPs <- merge(Bo.sumPP, Bo.sumPP2, by = 'TL')
names(Bo.sumPPs) <- c('TL', 'means','se')


###################
#modTCr
modTC4r <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(total.carbon))
names(rand.cat) <- c('Week', 'Tank', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modTC4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

## construct predicted intercepts for trophic groups
S$ITL.term <- ifelse(S$TL == 'P', S$X.Intercept., 0)
S$ITL.term <- ifelse(S$TL == 'PZ', (S$X.Intercept. + S$trophic.levelPZ), S$ITL.term)
S$ITL.term <- ifelse(S$TL == 'PZN', (S$X.Intercept. + S$trophic.levelPZN), S$ITL.term)

#add week ranefs
ranefs <- data.frame(ranef(modTC4r)$week)
ranefs$week <- as.numeric(rownames(ranefs))
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)
names(S2) <- c('Week',"Tank", "invT", "TL", "TC2","X.Intercept..x","I.invT...mean.invT...x", "trophic.levelPZ" ,"trophic.levelPZN","I.invT...mean.invT...trophic.levelPZ" ,"I.invT...mean.invT...trophic.levelPZN", "TL.term" ,"ITL.term","Int.ranef" ,"sl.ranef")

## summary with appropriate number of observations: 
TC.TLs <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
TC.TLi <- ddply(S2, .(Week, TL), summarize, mean(ITL.term))
TC.ranefs <- ddply(S2, .(Week, TL), summarize, mean(sl.ranef))
TC.ranefsi <- ddply(S2, .(Week, TL), summarize, mean(Int.ranef))
TC.S <- merge(TC.TLs, TC.ranefs, by = c('Week', 'TL'), all= FALSE)
TC.I <- merge(TC.TLi, TC.ranefsi, by = c('Week', 'TL'), all= FALSE)
TCs <- merge(TC.S, TC.I, by = c('Week', 'TL'), all = FALSE)
names(TCs) <- c("Week","TL","TL.term","sl.ranef","ITL.term","Int.ranef")
TCs$Ea <- TCs$TL.term + TCs$sl.ranef
TCs$Bo <- TCs$ITL.term + TCs$Int.ranef

Ea.sumTC <- ddply(TCs, .(TL), summarize, mean(Ea))
Ea.sumTC2 <- ddply(TCs, .(TL), summarize, se(Ea))
Ea.sumTCs <- merge(Ea.sumTC, Ea.sumTC2, by = 'TL')
names(Ea.sumTCs) <- c('TL', 'means','se')

Bo.sumTC <- ddply(TCs, .(TL), summarize, mean(Bo))
Bo.sumTC2 <- ddply(TCs, .(TL), summarize, se(Bo))
Bo.sumTCs <- merge(Bo.sumTC, Bo.sumTC2, by = 'TL')
names(Bo.sumTCs) <- c('TL', 'means','se')


### merge Eas
row <- c('','','')
Ea.sum <- rbind(Ea.sumTCs, row, Ea.sumPPs, row, Ea.sumERms, row, Ea.sumNPPms,  row, Ea.sumERs, row, Ea.sumNPPs, row)

### merge Bos
row <- c('','','')
Bo.sum <- rbind(Bo.sumTCs, row, Bo.sumPPs, row, Bo.sumERms, row, Bo.sumNPPms,  row, Bo.sumERs, row, Bo.sumNPPs, row)

### FIGURE
pdf(file = "figure 2.pdf", width = 7.5, height = 6)

#SLOPES
par(mar=(c(5,4,4,1)), mfrow = c(1,2)) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(-4, 3),                        	
     ylim = c(0, length(Ea.sum[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests.1 <- as.numeric(Ea.sum[,2]) #res.sl$trt == '1TL'
ses.1 <- as.numeric(Ea.sum[,3]) #
var.names <- Ea.sum$TL

for (i in 1:length(ests.1)) {                                            
  points(ests.1[i], i, pch = 19, cex = 1.2, col = 1)
  lines(c(ests.1[i] + 1.96*ses.1[i], ests.1[i] - 1.96*ses.1[i]), c(i, i), col = 1, lwd = 2)
  text(-4.5, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)      # add the variable names
  text(3, length(Ea.sum[,1]) + .3, 'A', cex = 1.2)
  text(-4, 23, 'NPP', adj = 0, cex = 0.7)
  text(-4, 19, 'ER', adj = 0, cex = 0.7)
  text(-4, 15, 'NPP/Mb', adj = 0, cex = 0.7)
  text(-4, 11, 'ER/Mt', adj = 0, cex = 0.7)
  text(-4, 7, 'Mb', adj = 0, cex = 0.7)
  text(-4, 3, 'Mt', adj = 0, cex = 0.7)
  #text(-2.5, 21, 'Total biomass', adj = 0, cex = 0.8)
}

# add axes and labels
axis(side = 1)                                                                                         
abline(v = 0, lty = 2, col = "grey40")     
abline(v = -0.65, lty = 1, col = "grey40") 
abline(v = -0.32, lty = 1, col = "grey40") 
abline(v = 0.65, lty = 1, col = "grey40") 
abline(v = 0.32, lty = 1, col = "grey40")
abline(h = 4, lty = 3, col = 'grey40')
abline(h = 8, lty = 3, col = 'grey40')
abline(h = 12, lty = 3, col = 'grey40')
abline(h = 16, lty = 3, col = 'grey40')
abline(h = 20, lty = 3, col = 'grey40')
mtext(side = 1, "Temperature depednence (B1)", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          

#dev.off()


### FIGURE
#pdf(file = "figure 2B.pdf", width = 7, height = 5)

#INTERCEPTS
#par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(-6, 7),                        	
     ylim = c(0, length(Bo.sum[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests.1 <- as.numeric(Bo.sum[,2]) #res.sl$trt == '1TL'
ses.1 <- as.numeric(Bo.sum[,3]) #
var.names <- Bo.sum$TL

for (i in 1:length(ests.1)) {                                            
  points(ests.1[i], i, pch = 19, cex = 1.2, col = 1)
  lines(c(ests.1[i] + 1.96*ses.1[i], ests.1[i] - 1.96*ses.1[i]), c(i, i), col = 1, lwd = 2)
  text(-7, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)      # add the variable names
  text(6, length(Ea.sum[,1]) + .3, 'B', cex = 1.2)
  text(-6, 23, 'NPP', adj = 0, cex = 0.7)
  text(-6, 19, 'ER', adj = 0, cex = 0.7)
  text(-6, 15, 'NPP/Mb', adj = 0, cex = 0.7)
  text(-6, 11, 'ER/Mt', adj = 0, cex = 0.7)
  text(-6, 7, 'Mb', adj = 0, cex = 0.7)
  text(-6, 3, 'Mt', adj = 0, cex = 0.7)
  #text(-2.5, 21, 'Total biomass', adj = 0, cex = 0.8)
}

# add axes and labels
axis(side = 1)                                                        
abline(h = 4, lty = 3, col = 'grey40')
abline(h = 8, lty = 3, col = 'grey40')
abline(h = 12, lty = 3, col = 'grey40')
abline(h = 16, lty = 3, col = 'grey40')
abline(h = 20, lty = 3, col = 'grey40')                   
mtext(side = 1, "Intercepts (B0)", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          

dev.off()

