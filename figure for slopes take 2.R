## making a figure for slopes only

## Figure X. Estimated activation energies for ecosystem carbon stocks and oxygen fluxes under three trophic conditions. Mean activation energies are estimated by taking the average of each weekly estimated activation energy (fixed effect + weekly random effect, so n = 6). Because the estimated activation energies are derived from the MEM, non-independence of treatments within weeks is considered in the random effect.

## So in practice: get the estimate (fixed + random) for each tank at each time; then take the means (and SEMs), and plot those in this figure.

############
#modNPP4r
modNPP <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)  #modNPP4r

rand.cat <- ddply(data, .(week, Tank, invT, trophic.level), summarize, mean(NPP2))
names(rand.cat) <- c('Week', 'Tank', 'invT','TL', 'NPP2')
Entry.coefs <- data.frame(coef(modNPP)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(modNPP)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## summary with appropriate number of observations: 
NPP.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
NPP.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.modNPP..week...2.))
NPP.S <- merge(NPP.TL, NPP.ranefs, by = c('Week', 'TL'), all= FALSE)
NPP.S$Ea <- NPP.S$..1.x + NPP.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumNPP <- ddply(NPP.S, .(TL), summarize, mean(Ea))
Ea.sumNPP2 <- ddply(NPP.S, .(TL), summarize, se(Ea))
Ea.sumNPPs <- merge(Ea.sumNPP, Ea.sumNPP2, by = 'TL')
names(Ea.sumNPPs) <- c('TL', 'means','se')


## repeat for ES
#################

#modERr
modER4r <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data[(data$ER2 >= 0),], REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, ID, invT, trophic.level), summarize, mean(ER2))
names(rand.cat) <- c('Week', 'Tank', 'ID', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modER4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(modER4r)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## summary with appropriate number of observations: 
ER.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
ER.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.modER4r..week...2.))
ER.S <- merge(ER.TL, ER.ranefs, by = c('Week', 'TL'), all= FALSE)
ER.S$Ea <- ER.S$..1.x + ER.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumER <- ddply(ER.S, .(TL), summarize, mean(Ea))
Ea.sumER2 <- ddply(ER.S, .(TL), summarize, se(Ea))
Ea.sumERs <- merge(Ea.sumER, Ea.sumER2, by = 'TL')


## repeat for NPPm
#################

#modNPPm4r
modNPPm4r <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, ID, invT, trophic.level), summarize, mean(NPP.mass))
names(rand.cat) <- c('Week', 'Tank', 'ID', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modNPPm4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(modNPPm4r)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## summary with appropriate number of observations: 
NPPm.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
NPPm.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.modNPPm4r..week...2.))
NPPm.S <- merge(NPPm.TL, NPPm.ranefs, by = c('Week', 'TL'), all= FALSE)
NPPm.S$Ea <- NPPm.S$..1.x + NPPm.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumNPPm <- ddply(NPPm.S, .(TL), summarize, mean(Ea))
Ea.sumNPPm2 <- ddply(NPPm.S, .(TL), summarize, se(Ea))
Ea.sumNPPms <- merge(Ea.sumNPPm, Ea.sumNPPm2, by = 'TL')


## repeat for ERm
#################

#modERm4r
modERm4r <- lmer(log(ER.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data[(data$ER2 >= 0),], REML = TRUE, na.action=na.omit) 

rand.cat <- ddply(data, .(week, Tank, ID, invT, trophic.level), summarize, mean(ER.mass))
names(rand.cat) <- c('Week', 'Tank', 'ID', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modERm4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(modERm4r)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## summary with appropriate number of observations: 
ERm.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
ERm.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.modERm4r..week...2.))
ERm.S <- merge(ERm.TL, ERm.ranefs, by = c('Week', 'TL'), all= FALSE)
ERm.S$Ea <- ERm.S$..1.x + ERm.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumERm <- ddply(ERm.S, .(TL), summarize, mean(Ea))
Ea.sumERm2 <- ddply(ERm.S, .(TL), summarize, se(Ea))
Ea.sumERms <- merge(Ea.sumERm, Ea.sumERm2, by = 'TL')


#modBBr
modPP4r <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, ID, invT, trophic.level), summarize, mean(PP.biomass))
names(rand.cat) <- c('Week', 'Tank', 'ID', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modPP4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(modPP4r)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## summary with appropriate number of observations: 
PBm.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
PBm.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.modPP4r..week...2.))
PBm.S <- merge(PBm.TL, PBm.ranefs, by = c('Week', 'TL'), all= FALSE)
PBm.S$Ea <- PBm.S$..1.x + PBm.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumPBm <- ddply(PBm.S, .(TL), summarize, mean(Ea))
Ea.sumPBm2 <- ddply(PBm.S, .(TL), summarize, se(Ea))
Ea.sumPBms <- merge(Ea.sumPBm, Ea.sumPBm2, by = 'TL')

###################
#modTCr
modTC4r <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

rand.cat <- ddply(data, .(week, Tank, ID, invT, trophic.level), summarize, mean(total.carbon))
names(rand.cat) <- c('Week', 'Tank', 'ID', 'invT','TL', 'ER2')
Entry.coefs <- data.frame(coef(modTC4r)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(modTC4r)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## summary with appropriate number of observations: 
TC.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
TC.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.modTC4r..week...2.))
TC.S <- merge(TC.TL, TC.ranefs, by = c('Week', 'TL'), all= FALSE)
TC.S$Ea <- TC.S$..1.x + TC.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sumTC <- ddply(TC.S, .(TL), summarize, mean(Ea))
Ea.sumTC2 <- ddply(TC.S, .(TL), summarize, se(Ea))
Ea.sumTCs <- merge(Ea.sumTC, Ea.sumTC2, by = 'TL')




### merge Eas
row <- c('','','')
Ea.sum <- rbind(Ea.sumTCs, row, Ea.sumPBms, row, Ea.sumERms, row, Ea.sumNPPms,  row, Ea.sumERs, row, Ea.sumNPPs, row)


### FIGURE
pdf(file = "figure 1B.pdf", width = 7, height = 5)

#SLOPES
par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
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
  text(3, length(Ea.sum[,1]) + .3, 'B', cex = 1.2)
  text(-4, 23, 'NPP', adj = 0, cex = 0.8)
  text(-4, 19, 'ER', adj = 0, cex = 0.8)
  text(-4, 15, 'NPP/Mb', adj = 0, cex = 0.8)
  text(-4, 11, 'ER/Mt', adj = 0, cex = 0.8)
  text(-4, 7, 'Phytoplankton biomass', adj = 0, cex = 0.8)
  text(-4, 3, 'Total biomass', adj = 0, cex = 0.8)
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
mtext(side = 1, "Activation energies (Ea)", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          

dev.off()



