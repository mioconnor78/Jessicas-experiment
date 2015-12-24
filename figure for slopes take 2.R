## making a figure for slopes only

## Figure X. Estimated activation energies for ecosystem carbon stocks and oxygen fluxes under three trophic conditions. Mean activation energies are estimated by taking the average of each weekly estimated activation energy (fixed effect + weekly random effect, so n = 6). Because the estimated activation energies are derived from the MEM, non-independence of treatments within weeks is considered in the random effect.

## So in practice: get the estimate (fixed + random) for each tank at each time; then take the means (and SEMs), and plot those in this figure.

# need an indentifier for each week x TL
data$ID <- paste(data$week, data$trophic.level, sep = '.')

#modNPP4r
modNPP <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1 + I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)  #modNPP4r

# get fixed effect coefficient estimates for each tank
mod <- modNPP

rand.cat <- ddply(data, .(week, Tank, ID, invT, trophic.level), summarize, mean(NPP2))
names(rand.cat) <- c('Week', 'Tank', 'ID', 'invT','TL', 'NPP2')
Entry.coefs <- data.frame(coef(mod)$week) #
Entry.coefs$week <- rownames(Entry.coefs)
S <- merge(rand.cat, Entry.coefs, by.x = 'Week', by.y = 'week', all = FALSE)

## construct predicted slopes for trophic groups
S$TL.term <- ifelse(S$TL == 'P', S$I.invT...mean.invT.., 0)
S$TL.term <- ifelse(S$TL == 'PZ', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZ), S$TL.term)
S$TL.term <- ifelse(S$TL == 'PZN', (S$I.invT...mean.invT.. + S$I.invT...mean.invT...trophic.levelPZN), S$TL.term)

#add week ranefs
ranefs <- data.frame(ranef(mod)$week[,2])
ranefs$week <- as.numeric(rownames(ranefs))+3
S2 <- merge(S, ranefs, by.x = 'Week', by.y = 'week', all= FALSE)

## still need this?
b <- as.numeric(fixef(mod)[2])

## estimating weekly slopes for models with no random slopes
S2$slope <- (S2$TL.term + S2$ranef.mod..week...2.)
S2$slope.se <- st.dev(ranefs[1])/sqrt(length(ranefs[1]))

summary(S2$slope)

## summary with appropriate number of observations: 
NPP.TL <- ddply(S2, .(Week, TL), summarize, mean(TL.term))
NPP.ranefs <- ddply(S2, .(Week, TL), summarize, mean(ranef.mod..week...2.))
NPP.S <- merge(NPP.TL, NPP.ranefs, by = c('Week', 'TL'), all= FALSE)
NPP.S$Ea <- NPP.S$..1.x + NPP.S$..1.y

se <- function(x) sd(x)/sqrt(length(x))

Ea.sum <- ddply(NPP.S, .(TL), summarize, mean(Ea))
Ea.sum2 <- ddply(NPP.S, .(TL), summarize, se(Ea))
Ea.sums <- merge(Ea.sum, Ea.sum2, by = 'TL')


### FIGURE
pdf(file = "figure 1B.pdf", width = 7.5, height = 4)

#SLOPES
par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(-2.5, 2),                        	
     ylim = c(0, length(res.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests.1 <- as.numeric(Ea.sums[,2]) #res.sl$trt == '1TL'
#ests.2 <- as.numeric(res.sl[res.sl$trt == '2TL',1])
#ests.3 <- as.numeric(res.sl[res.sl$trt == '3TL',1])
ses.1 <- as.numeric(Ea.sums[,3]) #
#ses.2 <- as.numeric(res.sl[res.sl$trt == '2TL',2])
#ses.3 <- as.numeric(res.sl[res.sl$trt == '3TL',2])
var.names <- Ea.sums$TL
#var.namesi <- rownames(res.int)

for (i in 1:length(ests.1)) {                                            
  points(ests.1[i], i, pch = 19, cex = 1.2, col = 1)
  #points(ests.2[i], i, pch = 19, cex = 1.2, col = 'gray80')
  #points(ests.3[i], i, pch = 19, cex = 1.2, col = 'gray60')
  lines(c(ests.1[i] + 1.96*ses.1[i], ests.1[i] - 1.96*ses.1[i]), c(i, i), col = 1, lwd = 2)
  #lines(c(ests.2[i] + 1.96*ses.2[i], ests.2[i] - 1.96*ses.2[i]), c(i, i), col = 1, lwd = 2)
  #lines(c(ests.3[i] + 1.96*ses.3[i], ests.3[i] - 1.96*ses.3[i]), c(i, i), col = 1, lwd = 2)
  text(-2.8, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)      # add the variable names
  text(2, length(res.sl[,1]) + .3, 'B', cex = 1.2)
  text(-2.5, 3, 'NPP', adj = 0, cex = 0.8)
  text(-2.5, 6, 'ER', adj = 0, cex = 0.8)
  text(-2.5, 9, 'NPP/Mb', adj = 0, cex = 0.8)
  text(-2.5, 12, 'ER/Mt', adj = 0, cex = 0.8)
  text(-2.5, 15, 'NEM', adj = 0, cex = 0.8)
  text(-2.5, 18, 'Phytoplankton biomass', adj = 0, cex = 0.8)
  text(-2.5, 21, 'Total biomass', adj = 0, cex = 0.8)
}

# add axes and labels
axis(side = 1)                                                                                         
abline(v = 0, lty = 2, col = "grey40")     
abline(v = -0.65, lty = 1, col = "grey40") 
abline(v = -0.32, lty = 1, col = "grey40") 
abline(v = 0.65, lty = 1, col = "grey40") 
abline(v = 0.32, lty = 1, col = "grey40")
abline(h = 3.5, lty = 3, col = 'grey40')
abline(h = 6.5, lty = 3, col = 'grey40')
abline(h = 9.5, lty = 3, col = 'grey40')
abline(h = 12.5, lty = 3, col = 'grey40')
abline(h = 15.5, lty = 3, col = 'grey40')
abline(h = 18.5, lty = 3, col = 'grey40')
mtext(side = 1, "Activation energies (Ea)", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                                          





### INTERCEPTS
par(mar=c(5,8,4,4))  #pin = c(2.3, 3.5)), but this doesn't seem to work with mar
plot(NULL,                                
     xlim = c(-6, 6),                          
     ylim = c(.7, length(est.B.int[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA)

# add the data
#est <- as.numeric(est.int[,1]) 
#se <- as.numeric(est.int[,2] )                                         
ests.B <- as.numeric(est.B.int[,1])
ses.B <- as.numeric(est.B.int[,2])
#ests.Ba <- as.numeric(est.Ba.int[,1])
#ses.Ba <- as.numeric(est.Ba.int[,2])
var.names<-rownames(est.B.int)

b <- 0
for (i in 1:length(ests.B)) {                                            
  #points(est[i], i, pch = 19, cex = 1.2)                              
  #lines(c(est[i] + 1.96*se[i], est[i] - 1.96*se[i]), c(i, i), lwd = 2)  
  points(ests.B[i], i+b, pch = 19, cex = 1.2, col = 1) 
  lines(c(ests.B[i] + 1.96*ses.B[i], ests.B[i] - 1.96*ses.B[i]), c(i+b, i+b), col = 1, lwd = 2)
  #lines(c(ests.Ba[i] + 1.96*ses.Ba[i], ests.Ba[i] - 1.96*ses.Ba[i]), c(i+2*b, i+2*b), col = 'gray50', lwd = 2)
  #points(ests.Ba[i], i+2*b, pch = 19, cex = 1.2, col = 'gray50')   # add 95% CIs
  text(-7, i, adj = c(1,0), var.namesi[i], xpd = T, cex = .8)        # add the variable names
  text(5.5, length(est.B.int[,1])+ 0.2, 'C', cex = 1.2)
}

# add axes and labels
axis(side = 1, at = c(-6, 0, 6))
#axis(side = 2, pos = -2)
abline(v = 0, lty = 3, col = "grey40")                                                                   
mtext(side = 1, "Intercept coefficients", line = 3)                                              
mtext(side = 3, "", line = 1, cex = 0.8)   # add title
box()                    

dev.off()

