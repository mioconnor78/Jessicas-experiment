## making a figure for slopes only

## maybe standardized coefficients? review what these are, and whether activation energy can be easily inferred. 
## if we wanted this plot to be activation energies with confidence intervals, then i have to get slopes for the T x TL interactions. not sure how to do that. could maybe do it by taking the mean and ci for each week... but that's no right.

#modNPP4r
modNPP <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = TRUE, na.action=na.omit)  #modNPP4r

# create data for best model estimates
est.B <- as.data.frame(as.numeric(round(fixef(modNPP),2)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(modNPP))),2))
names(est.B) <- c('est', 'se')
rownames(est.B) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.B$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.B$resp <- c('NPP', 'NPP','NPP','NPP','NPP','NPP')

#ER
modER4r <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data[(data$ER2 >= 0),], REML = TRUE, na.action=na.omit)
modER <- modER4r

# create data for best model estimates
est.ER <- as.data.frame(as.numeric(round(fixef(modER),2)))
est.ER$se <- as.numeric(round(sqrt(diag(vcov(modER))),2))
names(est.ER) <- c('est', 'se')
rownames(est.ER) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.ER$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.ER$resp <- c('ER', 'ER','ER','ER','ER','ER')


# NPPM
modNPPm <- lmer(log(NPP.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)   

est.NPPm <- as.data.frame(as.numeric(round(fixef(modNPPm),2)))
est.NPPm$se <- as.numeric(round(sqrt(diag(vcov(modNPPm))),2))
names(est.NPPm) <- c('est', 'se')
rownames(est.NPPm) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.NPPm$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.NPPm$resp <- c('NPPm', 'NPPm','NPPm','NPPm','NPPm','NPPm')



## ERm
modERm <- lmer(log(ER.mass) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data[(data$ER2 >= 0),], REML = TRUE, na.action=na.omit) 

est.ERm <- as.data.frame(as.numeric(round(fixef(modERm),2)))
est.ERm$se <- as.numeric(round(sqrt(diag(vcov(modERm))),2))
names(est.ERm) <- c('est', 'se')
rownames(est.ERm) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.ERm$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.ERm$resp <- c('ERm', 'ERm','ERm','ERm','ERm','ERm')

## NEM
modNEM <- lmer(log(NEM + .1) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = TRUE, na.action=na.omit)

est.NEM <- as.data.frame(as.numeric(round(fixef(modNEM),2)))
est.NEM$se <- as.numeric(round(sqrt(diag(vcov(modNEM))),2))
names(est.NEM) <- c('est', 'se')
rownames(est.NEM) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.NEM$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.NEM$resp <- c('NEM', 'NEM','NEM','NEM','NEM','NEM')

#PP
modPP <- lmer(log(PP.biomass) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

est.PP <- as.data.frame(as.numeric(round(fixef(modPP),2)))
est.PP$se <- as.numeric(round(sqrt(diag(vcov(modPP))),2))
names(est.PP) <- c('est', 'se')
rownames(est.PP) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.PP$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.PP$resp <- c('P.Biomass', 'P.Biomass','P.Biomass','P.Biomass','P.Biomass','P.Biomass')

#Total Carbon
modTC <- lmer(log(total.carbon) ~ 1 + I(invT-mean(invT))*trophic.level + (I(invT-mean(invT))|week), data=data, REML = TRUE, na.action=na.omit)

est.TC <- as.data.frame(as.numeric(round(fixef(modTC),2)))
est.TC$se <- as.numeric(round(sqrt(diag(vcov(modTC))),2))
names(est.TC) <- c('est', 'se')
rownames(est.TC) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.TC$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.TC$resp <- c('Total Biomass', 'Total Biomass','Total Biomass','Total Biomass','Total Biomass','Total Biomass')

## put them all together
res <- rbind(est.B, est.ER, est.NPPm, est.ERm, est.NEM, est.PP, est.TC)
res.sl <- res[res$slint == 'S',]
res.sl$trt <- rep(c('1TL', '2TL', '3TL'))
res.int <- res[res$slint == 'I',]

### two-paneled figure
pdf(file = "figure 1B.pdf", width = 7.5, height = 4)
par(
  family = "serif",  
  oma = c(0,0,0,0),  # Since it is a single plot, I set the outer margins to zero.
  #fin = c(7,5), pty = "m",
  mar = c(5,10,4,0),  # Inner margins are set through a little trial and error.
  mfcol = c(1,2)
)

#SLOPES
par(mar=(c(5,9,4,2))) #pin = c(2.3, 3.5), 
plot(NULL,                                
     xlim = c(-2.5, 2),                        	
     ylim = c(0, length(res.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests.1 <- as.numeric(res.sl[,1]) #res.sl$trt == '1TL'
#ests.2 <- as.numeric(res.sl[res.sl$trt == '2TL',1])
#ests.3 <- as.numeric(res.sl[res.sl$trt == '3TL',1])
ses.1 <- as.numeric(res.sl[,2]) #
#ses.2 <- as.numeric(res.sl[res.sl$trt == '2TL',2])
#ses.3 <- as.numeric(res.sl[res.sl$trt == '3TL',2])
var.names <- res.sl$trt
var.namesi <- rownames(res.int)

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
mtext(side = 1, "Slope coefficients", line = 3)                                              
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

