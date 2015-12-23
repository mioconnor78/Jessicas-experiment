## making a figure for slopes only

#modNPP4r
modNPP <- lmer(log(NPP2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data, REML = TRUE, na.action=na.omit)  #modNPP4r

# create data for best model estimates
est.B <- as.data.frame(as.numeric(round(fixef(modNPP),2)))
est.B$se <- as.numeric(round(sqrt(diag(vcov(modNPP))),2))
names(est.B) <- c('est', 'se')
rownames(est.B) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.B$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.B$resp <- c('NPP', 'NPP','NPP','NPP','NPP','NPP')
est.B.sl <- est.B[est.B$slint == 'S',]
est.B.int <- est.B[est.B$slint == 'I',]

#ER
modER4r <- lmer(log(ER2) ~ 1 + I(invT-mean(invT))*trophic.level + (1|week), data=data1, REML = TRUE, na.action=na.omit)
modER <- modER4r

# create data for best model estimates
est.ER <- as.data.frame(as.numeric(round(fixef(modER),2)))
est.ER$se <- as.numeric(round(sqrt(diag(vcov(modER))),2))
names(est.ER) <- c('est', 'se')
rownames(est.ER) <- c('Intercept', 'ln(T)', '2TL', '3TL', 'T*2TL', 'T*3TL')
est.ER$slint <- c('I',  'S', 'I', 'I', 'S', 'S')
est.ER$resp <- c('ER', 'ER','ER','ER','ER','ER')
est.ER.sl <- est.ER[est.ER$slint == 'S',]
est.ER.int <- est.ER[est.ER$slint == 'I',]

res <- rbind(est.ER, est.B)
res.sl <- res[res$slint == 'S',]
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
     xlim = c(-1.0, 1),                        	
     ylim = c(0, length(res.sl[,1]) + .3), 	
     axes = F, xlab = NA, ylab = NA, cex = 0.8)

# add the data
ests <- as.numeric(res.sl[,1])
ses <- as.numeric(res.sl[,2])
var.names<-rownames(res.sl)
var.namesi<-rownames(res.int)

for (i in 1:length(ests)) {                                            
  points(ests[i], i, pch = 19, cex = 1.2, col = 1) 
  lines(c(ests[i] + 1.96*ses[i], ests[i] - 1.96*ses[i]), c(i, i), col = 1, lwd = 2)
  #points(ests.Ba[i], i+2*b, pch = 19, cex = 1.2, col = 'gray50') 
  #lines(c(ests.Ba[i] + 1.96*ses.Ba[i], ests.Ba[i] - 1.96*ses.Ba[i]), c(i+2*b, i+2*b), col = 'gray50', lwd = 2)
  text(-1.2, i, adj = c(1,0), var.names[i], xpd = T, cex = .8)        # add the variable names
  text(1, length(res.sl[,1]) + .3, 'B', cex = 1.2)
}

# add axes and labels
axis(side = 1)                                                                                         
abline(v = 0, lty = 3, col = "grey40")     
abline(v = -0.65, lty = 3, col = "grey40") 
abline(v = -0.32, lty = 3, col = "grey40") 
abline(h = 3.5, lty = 3, col = 'grey40')
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

