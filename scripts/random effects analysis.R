#### Garzke et al SOM figure on random effects

mod <- modNPPm4r
ranefs <- data.frame(ranef(mod)$week)

pdf(file = "figure S2.pdf", width = 7.5, height = 6)
par(mfrow=c(2,2), omi=c(.1,.1,.1,.1))
plot(rownames(ranefs), ranefs[,1], pch = 19, cex = 2, xlab = 'Week', ylab = 'random effect mu.0j', ylim=c(-2, 2), main = 'Intercept, NPP.Mp')
plot(rownames(ranefs), ranefs[,2], pch = 19, cex = 2, xlab = 'Week', ylab = 'random effect mu.1j', ylim=c(-1, 1), main = 'Slope, NPP.Mp')
m1 <- lm(ranefs[,1]~as.numeric(rownames(ranefs)))
summary(m1)
m2 <- lm(ranefs[,2]~as.numeric(rownames(ranefs)))
summary(m2)


mod <- modPP4r
ranefs <- data.frame(ranef(mod)$week)
plot(rownames(ranefs), ranefs[,1], pch = 19, cex = 2, xlab = 'Week', ylab = 'random effect mu.0j', ylim=c(-2, 2), main = 'Intercept, Mp')
plot(rownames(ranefs), ranefs[,2], pch = 19, cex = 2, xlab = 'Week', ylab = 'random effect mu.1j', ylim=c(-0.6, 0.6), main = 'Slope, Mp')
m1 <- lm(ranefs[,1]~as.numeric(rownames(ranefs)))
summary(m1)

m2 <- lm(ranefs[,2]~as.numeric(rownames(ranefs)))
summary(m2)


dev.off()
