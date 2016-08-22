#### trying some new plots for dealing with the time effect

# run analysis

library(car)
library(rgl)

plot(log(data$calc.NPP)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)

scatter3d(data$invT, log(data$calc.NPP), data$week, point.col = as.numeric(data$trophic.level), surface = FALSE)


scatter3d(data$week, log(data$NPP2), data$invT, groups = data$trophic.level, point.col = as.numeric(data$trophic.level), surface.col = c("green", "red","black"), axis.col = c("black","black","black"), surface = FALSE, parallel = FALSE, ylab = 'ln(NPP.daily)', xlab='Week', zlab = 'T (1/kT)')
scatter3d(data$invT, log(data$NPP2), data$week, groups = data$trophic.level, surface.col = c("green", "red","black"), axis.col = c("black","black","black"), surface = TRUE, parallel = FALSE, ylab = 'ln(NPP)', zlab='Week', xlab = 'T (1/kT)') 

scatter3d(data$week, log(data$ER2), data$invT, point.col = as.numeric(data$trophic.level), surface = FALSE, ylab = 'ln(ER.daily)', xlab='Week', zlab = 'T (1/kT)')
scatter3d(data$invT, log(data$ER2), data$week, groups = data$trophic.level, surface.col = c("green", "red","black"), axis.col = c("black","black","black"), surface = TRUE, parallel = FALSE, ylab = 'ln(ER)', zlab='Week', xlab = 'T (1/kT)') 


scatter3d(data$week, data$NEM, data$invT, point.col = as.numeric(data$trophic.level), surface = FALSE, ylab = 'NEM', xlab='Week', zlab = 'T (1/kT)')
scatter3d(data$invT, data$NEM, data$week, groups = data$trophic.level, surface.col = c("green", "red","black"), axis.col = c("black","black","black"), surface = TRUE, parallel = FALSE, ylab = 'NEM', zlab='Week', xlab = 'T (1/kT)') 

scatter3d(data$week, log(data$PP.biomass), I(data$invT-mean(data$invT)), point.col = as.numeric(data$trophic.level), surface = TRUE)
scatter3d(I(data$invT-mean(data$invT)), log(data$PP.biomass), data$week, groups = data$trophic.level, surface.col = c("green", "red","black"), axis.col = c("black","black","black"), surface = TRUE, parallel = FALSE, ylab = 'ln(Bp)', zlab='Week', xlab = 'T (1/kT)') 

scatter3d(data$invT, log(data$total.carbon), data$week, groups = data$trophic.level, surface.col = c("green", "red","black"), axis.col = c("black","black","black"), surface = TRUE, parallel = FALSE, ylab = 'ln(Bt)', zlab='Week', xlab = 'T (1/kT)') 


scatter3d(dataz$invT, log(dataz$zoo.ug.carbon.liter + 0.01), dataz$week, groups = dataz$trophic.level, surface.col = c("red","black", "green"), axis.col = c("black","black","black"), surface = TRUE, parallel = FALSE, ylab = 'ln(Bz)', zlab='Week', xlab = 'T (1/kT)') 

#not working for some reason
scatter3d(dataz$week, log(dataz$zoo.ug.carbon.liter), dataz$invT, point.col = as.numeric(dataz$trophic.level), surface = TRUE)


savePlot(file = "JG3dNPP.pdf", type = "jpeg")


jpeg(file = "JG3dNPP.jpg", width = 5, height = 5)
scatter3d(data$week, log(data$calc.NPP), data$invT, point.col = as.numeric(data$trophic.level), surface = FALSE)
dev.off()
