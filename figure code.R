#### trying some new plots for dealing with the time effect

# run analysis

library(car)
library(rgl)

plot(log(data$calc.NPP)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)

scatter3d(data$invT, log(data$calc.NPP), data$week, point.col = as.numeric(data$trophic.level), surface = FALSE)


scatter3d(data$week, log(data$calc.NPP), data$invT, point.col = as.numeric(data$trophic.level), surface = TRUE)

scatter3d(data$week, log(data$ER), data$invT, point.col = as.numeric(data$trophic.level), surface = FALSE, ylab = 'ln(ER)', xlab='Week', zlab = 'T (1/kT)')

scatter3d(data$week, log(data$NEM), data$invT, point.col = as.numeric(data$trophic.level), surface = FALSE, ylab = 'ln(NEM)', xlab='Week', zlab = 'T (1/kT)')

scatter3d(data$invT, log(data$NEM), data$week, point.col = as.numeric(data$trophic.level), surface = FALSE, ylab = 'ln(NEM)', xlab='T (1/kT)', zlab = 'week')

scatter3d(data$week, log(data$PP.biomass), data$invT, point.col = as.numeric(data$trophic.level), surface = TRUE)
scatter3d(data$invT, log(data$PP.biomass), data$week, point.col = as.numeric(data$trophic.level), surface = FALSE)

#not working for some reason
scatter3d(dataz$week, log(dataz$zoo.ug.carbon.liter), dataz$invT, point.col = as.numeric(dataz$trophic.level), surface = TRUE)


savePlot(file = "JG3dNPP.pdf", type = "jpeg")


jpeg(file = "JG3dNPP.jpg", width = 5, height = 5)
scatter3d(data$week, log(data$calc.NPP), data$invT, point.col = as.numeric(data$trophic.level), surface = FALSE)
dev.off()
