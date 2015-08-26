#### trying some new plots for dealing with the time effect

# run analysis

library(car)
library(rgl)

plot(log(data$calc.NPP)~data$invT, pch = as.numeric(data$trophic.level)+14, col = data$trophic.level)

scatter3d(data$invT, log(data$calc.NPP), dadta$week)
