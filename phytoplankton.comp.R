setwd("/Users/maryo/Documents/temporary files/Jessicas experment")
pp <- read.csv("PPweek7.csv")
dim(pp)
names(pp)


plot(pp$Cyano.Diatom.Ratio~pp$average.temp, cex=1.5, pch='',  axes=FALSE, ylim=c(0,0.5), xlim=c(14,28), xlab='inv Temperature (C)', ylab='Diatom/cyano density') 
axis(1, at=c(14,18,22,24,28), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,.25,.5), pos=14, lwd=2, cex.lab=1.5)
points(pp[(pp$trophic.level=='P'),]$Cyano.Diatom.Ratio~pp[(pp$trophic.level=='P'),]$average.temp, pch=19, col = 'seagreen', cex = 1.5)
points(pp[(pp$trophic.level=='PZ'),]$Cyano.Diatom.Ratio~pp[(pp$trophic.level=='PZ'),]$average.temp, pch=15, col = 'brown', cex = 1.5)
points(pp[(pp$trophic.level=='PZN'),]$Cyano.Diatom.Ratio~pp[(pp$trophic.level=='PZN'),]$average.temp, pch=17, col = 'blue')
legend(26, 0.5, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))


plot(pp$proportion.cyano~pp$average.temp, cex=1.5, pch='',  axes=FALSE, ylim=c(0,1), xlim=c(14,28), xlab='inv Temperature (C)', ylab='Diatom/cyano density') 
axis(1, at=c(14,16,18,20,22,24,26,28), pos=0, lwd=2, cex.lab=1.5)
axis(2, at=c(0,.25,.5,0.75,1), pos=14, lwd=2, cex.lab=1.5)
points(pp[(pp$trophic.level=='P'),]$proportion.cyano~pp[(pp$trophic.level=='P'),]$average.temp, pch=19, col = 'seagreen', cex = 1.5)
points(pp[(pp$trophic.level=='PZ'),]$proportion.cyano~pp[(pp$trophic.level=='PZ'),]$average.temp, pch=15, col = 'brown', cex = 1.5)
points(pp[(pp$trophic.level=='PZN'),]$proportion.cyano~pp[(pp$trophic.level=='PZN'),]$average.temp, pch=17, col = 'blue')
legend(26, 0.5, c('1 TL', '2 TL','3 TL'), pch = c(19, 15, 17), col = c('seagreen', 'brown', 'blue'))
