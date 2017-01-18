### Mary's community composition analysis

data4 <- read.csv("./phytoplankton.speciesMO.csv")
View(data4)

library(vegan)


count.func <- function(x) {
  y = ifelse( x > 0, 1, 0)
  sum(y)
}

k <- 8.617342*10^-5  # eV/K
data4$invT <-  1/((data4$average.temp + 273)*k)
data4$invT <- as.numeric(as.character(data4$invT))

data1 <- merge(data.sum2, sites, by.x = 'site', by.y = 'site')

data4$H <- diversity(data4[,7:34], index = "shannon") 
data4$S <- diversity(data4[,7:34], index = "simpson")
data4$N <- apply(data4[,7:34], 1, function(x) sum(x))
data4$R <- apply(data4[,7:34], 1, function(x) count.func(x))
head(data4)

plot(data4$H ~ data4$invT)
plot(data4$S ~ data4$invT)
plot(log(data4$N) ~ data4$invT)
plot(data4$R ~ data4$invT)


mod1<-lm(data4$R ~ data4$invT)
summary(mod1)
mod2<-lm(log(data4$N+0.01) ~ data4$invT)
summary(mod2)
mod3<-lm(data4$H ~ data4$invT)
summary(mod3)
mod4<-lm(data4$S ~ data4$invT)
summary(mod4)
