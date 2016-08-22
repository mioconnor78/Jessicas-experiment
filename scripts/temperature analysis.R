#################################################
#### Garzke, O'Connor and Sommer temperature experiment
#### File for estimating average temperatures. Necessary for mixed effects models.
#### Code by Mary O'Connor, revised November 2015
#################################################


### estimating weekly mean temperatures

library(plyr)

### for each tank, each week, we will want: 

### weekly average of daily average temps

tdata <- read.csv("./data/avgtemps.csv")
head(tdata)
temp.data <- ddply(tdata, .(Week, Tank), summarise, mean(Temperature)) 
head(temp.data)
names(temp.data) <- c('Week', 'Tank', 'wklyTemp')

tank.means <- ddply(tdata, .(Tank), summarise, mean(Temperature)) 
head(tank.means)
names(tank.means) <- c("Tank", "TankTemp")

### this is the old way, don't need this now.
## bring in data file with temperatures at each sampling time
o.data <- read.csv("./oxygen_temp_temporal.csv")
o.data <- o.data[,-(4:14)]
o.data <- o.data[,-2]
head(o.data)
dim(o.data)
data2 <- merge(temp.data, o.data, by.x = c("Week", "Tank"), by.y = c("week", "Tank"))
head(data2)

data2$meanOtemp <- (data2$temp.dawn1 + data2$temp.dusk + data2$temp.dawn2)/3

head(data2)
plot(data2$..1 ~ data2$meanOtemp, xlim = c(15,35), ylim = c(15, 35), xlab = 'mean dawn-dusk-dawn', ylab = 'mean daily from dataloggers')
abline(0, 1)

plot(data2$..1 ~ data2$temp.dusk, xlim = c(15,35), ylim = c(15, 35), xlab = 'mean dawn-dusk-dawn', ylab = 'mean daily from dataloggers')
abline(0, 1)


# old code from analysis file to use these temps
## bring in data file with temperatures at each sampling time
o.data <- read.csv("./oxygen_temp_temporal.csv")
o.data <- o.data[,-(4:14)]
o.data <- o.data[,-2]
head(o.data)
dim(o.data)
data2 <- merge(data, o.data, by.x = c("week", "Tank"), by.y = c("week", "Tank"))
data <- data2

data3 <- merge(data, mass.data, by.x = c("week", "Tank"), by.y = c("week", "tank"), all=TRUE)
data3[is.na(data3)] <- 0
data <- data3

