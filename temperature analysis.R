#################################################
#### Garzke, O'Connor and Sommer temperature experiment
#### File for estimating average temperatures. Necessary for mixed effects models.
#### Code by Mary O'Connor, revised November 2015
#################################################


### estimating weekly mean temperatures

library(plyr)

### for each tank, each week, we will want: 

### weekly average of daily average temps

tdata <- read.csv("./avgtemps.csv")
head(tdata)

temp.data <- ddply(tdata, .(Week, Tank), summarise, mean(Temperature)) 

head(temp.data)
