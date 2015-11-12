<<<<<<< Updated upstream
data4 <- read.csv("./phytoplankton.species.csv")
View(data4)
attach(data4)

library(vegan)
=======
data4 <- read.csv("./Cyano.abundance.csv")
View(data4)
attach(data4)

>>>>>>> Stashed changes
library(nlme)
library(MuMIn)
library(lme4)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)

k <- 8.617342*10^-5  # eV/K
data4$invT <-  1/((data4$average.temp + 273)*k)
data4$invT <- as.numeric(as.character(data4$invT))

data4$Cyano.Rest.ratio <- (data4$Cyano.abundance/data4$Rest)
hist(data4$Cyano.Rest)
hist(log(data4$Cyano.Rest))

<<<<<<< Updated upstream
=======
hist(data4$Cyano.SumPhyto.ratio)
#hist(log(data4$Cyano.SumPhyto.ratio))


>>>>>>> Stashed changes
# basic plot
CyanoRest0 <- ggplot(data4, aes(x=invT, y=log(Cyano.Rest.ratio), colour= factor(trophic.level))) + 
  #geom_line()+
  geom_point(size = 3.5) +
  xlab('Temperature 1/k(T)') +
  ylab(expression('Cyanobacteria / Phytoplankton ratio')) +
  scale_y_continuous(limits= c(-4.5, 3.5), breaks=c(-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0))+
  scale_x_continuous(limits= c (38.0, 41.0), breaks=c(38.0, 39.0, 40.0, 41.0))+
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(-4.5, 3.5))+
  geom_abline(aes(intercept = 26.339954, slope = -0.689678), color='#000000', lty=1, size= 1, data=data4)+
  geom_text(aes(label = 'D', x= 40.8, y=-Inf), vjust= -0.5, size=10, colour='black')+
  theme_bw()+
  scale_colour_manual(values=c('#000000', '#666666', '#CCCCCC'), breaks=c('P', 'PZ', 'PZN'), labels=c('1-TL', '2-TL', '3-TL'))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face='bold', colour='black', size=16, angle=0, vjust=0),
        axis.text.x = element_text(face='bold', colour='black', size=16, angle=0, vjust=-1),
        axis.title.x = element_text(colour='black', size=18),
        axis.title.y = element_text(colour='black', size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=16, face='bold'),
        legend.key = element_rect(colour = 'white', fill = 'white')
  )
CyanoRest0

# plot with trasformed axis
CyanoRest1 <- ggplot(data4, aes(x=invT, y=log(Cyano.Rest.ratio), colour= factor(trophic.level))) + 
  #geom_line()+
  geom_point(size = 3.5) +
  xlab('Temperature C') +
  ylab(expression('Cyanobacteria / Phytoplankton ratio')) +
  scale_y_continuous(limits= c(-4.5, 3.5), breaks=c(-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0))+
  scale_x_continuous(limits=c(38, 41),breaks=c(38, 39, 40, 41),labels=c('30.0', '24.5', '16.2', '12.8')) +
  coord_cartesian(xlim=c(38, 41), ylim=c(-4.5, 3.5))+
  geom_abline(aes(intercept = 26.339954, slope = -0.689678), color='#000000', lty=1, size= 1, data=data4)+
  geom_text(aes(label = 'D', x= 40.8, y=-Inf), vjust= -0.5, size=10, colour='black')+
  theme_bw()+
  scale_colour_manual(values=c('#000000', '#666666', '#CCCCCC'), breaks=c('P', 'PZ', 'PZN'), labels=c('1-TL', '2-TL', '3-TL'))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face='bold', colour='black', size=16, angle=0, vjust=0),
        axis.text.x = element_text(face='bold', colour='black', size=16, angle=0, vjust=-1),
        axis.title.x = element_text(colour='black', size=18),
        axis.title.y = element_text(colour='black', size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=16, face='bold'),
        legend.key = element_rect(colour = 'white', fill = 'white')
  ) 
CyanoRest1

## extract gtable
A1 <- ggplot_gtable(ggplot_build(CyanoRest0)) #with invT
A2 <- ggplot_gtable(ggplot_build(CyanoRest1)) # with Celsius

## overlap the panel of the 2nd plot on that of the 1st plot
pp <- c(subset(A1$layout, name=="panel", se=t:r))

A <- gtable_add_grob(A1, A1$grobs[[which(A1$layout$name=="panel")]], pp$t, pp$l, pp$b, pp$l)

## steal axis from second plot (with Ceslius) and modify
iA <- which(A2$layout$name == "axis-b")
gA <- A2$grobs[[iA]]
Ax <- gA$children[[2]]

## switch position of ticks and labels
Ax$heights <- rev(Ax$heights)
Ax$grobs <- rev(Ax$grobs)
Ax$grobs[[2]]$y <- Ax$grobs[[2]]$y - unit(1, "npc") + unit(0.15, "cm")

## modify existing row to be tall enough for axis
A$heights[[2]] <- A$heights[A2$layout[iA,]$t]

## add new axis
A <- gtable_add_grob(A, Ax, 2, 4, 2, 4)

## add new row for upper axis label
A <- gtable_add_rows(A, A2$heights[1], 1)
A <- gtable_add_grob(A, A2$grob[[6]], 2, 4, 2, 4)

# draw it
grid.draw(A)

###### Analysis
modCyano0<-lme(log(Cyano.Rest.ratio) ~ 1, random=~1|week, method="ML", data=data4)
modCyano1<-lme(log(Cyano.Rest.ratio) ~ 1+I(invT-mean(invT)), random=~1|week, method="ML", data=data4)
modCyano2<-lme(log(Cyano.Rest.ratio) ~ 1+I(invT-mean(invT)) + trophic.level,  random=~1|week, method="ML", data=data4)
modCyano4<-lme(log(Cyano.Rest.ratio) ~ 1+I(invT-mean(invT)) * trophic.level, random=~1|week, method="ML", data=data4)

model.sel(modCyano0, modCyano1, modCyano2, modCyano4)
#Model selection table 
#(Int)   inT-men(inT) trp.lvl I(inT-men(inT)):trp.lvl df
#2 -0.5526 -0.7424                                      4 
#3 -0.4890 -0.7410      +                               6 
#1 -0.5526                                              3 
#4 -0.4889 -0.7905      +       +                       8 
#logLik   AICc  delta weight
#2 -190.857 390.1 0.00  0.819 
#3 -190.606 394.0 3.89  0.117 
#1 -195.145 396.5 6.43  0.033 
#4 -189.651 396.6 6.54  0.031 
#Random terms (all models): 
#  ‘1 | week’

anova(modCyano1, modCyano2)
#          Model df      AIC      BIC    logLik   Test
#modCyano1     1  4 389.7144 400.8644 -190.8572       
#modCyano2     2  6 393.2119 409.9368 -190.6060 1 vs 2
#L.Ratio p-value
#modCyano1                  
#modCyano2 0.5025182  0.7778

anova(modCyano2, modCyano0)
#          Model df      AIC      BIC    logLik   Test
#modCyano1     1  4 389.7144 400.8644 -190.8572       
#modCyano2     2  6 393.2119 409.9368 -190.6060 1 vs 2
#L.Ratio p-value
#modCyano1                  
#modCyano2 0.5025182  0.7778

anova(modCyano0, modCyano4)
#          Model df      AIC      BIC    logLik   Test
#modCyano0     1  3 396.2892 404.6517 -195.1446       
#modCyano4     2  8 395.3011 417.6010 -189.6505 1 vs 2
#L.Ratio p-value
#modCyano0                 
#modCyano4 10.98813  0.0516

summary(modCyano1)
#Linear mixed-effects model fit by maximum likelihood
#Data: data4 
#AIC      BIC    logLik
#389.7144 400.8644 -190.8572

#Random effects:
#  Formula: ~1 | week
#(Intercept) Residual
#StdDev:    0.339839 1.162244

#Fixed effects: log(Cyano.Rest.ratio) ~ 1 + I(invT - mean(invT)) 
#Value Std.Error  DF   t-value
#(Intercept)          -0.5525865 0.2020137 115 -2.735391
#I(invT - mean(invT)) -0.7424002 0.2204106 115 -3.368260
#p-value
#(Intercept)           0.0072
#I(invT - mean(invT))  0.0010
#Correlation: 
#  (Intr)
#I(invT - mean(invT)) 0     

#Standardized Within-Group Residuals:
#Min         Q1        Med         Q3        Max 
#-2.6106153 -0.6093780  0.0221220  0.5821335  2.6944966 

#Number of Observations: 120
#Number of Groups: 4 

# Model refitting
modCyano1<-lme(log(Cyano.Rest.ratio) ~ 1+invT, random=~1|week, method="REML", data=data4)
summary(modCyano1)
#Linear mixed-effects model fit by REML
#Data: data4 
#AIC      BIC   logLik
#392.054 403.1367 -192.027

#Random effects:
#  Formula: ~1 | week
#(Intercept) Residual
#StdDev:   0.4453598  1.16517

#Fixed effects: log(Cyano.Rest.ratio) ~ 1 + invT 
#Value Std.Error  DF   t-value p-value
#(Intercept) 26.339954  8.857701 115  2.973678  0.0036
#invT        -0.689678  0.227074 115 -3.037242  0.0030
#Correlation: 
#  (Intr)
#invT -1    

#Standardized Within-Group Residuals:
#  Min           Q1          Med           Q3 
#-2.582455325 -0.602564282  0.005854868  0.594342878 
#Max 
#2.691480593 

#Number of Observations: 120
<<<<<<< Updated upstream
#Number of Groups: 4 
=======
#Number of Groups: 4 


#### MARY ADDED THIS BIT ON OCT 20

#### analysis for cyanos as proportion of total - all weeks
modCprop0<-lme((Cyano.SumPhyto.ratio) ~ 1, random=~1|week, method="ML", data=data4)
modCprop1<-lme((Cyano.SumPhyto.ratio) ~ 1+I(invT-mean(invT)), random=~1|week, method="ML", data=data4)
modCprop2<-lme((Cyano.SumPhyto.ratio) ~ 1+I(invT-mean(invT)) + trophic.level,  random=~1|week, method="ML", data=data4)
modCprop4<-lme((Cyano.SumPhyto.ratio) ~ 1+I(invT-mean(invT)) * trophic.level, random=~1|week, method="ML", data=data4)

model.sel(modCprop0, modCprop1, modCprop2, modCprop4)

# basic plot
CyanoProp0 <- ggplot(data4, aes(x=invT, y=Cyano.SumPhyto.ratio, colour= factor(trophic.level))) + 
  geom_point(size = 3.5) +
  xlab('Temperature 1/k(T)') +
  ylab(expression('Cyanobacteria / Total Phytoplankton')) +
  scale_y_continuous(limits= c(-4.5, 3.5), breaks=c(0.0, 0.5, 1.0)) +
  scale_x_continuous(limits= c(38.0, 41.0), breaks=c(38.0, 39.0, 40.0, 41.0)) +
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(0,1)) +
  geom_abline(aes(intercept = (5.93299), slope = -0.142), color='#000000', lty=1, size= 1, data= data4) +
  geom_text(aes(label = 'D', x= 40.8, y=-Inf), vjust= -0.5, size=10, colour='black')+
  theme_bw()+
  scale_colour_manual(values=c('#000000', '#666666', '#CCCCCC'), breaks=c('P', 'PZ', 'PZN'), labels=c('1-TL', '2-TL', '3-TL'))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face='bold', colour='black', size=16, angle=0, vjust=0),
        axis.text.x = element_text(face='bold', colour='black', size=16, angle=0, vjust=-1),
        axis.title.x = element_text(colour='black', size=18),
        axis.title.y = element_text(colour='black', size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=16, face='bold'),
        legend.key = element_rect(colour = 'white', fill = 'white')
  )
CyanoProp0

############################################################
#### analysis for cyanos as proportion of total - weeks 7-9
data5 <- data4[(which(data4$week >=7)),]
modCprop0<-lme((Cyano.SumPhyto.ratio) ~ 1, random=~1|week, method="ML", data=data5)
modCprop1<-lme((Cyano.SumPhyto.ratio) ~ 1+I(invT-mean(invT)), random=~1|week, method="ML", data=data5)
modCprop2<-lme((Cyano.SumPhyto.ratio) ~ 1+I(invT-mean(invT)) + trophic.level,  random=~1|week, method="ML", data=data5)
modCprop4<-lme((Cyano.SumPhyto.ratio) ~ 1+I(invT-mean(invT)) * trophic.level, random=~1|week, method="ML", data=data5)

model.sel(modCprop0, modCprop1, modCprop2, modCprop4)

# basic plot
CyanoProp1 <- ggplot(data5, aes(x=invT, y=Cyano.SumPhyto.ratio, colour= factor(trophic.level))) + 
  geom_point(size = 3.5) +
  xlab('Temperature 1/k(T)') +
  ylab(expression('Cyanobacteria / Total Phytoplankton')) +
  scale_y_continuous(limits= c(-4.5, 3.5), breaks=c(0.0, 0.5, 1.0)) +
  scale_x_continuous(limits= c(38.0, 41.0), breaks=c(38.0, 39.0, 40.0, 41.0)) +
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(0,1)) +
  geom_abline(aes(intercept = (0.284-(-0.135*38.99289)), slope = -0.135), color='#000000', lty=1, size= 1, data= data4) +
  geom_text(aes(label = 'D', x= 40.8, y=-Inf), vjust= -0.5, size=10, colour='black')+
  theme_bw()+
  scale_colour_manual(values=c('#000000', '#666666', '#CCCCCC'), breaks=c('P', 'PZ', 'PZN'), labels=c('1-TL', '2-TL', '3-TL'))+
  theme(axis.line = element_line(colour = "black"),
        axis.text.y = element_text(face='bold', colour='black', size=16, angle=0, vjust=0),
        axis.text.x = element_text(face='bold', colour='black', size=16, angle=0, vjust=-1),
        axis.title.x = element_text(colour='black', size=18),
        axis.title.y = element_text(colour='black', size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black'),
        panel.background = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(color='black', size=16, face='bold'),
        legend.key = element_rect(colour = 'white', fill = 'white')
  )
CyanoProp1

>>>>>>> Stashed changes
