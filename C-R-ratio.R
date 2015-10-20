data1 <- read.csv("C:/Dropbox/Dropbox/OConnor Lab/manuscripts/Jessica Tank experiment 2012/DATA MARY SHOULD USE TODAY/Consumer.Resource.csv")
View(data1)
attach(data1)

k <- 8.617342*10^-5  # eV/K
data1$invT <-  1/((data1$average.temp + 273)*k)
data1$invT <- as.numeric(as.character(data1$invT))


hist(data1$Pb.Zb)
hist(log(data1$Pb.Zb))
hist(log(data1$Pb.Zb+1))

##Analysis Consumer.Resource Ratio = Phyto.biomass / Zoo.biomass
modCR0<-lme(log(Pb.Zb+1) ~ 1, random=~1|week, method="ML", data=data1, na.action=na.omit)
modCR1<-lme(log(Pb.Zb+1) ~ 1+I(invT-mean(invT)), random=~1|week, method="ML", data=data1, na.action=na.omit)
modCR2<-lme(log(Pb.Zb+1) ~ 1+I(invT-mean(invT)) + trophic.level,  random=~1|week, method="ML", data=data1, na.action=na.omit)
modCR4<-lme(log(Pb.Zb+1) ~ 1+I(invT-mean(invT)) * trophic.level, random=~1|week, method="ML", data=data1, na.action=na.omit)

model.sel(modCR0, modCR1, modCR2, modCR4)
#Model selection table 
#(Int) inT-men(inT) trp.lvl I(inT-men(inT)):trp.lvl df logLik   AICc  delta weight
#4 2.743 1.4000       +       +                       6   -97.460 208.4  0.00 0.995 
#3 2.638 0.4092       +                               5  -104.008 219.0 10.66 0.005 
#1 2.939                                              3  -109.196 224.8 16.42 0.000 
#2 2.978 0.5368                                       4  -108.104 224.9 16.51 0.000 
#Random terms (all models): 
#  ‘1 | week’

anova(modCR4, modCR2)
#Model df      AIC      BIC     logLik   Test  L.Ratio p-value
#modCR4     1  6 206.9201 219.9664  -97.46004                        
#modCR2     2  5 218.0157 228.8876 -104.00784 1 vs 2 13.09559   3e-04

anova(modCR2, modCR1)
#Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#modCR2     1  5 218.0157 228.8876 -104.0078                        
#modCR1     2  4 224.2074 232.9050 -108.1037 1 vs 2 8.191741  0.0042

# Model refitting
modCR4<-lme(log(Pb.Zb+1) ~ 1+invT * trophic.level, random=~1|week, method="REML", data=data1, na.action=na.omit)
anova(modCR1, modCR0)
#Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#modCR1     1  4 224.2074 232.9050 -108.1037                        
#modCR0     2  3 224.3929 230.9161 -109.1965 1 vs 2 2.185519  0.1393

summary(modCR4)
#Linear mixed-effects model fit by maximum likelihood
#Data: data1 
#AIC      BIC    logLik
#206.9201 219.9664 -97.46004

#Random effects:
#  Formula: ~1 | week
#(Intercept)  Residual
#StdDev:   0.7478776 0.9914631

#Fixed effects: log(Pb.Zb + 1) ~ 1 + I(invT - mean(invT)) * trophic.level 
#Value Std.Error DF   t-value p-value
#(Intercept)                            2.7434007 0.3594796 56  7.631590   0e+00
#I(invT - mean(invT))                   1.4002288 0.3802251 56  3.682631   5e-04
#trophic.levelPZN                       0.9490274 0.2678273 56  3.543430   8e-04
#I(invT - mean(invT)):trophic.levelPZN -1.9163147 0.5105603 56 -3.753356   4e-04
#Correlation: 
#  (Intr) I(nT-m(T)) tr.PZN
#I(invT - mean(invT))                   0.117                  
#trophic.levelPZN                      -0.292 -0.096           
#I(invT - mean(invT)):trophic.levelPZN -0.072 -0.655     -0.068

#Standardized Within-Group Residuals:
#  Min         Q1        Med         Q3        Max 
#-2.5972347 -0.6744163 -0.1839581  0.8093832  1.8700109 

#Number of Observations: 65
#Number of Groups: 6 

#Model refitting
modCR4<-lme(log(Pb.Zb+1) ~ 1+invT * trophic.level, random=~1|week, method="REML", data=data1, na.action=na.omit)
summary(modCR4)
#Linear mixed-effects model fit by REML
#Data: data1 
#AIC      BIC    logLik
#208.3081 220.9734 -98.15407

#Random effects:
#  Formula: ~1 | week
#(Intercept) Residual
#StdDev:   0.8341903 1.017091

#Fixed effects: log(Pb.Zb + 1) ~ 1 + invT * trophic.level 
#Value Std.Error DF   t-value p-value
#(Intercept)           -52.98823 14.822879 56 -3.574760   7e-04
#invT                    1.42038  0.378701 56  3.750656   4e-04
#trophic.levelPZN       76.48350 19.936398 56  3.836375   3e-04
#invT:trophic.levelPZN  -1.92489  0.507561 56 -3.792427   4e-04
#Correlation: 
#  (Intr) invT   tr.PZN
#invT                  -1.000              
#trophic.levelPZN      -0.652  0.652       
#invT:trophic.levelPZN  0.654 -0.654 -1.000

#Standardized Within-Group Residuals:
#  Min         Q1        Med         Q3        Max 
#-2.5040691 -0.6452062 -0.1711511  0.7826078  1.8308283 

#Number of Observations: 65
#Number of Groups: 6 

#### Figures###
# basic plot
Pb.Zb0 <- ggplot(data1, aes(x=invT, y=log(Pb.Zb + 1), colour= factor(trophic.level))) + 
  #geom_line()+
  geom_point(size = 3.5) +
  xlab('Temperature 1/k(T)') +
  ylab(expression('Pb / Zb ratio')) +
  scale_y_continuous(limits= c(0.0, 6.5), breaks=c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0))+
  scale_x_continuous(limits= c (38.0, 41.0), breaks=c(38.0, 39.0, 40.0, 41.0))+
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(0.0, 6.5))+
  geom_abline(aes(intercept = -52.98823, slope = 1.42038), color='#666666', lty=1, size= 1, data=data1)+
  geom_abline(aes(intercept = (-52.98823+76.48350), slope = (1.42038-1.92489)), color='#CCCCCC', lty=1, size= 1, data=data1)+
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
Pb.Zb0

# plot with trasformed axis
Pb.Zb1 <- ggplot(data1, aes(x=invT, y=log(Pb.Zb + 1), colour= factor(trophic.level))) + 
  #geom_line()+
  geom_point(size = 3.5) +
  xlab('Temperature C') +
  ylab(expression('Pb / Zb ratio')) +
  scale_y_continuous(limits= c(0.0, 6.5), breaks=c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0))+
  scale_x_continuous(limits=c(38, 41),breaks=c(38, 39, 40, 41),labels=c('30.0', '24.5', '16.2', '12.8')) +
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(0.0, 6.5))+
  geom_abline(aes(intercept = -52.98823, slope = 1.42038), color='#666666', lty=1, size= 1, data=data1)+
  geom_abline(aes(intercept = (-52.98823+76.48350), slope = (1.42038-1.92489)), color='#CCCCCC', lty=1, size= 1, data=data1)+
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
Pb.Zb1

## extract gtable
A1 <- ggplot_gtable(ggplot_build(Pb.Zb0)) #with invT
A2 <- ggplot_gtable(ggplot_build(Pb.Zb1)) # with Celsius

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


########################## Zb + Noto carbon / Pb
hist(data1$Zb.Nb.Pb)
hist(log(data1$Zb.Nb.Pb))
hist(log(data1$Zb.Nb.Pb+1))

##Analysis Consumer.Resource Ratio = Phyto.biomass / Zoo.biomass
modZbNb.Pb0<-lme(log(Zb.Nb.Pb+1) ~ 1, random=~1|week, method="ML", data=data1, na.action=na.omit)
modZbNb.Pb1<-lme(log(Zb.Nb.Pb+1) ~ 1+I(invT-mean(invT)), random=~1|week, method="ML", data=data1, na.action=na.omit)
modZbNb.Pb2<-lme(log(Zb.Nb.Pb+1) ~ 1+I(invT-mean(invT)) + trophic.level,  random=~1|week, method="ML", data=data1, na.action=na.omit)
modZbNb.Pb4<-lme(log(Zb.Nb.Pb+1) ~ 1+I(invT-mean(invT)) * trophic.level, random=~1|week, method="ML", data=data1, na.action=na.omit)

model.sel(modZbNb.Pb0, modZbNb.Pb1, modZbNb.Pb2)
#Model selection table 
#(Int) inT-men(inT) trp.lvl I(inT-men(inT)):trp.lvl df logLik   AICc  delta weight
#4 2.743 1.4000       +       +                       6   -97.460 208.4  0.00 0.995 
#3 2.638 0.4092       +                               5  -104.008 219.0 10.66 0.005 
#1 2.939                                              3  -109.196 224.8 16.42 0.000 
#2 2.978 0.5368                                       4  -108.104 224.9 16.51 0.000 
#Random terms (all models): 
#  ‘1 | week’

anova(modCR4, modCR2)
#Model df      AIC      BIC     logLik   Test  L.Ratio p-value
#modCR4     1  6 206.9201 219.9664  -97.46004                        
#modCR2     2  5 218.0157 228.8876 -104.00784 1 vs 2 13.09559   3e-04

anova(modCR2, modCR1)
#Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#modCR2     1  5 218.0157 228.8876 -104.0078                        
#modCR1     2  4 224.2074 232.9050 -108.1037 1 vs 2 8.191741  0.0042

# Model refitting
modCR4<-lme(log(Pb.Zb+1) ~ 1+invT * trophic.level, random=~1|week, method="REML", data=data1, na.action=na.omit)
anova(modCR1, modCR0)
#Model df      AIC      BIC    logLik   Test  L.Ratio p-value
#modCR1     1  4 224.2074 232.9050 -108.1037                        
#modCR0     2  3 224.3929 230.9161 -109.1965 1 vs 2 2.185519  0.1393

summary(modCR4)
#Linear mixed-effects model fit by maximum likelihood
#Data: data1 
#AIC      BIC    logLik
#206.9201 219.9664 -97.46004

#Random effects:
#  Formula: ~1 | week
#(Intercept)  Residual
#StdDev:   0.7478776 0.9914631

#Fixed effects: log(Pb.Zb + 1) ~ 1 + I(invT - mean(invT)) * trophic.level 
#Value Std.Error DF   t-value p-value
#(Intercept)                            2.7434007 0.3594796 56  7.631590   0e+00
#I(invT - mean(invT))                   1.4002288 0.3802251 56  3.682631   5e-04
#trophic.levelPZN                       0.9490274 0.2678273 56  3.543430   8e-04
#I(invT - mean(invT)):trophic.levelPZN -1.9163147 0.5105603 56 -3.753356   4e-04
#Correlation: 
#  (Intr) I(nT-m(T)) tr.PZN
#I(invT - mean(invT))                   0.117                  
#trophic.levelPZN                      -0.292 -0.096           
#I(invT - mean(invT)):trophic.levelPZN -0.072 -0.655     -0.068

#Standardized Within-Group Residuals:
#  Min         Q1        Med         Q3        Max 
#-2.5972347 -0.6744163 -0.1839581  0.8093832  1.8700109 

#Number of Observations: 65
#Number of Groups: 6 

#Model refitting
modCR4<-lme(log(Pb.Zb+1) ~ 1+invT * trophic.level, random=~1|week, method="REML", data=data1, na.action=na.omit)
summary(modCR4)
#Linear mixed-effects model fit by REML
#Data: data1 
#AIC      BIC    logLik
#208.3081 220.9734 -98.15407

#Random effects:
#  Formula: ~1 | week
#(Intercept) Residual
#StdDev:   0.8341903 1.017091

#Fixed effects: log(Pb.Zb + 1) ~ 1 + invT * trophic.level 
#Value Std.Error DF   t-value p-value
#(Intercept)           -52.98823 14.822879 56 -3.574760   7e-04
#invT                    1.42038  0.378701 56  3.750656   4e-04
#trophic.levelPZN       76.48350 19.936398 56  3.836375   3e-04
#invT:trophic.levelPZN  -1.92489  0.507561 56 -3.792427   4e-04
#Correlation: 
#  (Intr) invT   tr.PZN
#invT                  -1.000              
#trophic.levelPZN      -0.652  0.652       
#invT:trophic.levelPZN  0.654 -0.654 -1.000

#Standardized Within-Group Residuals:
#  Min         Q1        Med         Q3        Max 
#-2.5040691 -0.6452062 -0.1711511  0.7826078  1.8308283 

#Number of Observations: 65
#Number of Groups: 6 

#### Figures###
# basic plot
Pb.Zb0 <- ggplot(data1, aes(x=invT, y=log(Pb.Zb + 1), colour= factor(trophic.level))) + 
  #geom_line()+
  geom_point(size = 3.5) +
  xlab('Temperature 1/k(T)') +
  ylab(expression('Pb / Zb ratio')) +
  scale_y_continuous(limits= c(0.0, 6.5), breaks=c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0))+
  scale_x_continuous(limits= c (38.0, 41.0), breaks=c(38.0, 39.0, 40.0, 41.0))+
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(0.0, 6.5))+
  geom_abline(aes(intercept = -52.98823, slope = 1.42038), color='#666666', lty=1, size= 1, data=data1)+
  geom_abline(aes(intercept = (-52.98823+76.48350), slope = (1.42038-1.92489)), color='#CCCCCC', lty=1, size= 1, data=data1)+
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
Pb.Zb0

# plot with trasformed axis
Pb.Zb1 <- ggplot(data1, aes(x=invT, y=log(Pb.Zb + 1), colour= factor(trophic.level))) + 
  #geom_line()+
  geom_point(size = 3.5) +
  xlab('Temperature C') +
  ylab(expression('Pb / Zb ratio')) +
  scale_y_continuous(limits= c(0.0, 6.5), breaks=c(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0))+
  scale_x_continuous(limits=c(38, 41),breaks=c(38, 39, 40, 41),labels=c('30.0', '24.5', '16.2', '12.8')) +
  coord_cartesian(xlim=c(38.0, 41.0), ylim=c(0.0, 6.5))+
  geom_abline(aes(intercept = -52.98823, slope = 1.42038), color='#666666', lty=1, size= 1, data=data1)+
  geom_abline(aes(intercept = (-52.98823+76.48350), slope = (1.42038-1.92489)), color='#CCCCCC', lty=1, size= 1, data=data1)+

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
Pb.Zb1

## extract gtable
A1 <- ggplot_gtable(ggplot_build(Pb.Zb0)) #with invT
A2 <- ggplot_gtable(ggplot_build(Pb.Zb1)) # with Celsius

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