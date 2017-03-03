data1 <-data[(data$total.zoo.abundance.liter > 0),]

hist(log(data$total.zoo.abundance.liter))

## model selection
modZD0 <- lme(log(total.zoo.abundance.liter) ~ 1, random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")  
modZD1 <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, na.action=na.omit, method="ML")
modZD2 <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit)  
modZD4 <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modZD5 <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT)*trophic.level + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modZD6 <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 
modZD7 <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="ML", na.action=na.omit) 

model.sel(modZD0, modZD2, modZD4,modZD1, modZD5, modZD6, modZD7)

modZD2r <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level, random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit)
summary(modZD2r)
intervals(modZD2r)

#equation for line: 
intPZ <-  fixef(modZD2r)[1] - fixef(modZD2r)[3]*mean(data1$invTT)
intPZN <-  fixef(modZD2r)[1] + fixef(modZD2r)[4] - fixef(modZD2r)[3]*mean(data1$invTT)

# just checking if coefs are very different, since the model has such a similar aic
modZD1r <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data1, method="REML", na.action=na.omit)
summary(modZD1r)
intervals(modZD1r)

#Linear mixed-effects model fit by REML
#Data: data1 
#AIC      BIC    logLik
#219.1018 232.4189 -103.5509

#Random effects:
#  Formula: ~1 | Tank
#(Intercept) Residual
#StdDev:   0.2039838  1.01234

#Fixed effects: log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + I(invTT -      mean(invTT)) + trophic.level 
#                            Value Std.Error DF   t-value p-value
#(Intercept)            -0.0912853 0.1962011 51 -0.465264  0.6437
#I(invTi - invTT)        2.0201103 0.4653225 51  4.341312  0.0001
#I(invTT - mean(invTT))  0.7496764 0.5154921 17  1.454293  0.1641
#trophic.levelPZN       -0.4655405 0.2620151 17 -1.776770  0.0935
#
#Correlation: 
#  (Intr) I(T-iT I(TT-m
#                  I(invTi - invTT)       -0.456              
#                  I(invTT - mean(invTT))  0.076 -0.037       
#                  trophic.levelPZN       -0.541 -0.120 -0.093
#                  
#                  Standardized Within-Group Residuals:
#                    Min          Q1         Med          Q3         Max 
#                  -2.35505317 -0.74770896  0.04994241  1.00132387  1.42720663 
#                  
#                  Number of Observations: 72
#                  Number of Groups: 20 

mod.coefs <- augment(modZD2r, effect = "random")

### i changed these functions, they were not right. Because the data here are just for the two trophic levels (PZ and PZN), we can't just copy the expression from other code used on data with 3 trophic levels. the numbers for the coefficients, and even their formulas, didn't align because the 'base' trophic level here is the PZ level, not the P. 
ZD.funcPZ <- function(x) { (fixef(modZD2r)[1] - fixef(modZD2r)[3]*mean(mod.coefs$invTT)) + (fixef(modZD2r)[3])*x } # for trophic level 2
yvalsPZ <- ZD.funcPZ(mod.coefs$invTT)

ZD.funcPZN <- function(x) { (fixef(modZD2r)[1] + fixef(modZD2r)[4] - fixef(modZD2r)[3]*mean(mod.coefs$invTT)) + (fixef(modZD2r)[3])*x } # for trophic level 3
yvalsPZN <- ZD.funcPZN(mod.coefs$invTT)


#### Community Size
data1 <-data[(data$size > 0),]
hist(data$size)

modCS6r <- lme(log(size) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)) + I(invTT - mean(invTT))*trophic.level, random = ~ 1 | Tank, data=data, method="REML", na.action=na.omit)
summary(modCS6r)

#Linear mixed-effects model fit by REML
#Data: data 
#AIC    BIC    logLik
#419.291 454.92 -201.6455

#Random effects:
#  Formula: ~1 | Tank
#(Intercept)  Residual
#StdDev:   0.1223552 0.3213449

#Fixed effects: log(size) ~ 1 + I(invTi - invTT) * I(invTT - mean(invTT)) + I(invTT -      mean(invTT)) * trophic.level 
#Value Std.Error  DF   t-value p-value
#(Intercept)                             -0.3818456 0.0470874 621 -8.109297  0.0000
#I(invTi - invTT)                         0.2188096 0.0728203 621  3.004789  0.0028
#I(invTT - mean(invTT))                   0.3872953 0.2028975  14  1.908823  0.0770
#trophic.levelPZN                        -0.2936854 0.0682367  14 -4.303922  0.0007
#I(invTi - invTT):I(invTT - mean(invTT)) -0.5281246 0.3255477 621 -1.622265  0.1053
#I(invTT - mean(invTT)):trophic.levelPZN -0.2216658 0.2644378  14 -0.838253  0.4160
#Correlation: 
#  (Intr) I(T-iT I(nTT-m(TT)) tr.PZN I(-i-m
#                                      I(invTi - invTT)                        -0.332                                  
#                                      I(invTT - mean(invTT))                  -0.014  0.164                           
#                                      trophic.levelPZN                        -0.541 -0.224 -0.085                    
#                                      I(invTi - invTT):I(invTT - mean(invTT))  0.098 -0.377 -0.486        0.144       
#                                      I(invTT - mean(invTT)):trophic.levelPZN -0.052  0.104 -0.518        0.045 -0.142
#                                      
#                                      Standardized Within-Group Residuals:
#                                        Min         Q1        Med         Q3        Max 
#                                      -2.5933429 -0.6307137 -0.1008348  0.4538367  3.7263461 
#                                      
#                                      Number of Observations: 641
#                                      Number of Groups: 18 
#

# PP coefs
z <- 0.5 #invTi - invTT for each tank;

CS.PP.func <- function(x) { (fixef(modCS6r)[1] - fixef(modCS6r)[3]*mean(data$invTT) - fixef(modCS6r)[5]*z*mean(mod.coefs$invTT)) + (fixef(modCS6r)[3] + fixef(modCS6r)[5]*z)*x } #x = invTT # slope = -2.47
# use this function to compute yvals for plotting.
yvalsPZ <- CS.PP.func(mod.coefs$invTT)

# ZP coefs
CS.ZP.func <- function(x) { (fixef(modCS6r)[1] - fixef(modCS6r)[3]*mean(data$invTT) - fixef(modCS6r)[5]*z*mean(mod.coefs$invTT) + fixef(modCS6r)[4] - fixef(modCS6r)[6]*mean(mod.coefs$invTT)) + (fixef(modCS6r)[3] + fixef(modCS6r)[5]*z + fixef(modCS6r)[6])*x } #x = invTT #slope = -3.98
# use this function to compute yvals for plotting.
yvalsPZN <- CS.ZP.func(mod.coefs$invTT)


CS1 <-ggplot(data = data, aes(x = average.temp, y = log(size), ymin = 2)) + 
  theme_bw(22) +
  theme(legend.position = c(0.9, 0.9),
        legend.title=element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank()) +
  geom_point(aes(group=Tank, colour=trophic.level, shape=trophic.level), alpha=0.75, size=2) +
  xlab("Temperature [°C]") +
  ylab("Average community size log(mm)") +
  scale_shape_manual(values=c(15, 17)) +
  scale_colour_manual(limits=c("PZ", "PZN"), breaks=c("PZ", "PZN"), values=c("black", "grey60")) +
  scale_x_continuous(limits=c(15, 30),breaks=c(15, 20, 25, 30))
CS1

CS1 <- CS1 +
  geom_smooth(data=mod.coefs, aes(x = average.temp , y = yvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = average.temp, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey60', size = 1.5) +
  annotate("text", x=15, y=2, label="A", size=9)
CS1

CS2<-ggplot(data = data, aes(x = average.temp, y = log(size), ymin = 2)) + 
  theme_bw(22) +
  theme(legend.position = c(0.9, 0.9),
        legend.title=element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank()) +
  geom_point(aes(group=Tank, colour=trophic.level, shape=trophic.level), alpha=0.75, size=2) +
  xlab("Temperature [1/kT]") +
  ylab("Average community size log(mm)") +
  scale_shape_manual(values=c(15, 17)) +
  scale_colour_manual(limits=c("PZ", "PZN"), breaks=c("PZ", "PZN"), values=c("black", "grey60")) +
  scale_x_continuous(limits=c(15, 30),breaks=c(15, 20, 25, 30),labels=c('40.3', '39.6', '39.9', '38.3'))
CS2

CS2 <- CS2 +
  geom_smooth(data=mod.coefs, aes(x = average.temp , y = yvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = mod.coefs, aes(x = average.temp, y = yvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'grey60', size = 1.5) +
  annotate("text", x=15, y=2, label="A", size=9)
CS2

## extract gtable
A1 <- ggplot_gtable(ggplot_build(CS1)) #with invT
A2 <- ggplot_gtable(ggplot_build(CS2)) # with Celsius

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
Ax$grobs[[2]]$y <- Ax$grobs[[2]]$y - unit(1, "npc") + unit(0.13, "cm")

## modify existing row to be tall enough for axis
A$heights[[2]] <- A$heights[A2$layout[iA,]$t]

## add new axis
A <- gtable_add_grob(A, Ax, 2, 4, 2, 4)

## add new row for upper axis label
A <- gtable_add_rows(A, A2$heights[1], 1)
A <- gtable_add_grob(A, A2$grob[[6]], 1, 4, 2, 4)

# draw it
grid.draw(A)

### Average Community density
