## make data file

## objective here is to test whether zooplankton abundance, body size distribution or compositions shifted with temperature or predator treatment.


# data preparation --------------------------------------------------------

sizes <- read.csv("./data/CommunitySizes.csv")
sizes$Tank <- as.character(sizes$Tank)
sizes$week <- sizes$Week

## run data prep first from other file to get temps.wk and temps.mn
sizes.t <- left_join(sizes, temps.wk, by = c("week", "Tank")) # add weekly temps to sizes file
sizes.t <- sizes.t[,-c(1, 6:7)]

sizes.t <- left_join(sizes.t, temps.Tmn, by = c("Tank")) 

sizes.t2 <- as.tibble(sizes.t) %>%
  group_by(Tank, trophic.level) %>%
  summarise(., mean.temp = mean(temp.wk)) %>%
  arrange(trophic.level, mean.temp) %>%
  select(-(mean.temp)) 

sizes.t2$Tankn <- rep(c(1:9), 2)
sizes.t3 <- left_join(sizes.t, sizes.t2, by = c("Tank", "trophic.level")) 
sizes.t <- sizes.t3

sizes.t$invTi <-  1/((sizes.t$temp.wk + 273)*k) # average temp of the tank each week
sizes.t$invTT <-  1/((sizes.t$temp.Tmn + 273)*k) # average temp of the tank over all weeks
sizes.t$invTi <- as.numeric(as.character(sizes.t$invTi))
sizes.t$invTT <- as.numeric(as.character(sizes.t$invTT))

#sizes.t2 <- sizes.t[sizes.t$week >= '4',]

data.N <- as.tibble(data) %>%
  group_by(Tank, trophic.level, week) %>%
  summarise(., N = mean(total.zoo.abundance.liter)) %>%
  arrange(trophic.level, week) 

data$N1 <- 300*(data$total.zoo.abundance.liter) + 1
data$N <- 300*(data$total.zoo.abundance.liter)


# visualizing data --------------------------------------------------------

## abundance - waiting for raw abundance data to describe the distribution. it's count data, so makes me think we would be using a Poisson or negbinom model. But these have been transformed to density, so not whole numbers now. this is making it hard for me to test the distribution.

## from the data on abundance / L...
qqp(data[(data$trophic.level != "P"),]$N1, "lnorm")
qqp(data[(data$trophic.level != "P"),]$N1, "norm")
nbinom <- fitdistr(data[(data$trophic.level != "P"),]$N1, "Negative Binomial")
qqp(data[(data$trophic.level != "P"),]$N1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

gamma <- fitdistr(data[(data$trophic.level != "P"),]$N1, "gamma")


## how are the size data distributed?
sizes <- ggplot(data = sizes.t3, aes(x = log(size))) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  facet_grid(week~Tank)

sizes # pretty close to log normal
ggsave("sizes.png")
## seeing that we only have good coverage for week 8, so let's just use that.

data1 <- sizes.t[(sizes.t$week == "8"),] %>%
  mutate(taxon = species) %>%
  mutate(taxon = replace(taxon, species == "calanoid", "Copepod")) %>%
mutate(taxon = replace(taxon, species == "cyclopoid", "Copepod"))


sizes8 <- ggplot(data = data1, aes(x = log(size))) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  facet_grid(week~Tank)

sizes8

size.plot <- ggplot(data = data1, aes(x = -invTT, y = log(size))) + #, ymax = 1.2)
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(trophic.level~taxon) + ## this sets it up as facets
  geom_point(aes(group = taxon, shape = taxon), size = 2) + 
  scale_alpha("Tankn", guide = "none") +
  #theme(legend.position = c(0.88, 1), legend.text=element_text(size=6)) + 
  scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  scale_x_continuous("Celcius", sec.axis = sec_axis(~((1/(k*-.))-273))) +
  xlab("Temperature 1/kTi") +
  ylab("ln(size)") +
  geom_smooth(data = data1, method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = log(size)),  size = .8, color = alpha("steelblue", 0.5))
  
size.plot
ggsave("sizeswk8.png")
View(data1)

# Body Size Analysis ------------------------------------------------------

head(sizes.t3)
hist(sizes.t3$size)
hist(log(sizes.t3$size + 1))
hist(sizes.t3$week)
hist(as.numeric(sizes.t3$Tank))
hist(sizes.t3$temp.wk)
plot(sizes.t3$size ~ sizes.t3$week, col = sizes.t3$Tank)

hist((data$total.zoo.abundance.liter))

# ok, rather than plot every body size, I think we want to characterize the size distribution and plot that?

# our questions were: do temperature or predation change the zooplankton community. this can be assessed by looking at zooplankton density, species composition, and size distributions. 
# i need to think more about how to handle the temporal variation here.

# SIZE candidate model set -------------------------------------------------
## analyzing size for week 8
## suggests that predators reduce body sizes, across all measured individuals. temp too, but not significant. pooled across all species.
## with the interaction not significant, it's 'ignoring' the other trophic treatment.

modSF <- lm(log(size) ~ 1 + trophic.level + I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modS8 <- lm(log(size) ~ 1 + trophic.level + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modS7 <- lm(log(size) ~ 1 + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modS6 <- lm(log(size) ~ 1 + trophic.level, data=data1, na.action=na.omit)
modS5 <- lm(log(size) ~ 1, data=data1, na.action=na.omit)

model.sel(modS5, modS6, modS7, modS8, modSF)

## or, what if we use an averaged model: 
m.avgS <- model.avg(modS6, modS8)
coefficients(m.avgS)
confint(m.avgS)

size.funcZ <- function(x) { coefficients(m.avgS)[1] + coefficients(m.avgS)[3]*x }
size.funcZN <- function(x) { coefficients(m.avgS)[1] + coefficients(m.avgS)[2] + coefficients(m.avgS)[3]*x }
z.vals <- size.funcZ(data1[(data1$trophic.level =='PZ'),]$invTT - mean(data1$invTT))
zn.vals <- size.funcZN(data1[(data1$trophic.level =='PZN'),]$invTT - mean(data1$invTT))



Fig3B <- ggplot(data = data1, aes(x = -invTT, y = log(size))) + #, ymax = 1.2)
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(trophic.level~.) + ## this sets it up as facets
  geom_point(aes(group = taxon, shape = taxon), size = 2) + 
  scale_alpha("Tankn", guide = "none") +
  #theme(legend.position = c(0.88, 1), legend.text=element_text(size=6)) + 
  scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  scale_x_continuous("Temperature", sec.axis = sec_axis(~((1/(k*-.))-273))) +
  #xlab("Temperature 1/kTi") +
  ylab("ln(size)") +
  geom_smooth(data = data1[(data1$trophic.level =='PZ'),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = z.vals),  size = .8, color = alpha("steelblue", 0.5)) +
  geom_smooth(data = data1[(data1$trophic.level =='PZN'),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = zn.vals),  size = .8, color = alpha("steelblue", 0.5))

Fig3B

ggsave("Fig 3B.png")


## daphnia sizes
## i think we want a figure like this, but with ses on the datapoints. 
## start with the zooplankton size data i just added to the data folder. 

data1 <- data[(data$week == "8"),] #from main datafile
Fig3D <- ggplot(data = data1[(data1$trophic.level != "P"),], aes(x = -invTT, y = (Daphnia.mature.size), ymax = .8, ymin = -0.2)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(trophic.level~.) + ## this sets it up as facets
  geom_point(size = 2) + 
  geom_point(aes(x = -invTT, y = (community.size)), shape = 2) +
  scale_alpha("Tankn", guide = "none") +
  #theme(legend.position = c(0.88, 1), legend.text=element_text(size=6)) + 
  #scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  scale_x_continuous("Temperature", sec.axis = sec_axis(~((1/(k*-.))-273))) +
  #xlab("Temperature 1/kTi") +
  ylab("Mature Daphnia ln(length)") +
  geom_smooth(data = data1[(data1$trophic.level =='PZ'),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = z.vals),  size = .8, color = alpha("steelblue", 0.5)) +
  geom_smooth(data = data1[(data1$trophic.level =='PZN'),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = zn.vals),  size = .8, color = alpha("steelblue", 0.5))

Fig3D


## abundance analysis
require(car)
require(MASS)

### used a fixed effects model only because i'm not aware of methods for the mixed effects models on an poisson distribution
### following http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
### mean < 5, so...
library("mlmRev")
library(lme4)
GHQ <- glmer(N1 ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)) + (1 | Tank), data=data[(data$trophic.level != "P"),], family = binomial(link = "logit"), nAGQ = 25)
summary(PQL)
summary(GHQ)

## check y vals... refresh on what this method is even doing...



### functionf or checking for overdispersion
model<-GHQ
overdisp_fun <- function(model) {
  ## number of variance parameters in an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m) * (nrow(m) + 1)/2
  }
  # The next two lines calculate the residual degrees of freedom
  model.df <- sum(sapply(VarCorr(model), vpars)) + length(fixef(model))
  rdf <- nrow(model.frame(model)) - model.df
  # extracts the Pearson residuals
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  # Generates a p-value. If less than 0.05, the data are overdispersed.
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = prat, rdf = rdf, p = pval)
}
# the data are overdisperssed. crap.


Fig3A <- 
  size.plot + 
  geom_smooth(data = subset(sizes.t3, week == "9"), method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTi, y = log(size + 1), group = Tank),  size = .8, color = alpha("steelblue", 0.5)) + 
  geom_smooth(data = subset(data1, trophic.level == "P"), aes(x = -invTT, y = NvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data1, trophic.level == "PZ"), aes(x = -invTT, y = NvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data1, trophic.level == "PZN"), aes(x = -invTT, y = NvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) 


# Possibly recovered code -------------------------------------------------
#####
#### possibly recovered code?

# FIGURE 3 CD -------------------------------------------------------------
### THIS CODE MAKES FIGURE 3 IN MANUSCRIPT
## use 'final size' spreadsheet
sz <- read.csv("./data/Finalzpsize.csv")
sz$Tank <- as.character(sz$Tank)
#sz$week <- sz$Week
#View(sz)


sz2 <- sz %>%
  mutate(taxon = species) %>%
  mutate(taxon = replace(taxon, species == "calanoid", "Copepod")) %>%
  mutate(taxon = replace(taxon, species == "cyclopoid", "Copepod")) %>%
  mutate(trophic.level = as.character(treatment)) %>%
  mutate(trophic.level = replace(trophic.level, treatment == "PP+Z", "PZ")) %>%
  mutate(trophic.level = replace(trophic.level, treatment == "PP+Z+N", "PZN")) %>%
  filter(stage != "larvae") %>%
  filter(stage != "egg") %>%
  filter(size.in.cm > 0)

View(sz2)

testF <- sz2 %>%
  group_by(., Tank, week, taxon) %>%
  summarise(., length(size.in.cm)) %>%
  arrange(as.numeric(week), as.numeric(Tank)) %>%
  filter(taxon == "Daphnia")
dim(testF)

sz2.t <- left_join(sz2, temps.wk, by = c("week", "Tank")) # add weekly temps to sizes file
sz2.t <- sz2.t[,-c(1:2, 5:6)]

sz2.t1 <- left_join(sz2.t, temps.Tmn, by = c("Tank")) 

sz2.t2 <- as.tibble(sz2.t1) %>%
  group_by(Tank, trophic.level) %>%
  summarise(., mean.temp = mean(temp.wk)) %>%
  arrange(trophic.level, mean.temp)
  #select(-mean.temp) 

sz2.t2$Tankn <- rep(c(1:9), 2)
sz2.t3 <- left_join(sz2.t1, sz2.t2, by = c("Tank", "trophic.level")) 

sz2.t3$invTi <-  1/((sz2.t3$temp.wk + 273)*k) # average temp of the tank each week
sz2.t3$invTT <-  1/((sz2.t3$temp.Tmn + 273)*k) # average temp of the tank over all weeks
sz2.t3$invTi <- as.numeric(as.character(sz2.t3$invTi))
sz2.t3$invTT <- as.numeric(as.character(sz2.t3$invTT))

sizes <- ggplot(data = sz2.t3, aes(x = log(size.in.cm))) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  facet_grid(week~Tank)

data1 <- sz2.t3[(sz2.t3$week == "8"),]
sizes8 <- ggplot(data = data1, aes(x = log(size.in.cm))) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  facet_grid(week~Tank)
sizes8

size.plot <- ggplot(data = data1, aes(x = -invTT, y = log(size.in.cm))) + #, ymax = 1.2)
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(trophic.level~taxon) + ## this sets it up as facets
  facet_grid(trophic.level~.) + ## this sets it up as facets
  geom_point(aes(group = taxon, shape = taxon), size = 2) + 
  scale_alpha("Tankn", guide = "none") +
  #theme(legend.position = c(0.88, 1), legend.text=element_text(size=6)) +
  scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  scale_x_continuous("Celcius", sec.axis = sec_axis(~((1/(k*-.))-273))) +
  xlab("Temperature 1/kTi") +
  ylab("Length ln(cm)")

#size.plot

data1$size <- data1$size.in.cm
## this is what I want, I think [sept 27 2017]
Fig3B <- ggplot(data = data1, aes(x = -invTT, y = log(size))) + #, ymax = 1.2)
    theme_bw() +
    theme(legend.position = "none") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(strip.background = element_rect(colour="white", fill="white")) +
    facet_grid(.~trophic.level, labeller=labeller(trophic.level = labels)) + ## this sets it up as facets
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    geom_point(aes(group = taxon, shape = taxon, color = taxon), size = 2, alpha = 0.5) + 
    scale_alpha("Tankn", guide = "none") +
    scale_colour_grey(start = 0, end = 0.6, guide = "none") +
    scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~((1/(k*-.))-273))) +
    ylab("ln(size)") + #, name = "deg Celcius"
    ylab("Length ln(cm)") 
    #geom_smooth(data = data1[(data1$trophic.level=="PZ"),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = log(size)),  size = .8, color = alpha("steelblue", 0.5)) +
    #geom_smooth(data = data1[(data1$trophic.level=="PZN"),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = SvalsPZN),  size = .8, color = alpha("steelblue", 0.5))
Fig3B
  ggsave("Fig3B.png", device = "png", width =4, height = 3) 
  
  ## stats for size data in figure 3:
  modSF <- lm(log(size) ~ 1 + trophic.level*taxon + I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
  modS8 <- lm(log(size) ~ 1 + trophic.level*taxon + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
  modS7 <- lm(log(size) ~ 1 + I(invTT - mean(invTT))*taxon, data=data1, na.action=na.omit)
  modS6 <- lm(log(size) ~ 1 + trophic.level*taxon, data=data1, na.action=na.omit)
  modS5 <- lm(log(size) ~ 1 + taxon, data=data1, na.action=na.omit)
  modSF2 <- lm(log(size) ~ 1 + trophic.level + I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
  modS4 <- lm(log(size) ~ 1 + trophic.level + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
  modS3 <- lm(log(size) ~ 1 + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
  modS2 <- lm(log(size) ~ 1 + trophic.level , data=data1, na.action=na.omit)
  modS1 <- lm(log(size) ~ 1, data=data1, na.action=na.omit)
  
  
  model.sel(modS5, modS6, modS7, modS8, modSF, modS1, modS2, modS3, modS4, modSF2)
  
  ## or, what if we use an averaged model: 
  m.avgS <- model.avg(modS6, modS8)
  coefficients(m.avgS)
  confint(m.avgS)
  
 ## size figure without temperature: 
  
  labels2 <- c("copepod", "daphnia")
  xlabels <- c("AG", "AGP")
  
Fig3B <- ggplot(data = data1, aes(x = trophic.level, y = log(size))) + #, ymax = 1.2)
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          strip.background = element_rect(colour="white", fill="white"),
          #strip.background = element_blank(), 
          strip.text = element_blank()) +
    geom_boxplot() +
    facet_grid(.~taxon, labeller=labeller(taxon = labels2)) + ## this sets it up as facets
    scale_alpha("Tankn", guide = "none") +
    scale_colour_grey(start = 0, end = 0.6, guide = "none") +
    scale_x_discrete(labels= xlabels) +
    ylab("Length ln(cm)") +
    xlab("Food Chain Length")

# Size figure box plot ----------------------------------------------------


  Fig3B
  ggsave("Fig3B.png", device = "png", width =4, height = 3) 

### recovered from github: 

  # ABUNDANCE ---------------------------------------------------------------
    ## abundance - waiting for raw abundance data to describe the distribution. it's count data, so makes me think we would be using a Poisson or negbinom model.
    
qqp(data[(data$trophic.level != "P"),]$N1, "lnorm")
hist(log((data[(data$trophic.level != "P"),]$N1)))

# Dec 2017: Neg binomial
mod1 <- glm.nb(log(N1) ~ 1 + (trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT))), data=data[(data$trophic.level != "P"),]) 
    
## i think we're going to go with the lognormal distribution
### might just ax the random int...or test for it: 
modNF <- lme(log(N1) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], method="ML", na.action=na.omit) 

modNa <- lm(log(N1) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data[(data$trophic.level != "P"),], na.action=na.omit)
                             
anova(modNF, modNa) 
                             
modNF <- lme(log(N1) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN8 <- lme(log(N1) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN7 <- lme(log(N1) ~ 1 + I(invTi - invTT) + trophic.level * I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN6 <- lme(log(N1) ~ 1 + trophic.level*I(invTi - invTT), data=data[(data$trophic.level != "P"),], random = ~ 1 | Tank, na.action=na.omit)
modN5 <- lme(log(N1) ~ 1 + I(invTi - invTT) + trophic.level, random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
#modN4 <- lm(log(N1) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN3 <- lme(log(N1) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN2 <- lme(log(N1) ~ 1 + I(invTi - invTT), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN1 <- lme(log(N1) ~ 1 + trophic.level, random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN0 <- lme(log(N1) ~ 1, random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)

model.sel(modN0, modN1, modN2, modN3, modN5, modN6, modN7, modN8, modNF)

# with random effect included, modNF is best
coefficients(modNF)
intervals(modNF)
                             
## or, what if we use an averaged model: 
m.avgN <- model.avg(modNF, modN8)
coefficients(m.avgN)
confint(m.avgN)


# figure 3A Density week 8 ------------------------------------------------


## analysis for just week 8:
data1 <- data[(data$trophic.level != "P" & data$week == "8"),]
modN3 <- lm(log(N1) ~ 1 + invTi*trophic.level, data=data1, na.action=na.omit)
modN2 <- lm(log(N1) ~ 1 + invTi, data=data1, na.action=na.omit)
modN1 <- lm(log(N1) ~ 1 + trophic.level, data=data1, na.action=na.omit)
modN0 <- lm(log(N1) ~ 1, data=data1, na.action=na.omit)

model.sel(modN0, modN1, modN2, modN3)
coef(modN2)
confint(modN2)

mod.coefs <- augment(modN2)  
mod.coefs1 <- left_join(data1, mod.coefs)
          
labels <- c(P = "A", PZ = "AG", PZN = "AGP")
xlab <- expression(paste('Temperature (',~degree,'C)',sep=''))

N.plot <- ggplot(data=mod.coefs1, aes(x = -invTi, y = log.N1., ymin = 0, ymax = 9)) + #
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  #facet_grid(.~trophic.level, labeller=labeller(trophic.level = labels)) + ## this sets it up as facets
  #geom_ribbon(aes(ymin = .fitted - .se.fit, ymax = .fitted + .se.fit), fill = "grey70") +
  geom_point(aes(group = as.character(Tankn), color = as.character(Tankn), shape = as.factor(trophic.level), alpha = Tankn), size = 4) + 
  geom_smooth(aes(x=-invTi, y = log.N1.), method = "lm") +
  geom_line(aes(x=-invTi, y = .fitted)) +
  scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none") +
  scale_alpha("Tankn", guide = "none") +
  scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~((1/(k*-.))-273), name = xlab)) +
  #theme(legend.position = c(0.88, 0.15), legend.text=element_text(size=6)) + 
  #scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  xlab("Temperature 1/kTi") +
  ylab("Density ln((N+1) / 300 L)")


N.plot

ggsave("Nplot.png", device = "png", width = 4, height = 3) # save for appendix
                             
N.PZ <- function(x) { coefficients(m.avgN)[1] + coefficients(m.avgN)[3]*x}
NvalsPZ <- N.PZ(data[(data$trophic.level=="PZ"),]$invTT - mean(data$invTT))
N.PZN <- function(x) { (coefficients(m.avgN)[1] + coefficients(m.avgN)[3]) + (coefficients(m.avgN)[3])*x}
NvalsPZN <- N.PZN(data[(data$trophic.level=="PZN"),]$invTT- mean(data$invTT))

Fig3A <- 
N.plot +
  geom_line(aes(x = fitted(modN2), y = fitted(modN2)))
  
## just fit a line to the fitted values, and then use geom_ribbon for CIs.


  geom_smooth(method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTi, y = log(N1), group = week, colour = "week")) #color = alpha("steelblue", 0.5)
# scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none")
# stop here - no 'among tank' temp term in best models
#geom_smooth(method = "lm", se = FALSE, aes(group = Tank), color = "gray40", alpha = 0.23, size = .8) +
#geom_smooth(data = subset(data, trophic.level == "PZ"), aes(x = -invTT, y = NvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 1, size = 1.5) +
#geom_smooth(data = subset(data, trophic.level == "PZN"), aes(x = -invTT, y = NvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', linetype = 1, size = 1.5) 

ggsave("Fig3A.png", device = "png", width = 4, height = 3)     

### abundance by species

data1 <- data[(data$trophic.level != "P" & data$week == "8"),]
hist(data1$abundance.Daphnia)
modD3 <- lm(log(abundance.Daphnia) ~ 1 + invTi*trophic.level, data=data1, na.action=na.omit)
modD2 <- lm(log(abundance.Daphnia) ~ 1 + invTi, data=data1, na.action=na.omit)
modD1 <- lm(log(abundance.Daphnia) ~ 1 + trophic.level, data=data1, na.action=na.omit)
modD0 <- lm(log(abundance.Daphnia) ~ 1, data=data1, na.action=na.omit)

results <- model.sel(modD3, modD2, modD1, modD0)

D.plot <- ggplot(data=data[(data$trophic.level != "P"),], aes(x = -invTi, y = log.N1., ymin = 0, ymax = 9)) + #
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +
  #facet_grid(.~trophic.level, labeller=labeller(trophic.level = labels)) + ## this sets it up as facets
  #geom_ribbon(aes(ymin = .fitted - .se.fit, ymax = .fitted + .se.fit), fill = "grey70") +
  geom_point(aes(group = as.character(Tankn), color = as.character(Tankn), shape = as.factor(trophic.level), alpha = Tankn), size = 4) + 
  geom_smooth(aes(x=-invTi, y = log.N1.), method = "lm") +
  geom_line(aes(x=-invTi, y = .fitted)) +
  scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none") +
  scale_alpha("Tankn", guide = "none") +
  scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~((1/(k*-.))-273), name = xlab)) +
  #theme(legend.position = c(0.88, 0.15), legend.text=element_text(size=6)) + 
  #scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  xlab("Temperature 1/kTi") +
  ylab("Density ln((N+1) / 300 L)")


N.plot