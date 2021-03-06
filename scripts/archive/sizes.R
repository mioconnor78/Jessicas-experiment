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

sizes.t2 <- sizes.t[sizes.t$week >= '4',]

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



## abundance analysis
require(car)
require(MASS)

### used a fixed effects model only because i'm not aware of methods for the mixed effects models on an poisson distribution
### following http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
### mean < 5, so...
library("mlmRev")
library(lmer)
GHQ <- glmer(N1 ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)) + (1 | Tank), data=data[(data$trophic.level != "P"),], family = binomial(link = "logit"), nAGQ = 25)
summary(PQL)



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


modNF <- lme(log(total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, quasipoisson(link = "log"), data=data[(data$trophic.level != "P"),], method="ML", na.action=na.omit) 
modNa <- glm((total.zoo.abundance.liter) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), quasipoisson(link = "log"), data=data[(data$trophic.level != "P"),], na.action=na.omit)

modNPPF <- lm(log(NPP2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modNPP8 <- lm(log(NPP2) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modNPP7 <- lm(log(NPP2) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modNPP6 <- lm(log(NPP2) ~ 1 + trophic.level*I(invTi - invTT), data=data1, na.action=na.omit)
modNPP5 <- lm(log(NPP2) ~ 1 + I(invTi - invTT) + trophic.level, data=data1, na.action=na.omit)
modNPP4 <- lm(log(NPP2) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modNPP3 <- lm(log(NPP2) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modNPP2 <- lm(log(NPP2) ~ 1 + I(invTi - invTT), data=data1, na.action=na.omit)
modNPP1 <- lm(log(NPP2) ~ 1 + trophic.level, data=data1, na.action=na.omit)
modNPP0 <- lm(log(NPP2) ~ 1, data=data1, na.action=na.omit)

model.sel(modNPP0, modNPP1, modNPP2, modNPP3, modNPP4, modNPP5, modNPP6, modNPP7, modNPP8, modNPPF)

## or, what if we use an averaged model: 
m.avgN <- model.avg(modNPP8, modNPP3)
confint(m.avgN)


Fig3A <- 
  size.plot + 
  geom_smooth(data = subset(sizes.t3, week == "9"), method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTi, y = log(size + 1), group = Tank),  size = .8, color = alpha("steelblue", 0.5)) + 
  geom_smooth(data = subset(data1, trophic.level == "P"), aes(x = -invTT, y = NvalsP), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data1, trophic.level == "PZ"), aes(x = -invTT, y = NvalsPZ), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) +
  geom_smooth(data = subset(data1, trophic.level == "PZN"), aes(x = -invTT, y = NvalsPZN), method = "lm", se = FALSE, inherit.aes = FALSE, formula = y ~ x, color = 'black', size = 1.5) 


