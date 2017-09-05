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
  select(-mean.temp) 

sizes.t2$Tankn <- rep(c(1:9), 2)
sizes.t3 <- left_join(sizes.t, sizes.t2, by = c("Tank", "trophic.level")) 

sizes.t$invTi <-  1/((sizes.t$temp.wk + 273)*k) # average temp of the tank each week
sizes.t$invTT <-  1/((sizes.t$temp.Tmn + 273)*k) # average temp of the tank over all weeks
sizes.t$invTi <- as.numeric(as.character(sizes.t$invTi))
sizes.t$invTT <- as.numeric(as.character(sizes.t$invTT))

sizes.t2 <- sizes.t[sizes.t$week >= '4',]


---------------
## daphnia sizes
daph <- read.csv("./data/zooplankton size for correction no eggs_new.csv")
daph$Tank <- as.character(daph$tank)

## curious if the numbers of daphnia match...
## they don't, so for now use the data in the community size file
test2 <- daph %>%
  filter(stage == "mature") %>%
  group_by(., Tank, week) %>%
  summarise(., length(size.in.mm)) %>%
  arrange(as.numeric(week), as.numeric(Tank))

test <- sizes.t3 %>%
  mutate(taxon = species) %>%
  mutate(taxon = replace(taxon, species == "calanoid", "Copepod")) %>%
  mutate(taxon = replace(taxon, species == "cyclopoid", "Copepod")) %>%
  group_by(., Tank, week, taxon) %>%
  summarise(., length(size)) %>%
  arrange(as.numeric(week), as.numeric(Tank)) %>%
  filter(taxon == "Daphnia")
------------------

data.N <- as.tibble(data) %>%
  group_by(Tank, trophic.level, week) %>%
  summarise(., N = mean(total.zoo.abundance.liter)) %>%
  arrange(trophic.level, week) 

data$N1 <- 5*(data$total.zoo.abundance.liter) + 1
data$N <- 5*(data$total.zoo.abundance.liter)


# visualizing data --------------------------------------------------------

## how are the size data distributed?
sizes <- ggplot(data = sizes.t3, aes(x = log(size))) +
  geom_histogram(binwidth=.5, colour="black", fill="white") +
  facet_grid(week~Tank)

sizes # pretty close to log normal
ggsave("sizes.png")
## seeing that we only have good coverage for week 8, so let's just use that.
## because we have uneven samples among groups, i think we want regress the means against temperature, not the full samples...(right?)
data1 <- sizes.t3[(sizes.t3$week == "8"),] %>%
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

modSF <- lm(log(size) ~ 1 + trophic.level*I(invTT - mean(invTT)) + trophic.level*I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modS8 <- lm(log(size) ~ 1 + trophic.level + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modS7 <- lm(log(size) ~ 1 + trophic.level, data=data1, na.action=na.omit)
modS6 <- lm(log(size) ~ 1 + I(invTT - mean(invTT)), data=data1, na.action=na.omit)
modS5 <- lm(log(size) ~ 1, data=data1, na.action=na.omit)

model.sel(modS5, modS6, modS7, modS8, modSF)

## or, what if we use an averaged model: 
m.avgS <- model.avg(modS7, modS8)
coefficients(m.avgS)
confint(m.avgS)


S.PZ <- function(x) { coefficients(m.avgS)[1] + coefficients(m.avgS)[3]*x}
SvalsPZ <- S.PZ(data1[(data1$trophic.level=="PZ"),]$invTT - mean(data1$invTT))

S.PZN <- function(x) { (coefficients(m.avgS)[1] + coefficients(m.avgS)[2]) + coefficients(m.avgS)[3]*x}
SvalsPZN <- S.PZN(data1[(data1$trophic.level=="PZN"),]$invTT- mean(data1$invTT))

# because there is no trend with tempeature, maybe don't show this plot? unles we want to talk about the species composition shifts (daphnia really only at warm temps in PZN treatments)
Fig3A <- ggplot(data = data1, aes(x = -invTT, y = log(size))) + #, ymax = 1.2)
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(trophic.level~.) + ## this sets it up as facets
  geom_point(aes(group = taxon, shape = taxon, color = taxon), size = 2, alpha = 0.5) + 
  scale_alpha("Tankn", guide = "none") +
  scale_colour_grey(start = 0, end = 0.6, guide = "none") +
  scale_x_continuous("Temperature (1/kTi)", sec.axis = sec_axis(~((1/(k*-.))-273), name = "deg Celcius")) +
  #xlab("Temperature 1/kTi") +
  ylab("ln(size)") +
  geom_smooth(data = data1[(data1$trophic.level=="PZ"),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = log(size)),  size = .8, color = alpha("steelblue", 0.5)) +
  geom_smooth(data = data1[(data1$trophic.level=="PZN"),], method = "lm", se = FALSE, inherit.aes = FALSE, aes(x = -invTT, y = SvalsPZN),  size = .8, color = alpha("steelblue", 0.5))

Fig3A

ggsave("Fig3A sizes wk 8.png", device = "png", width =4, height = 3) 

## abundance analysis
require(car)
require(MASS)

## abundance - waiting for raw abundance data to describe the distribution. it's count data, so makes me think we would be using a Poisson or negbinom model.

qqp(data[(data$trophic.level != "P"),]$N1, "lnorm")
qqp(data[(data$trophic.level != "P"),]$N1, "norm")
nbinom <- fitdistr(data[(data$trophic.level != "P"),]$N1, "Negative Binomial")
qqp(data[(data$trophic.level != "P"),]$N1, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
hist(data[(data$trophic.level != "P"),]$N)
plot(data[(data$trophic.level != "P"),]$N ~ data[(data$trophic.level != "P"),]$Tank)

## i think we're going to go with the lognormal distribution
### might just ax the random int...or test for it: 
modNF <- lme(log(N1) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], method="ML", na.action=na.omit) 
modNa <- lm(log(N1) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data[(data$trophic.level != "P"),], na.action=na.omit)

anova(modNF, modNa)

modNF <- lm(log(N1) ~ 1 + I(invTi - invTT) + trophic.level + trophic.level*I(invTi - invTT) + I(invTT - mean(invTT)) + I(invTi - invTT)*I(invTT - mean(invTT)), data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN8 <- lm(log(N1) ~ 1 + trophic.level*I(invTi - invTT) + trophic.level*I(invTT - mean(invTT)), data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN7 <- lm(log(N1) ~ 1 + I(invTi - invTT) + trophic.level * I(invTT - mean(invTT)), data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN6 <- lm(log(N1) ~ 1 + trophic.level*I(invTi - invTT), data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN5 <- lm(log(N1) ~ 1 + I(invTi - invTT) + trophic.level, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN4 <- lm(log(N1) ~ 1 + I(invTi - invTT)*I(invTT - mean(invTT)), random = ~ 1 | Tank, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN3 <- lm(log(N1) ~ 1 + I(invTi - invTT) + I(invTT - mean(invTT)), data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN2 <- lm(log(N1) ~ 1 + I(invTi - invTT), data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN1 <- lm(log(N1) ~ 1 + trophic.level, data=data[(data$trophic.level != "P"),], na.action=na.omit)
modN0 <- lm(log(N1) ~ 1, data=data[(data$trophic.level != "P"),], na.action=na.omit)

model.sel(modN0, modN1, modN2, modN3, modN4, modN5, modN6, modN7, modN8, modNF)

## or, what if we use an averaged model: 
m.avgN <- model.avg(modN5, modN6)
confint(modNF)



N.plot <- ggplot(data=data[(data$trophic.level != "P"),], aes(x = -invTi, y = log(N1))) + #, ymin = -2, ymax = 6
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  facet_grid(trophic.level~.) + ## this sets it up as facets
  geom_point(aes(group = as.character(Tankn), color = as.character(Tankn), shape = as.factor(week), alpha = Tankn), size = 2) + 
  scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none") +
  scale_alpha("Tankn", guide = "none") +
  #theme(legend.position = c(0.88, 0.15), legend.text=element_text(size=6)) + 
  #scale_shape(name = "Week", guide = guide_legend(ncol = 2, size = 6)) +
  xlab("") + #xlab("Temperature 1/kTi") +
  ylab("ln(No. / 5 L)")

N.plot
ggsave("Nplot.png", device = "png", width = 4, height = 3) # save for appendix





### used a fixed effects model only because i'm not aware of methods for the mixed effects models on an poisson distribution
### following http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
### mean < 5, so...
library("mlmRev")
#library(lme4)
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


