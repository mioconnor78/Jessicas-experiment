#### code to accompany 2016Aug temporal analysis file at line 108


#------
## this is for checking for weird temp readings
temps4a <- temps3 %>% 
  unite(date_complete, Year, Month, Day, time, sep = "-") %>%
  mutate(date_formatted_h = ymd_h(date_complete)) %>% 
  filter(!is.na(date_formatted_h))
View(temps4a)


## looking at daily temp changes; ok I'm satisfied that the temperatures look ok.  
temps.daily <-
  temps4 %>%
  group_by(Tank, date_formatted) %>%
  summarise(., min(temp), max(temp), range = (max(temp)-min(temp)), mean(temp)) %>%
  arrange(as.numeric(date_formatted), as.numeric(Tank))

names(temps.daily) <- c("Tank", "date_formatted", "min", "max", "range", "mean")
View(temps.daily)
plot(temps.daily$range ~ temps.daily$mean)
summary(lm((temps.daily$range ~ temps.daily$mean)))

hist(temps.daily$range)
min(temps.daily$range)
max(temps.daily$range)
mean(temps.daily$range)
sd(temps.daily$range)


temp.plot <- ggplot(data = temps4a, aes(x = date_formatted_h, y = temp)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  geom_point() +
  facet_wrap(~Tank, as.table = TRUE) #+
#geom_vline(xintercept = 2012-07-08 00:00:00)

ggsave("temp.plot.png", device = "png", width = 14, height = 14)

## ok, I can't think of a reason to remove that July 8th blip. So I'm leaving it.

temp.plot <- ggplot(data = temps4a, aes(x = date_formatted_h, y = temp)) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(strip.background = element_rect(colour="white", fill="white")) +
  geom_point(data = subset(temps4a, week == "2")) +
  facet_wrap(~Tank, as.table = TRUE) +
  geom_hline(yintercept = 30)

#aes(group = as.character(Tank), color = as.character(Tank),
#aes(group = as.character(Tank), color = as.character(Tank))

#scale_colour_grey(start = 0, end = 0.6, name = "Tank", guide = "none")

#----

## old and somewhat mangled code for correcting for physical fluxes

## correct for water-air oxygen flux: http://www.waterontheweb.org/under/waterquality/oxygen.html
### adjusting for effects of temperature on physical exchange
C.star <- function(T) exp(7.7117 - 1.31403*log(T+45.93)) 
#- 0.035
# the -0.035 is an approximate adjustment for elevation
C.star(T) # yields the oxygen concentration expected at a given temperature (T in C).
## 7.7117
## 1.314

# height above sea level: 
h <- 0.2
# pressure
P <- 5.25 * log(1 - (h/44.3))
P.wv <- function(T) exp(11.8571 - (3840.70*(T+273)/(T+273)) - 216,961/(T+273)^2)
theta <- function(T) 0.000975 - (0.00001426*T) + (6.436x10^-8 * T^2)
C.p <- function(T) C.star*P*[((1-P.wv/P)(1-theta*P))/((1-P.wv)(1-theta))] 
