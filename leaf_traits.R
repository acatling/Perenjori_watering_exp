###### Leaf traits analysis of WA 2019 Perenjori Water Experiment
## Alexandra Catling
## March and June 2021

library(tidyverse)
library(ggplot2)
#Data imported from data preparation sheet
source("data_preparation.R")
#Kept the plots here, all needed files are in data prep sheet

#### ISOTOPE DATA
#Plotting by species
ggplot(isotopedata, aes(y = delta, x = Species))+
  geom_boxplot()+
  geom_point(color = "dodgerblue", position = (position_jitter(width = .1)))+
  theme_classic()+
  my_theme
#Merging with data on Cover
#Note that we don't have this data for HYGL (merged across all sites)
ggplot(coverisotope, aes(y = delta, x = Species))+
  geom_boxplot()+
  geom_point(color = "dodgerblue", position = (position_jitter(width = .1)))+
  theme_classic()
ggplot(coverisotope, aes(x = Cover, y = delta))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species))
#Smaller sample sizes for above because I merged across individuals to have enough leaves (and expensive!)

#Models of isotope data
#Need to import information on sun vs. shade ****** 
model1 <- aov(D13C ~ Species, data = isotopedata)
summary(model1)
bartlett.test(D13C ~ Species, data = isotopedata)
TukeyHSD(model1)
#HYGL is significantly different to everything (ARCA, LARO, PEAI, PLDE, POLE, TRCY, TROR and VERO)
#LARO is significantly different to ARCA, HYGL, TRCY, TROR and VERO
#POLE is significantly different to ARCA, HYGL, TRCY and VERO

#### SLA and LDMC
#Need to summarise data, raw data is three leaves per individual (not independent)
# Need a mean per individual with standard error/sd
#Then need to look at distribution of data (probably needs to be logged)
ggplot(leafsimple, aes(x = Species, y = SLA))+
  geom_boxplot()+
  geom_point(color = "dodgerblue", position = (position_jitter(width = .1)), alpha = 0.5)+
  theme_classic()+
  my_theme

ggplot(leafsimple, aes(x = Species, y = log(LDMC)))+
  geom_boxplot()+
  geom_point(color = "dodgerblue", position = (position_jitter(width = .1)), alpha = 0.5)+
  theme_classic()+
  my_theme

#Is survival rate related to SLA, LDMC or D13?
#Haven't brought in survival info from other file.
# survivaltraits <- merge(survivalsp, alltraits)
# ggplot(survivaltraits, aes(x = log(SLA), y = survival_rate))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   theme_classic()+
#   my_theme
# ggplot(survivaltraits, aes(x = log(LDMC), y = survival_rate))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   theme_classic()+
#   my_theme
# ggplot(survivaltraits, aes(x = mean_D13C, y = survival_rate))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   theme_classic()+
#   my_theme

#Seeing how leaf traits are correlated
ggplot(alltraits, aes(x = log(SLA), y = mean_D13C))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
#By sun/shade
ggplot(leafsimplecover, aes(x = Cover, y = log(SLA)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species))
ggplot(leafsimplecover, aes(x = Cover, y = log(LDMC)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species))

#####Vero - can look within plots 
ggplot(verotreatments, aes(x = Cover, y = area))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verotreatments, aes(x = Treatment, y = area))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verotreatments, aes(x = Cover, y = log(SLA)))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verotreatments, aes(x = Treatment, y = log(SLA)))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verotreatments, aes(x = Cover, y = log(LDMC)))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verotreatments, aes(x = Treatment, y = log(LDMC)))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
#Only for vero, can look at how relationship between fecundity/survival and leaf traits 
#at the plot/site level
ggplot(verotreatments, aes(x = Collection_site, y = log(SLA)))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
#Survival data
# ggplot(verostplot, aes(x = Site, y = survival_rate))+
#   geom_boxplot()+
#   geom_jitter()+
#   theme_classic()+
#   my_theme
# ggplot(verostplot, aes(x = Cover, y = survival_rate))+
#   geom_boxplot()+
#   geom_jitter()+
#   theme_classic()+
#   my_theme
# ggplot(verostplot, aes(x = log(LDMC), y = survival_rate))+
#   geom_jitter()+
#   theme_classic()+
#   my_theme
#Plotting traits against canopy cover as continuous variable
ggplot(verocover, aes(x = cc_percentage, y = log(SLA)))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
ggplot(verocover, aes(x = cc_percentage, y = log(LDMC)))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
ggplot(verocover, aes(x = cc_percentage, y = area))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
#And isotope data
ggplot(verocoverisotope, aes(x = cc_percentage, y = delta))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme





