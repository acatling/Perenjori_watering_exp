###### Leaf traits analysis of WA 2019 Perenjori Water Experiment
## Alexandra Catling
## March and June 2021

library(tidyverse)
library(ggplot2)

#Theme for plotting
my_theme <- theme(axis.title.x = element_text(size = 16),
axis.title.y = element_text(size = 16),
axis.text = element_text(size = 16),
strip.text.x = element_text(size = 16),
legend.text = element_text(size = 14),
legend.title = element_blank())

## Watering and shade/sun treatment data
treatments <- read_csv("Data/treatments_meta.csv")

##Isotope data
isotopedataraw <- read_csv("Data/leaf_isotope_data.csv")
IDkey <- read_csv("Data/leaf_isotope_ID_key.csv")

isotopedata <- cbind(isotopedataraw, IDkey)
head(isotopedata)
#Transformation that John used to calculate delta from raw values (x):
# (-8-x)/(1+x/1000)
isotopedata <- isotopedata %>% mutate(delta = (-8-D13C)/(1+D13C/1000))

##Leaf area, SLA, LDMC -- leaf traits
leafdataraw <- read_csv("Data/wa_2020_leaf_traits.csv")
str(leafdataraw)
#Calculating LDMC (dry mass/fresh mass)
leafdata <- leafdataraw %>% filter(!is.na(Dry_mass_mg)) %>% mutate(LDMC = Dry_mass_mg/Fresh_mass_mg)
#Calculating SLA
leafdata <- leafdata %>% filter(!is.na(Leaf_area_cm2)) %>% mutate(SLA = Leaf_area_cm2/Dry_mass_mg)

#### ISOTOPE DATA
#Plotting by species
ggplot(isotopedata, aes(y = delta, x = Species))+
  geom_boxplot()+
  geom_point(color = "dodgerblue", position = (position_jitter(width = .1)))+
  theme_classic()+
  my_theme
#Merging with data on Cover
#Note that we don't have this data for HYGL (merged across all sites)
traitkey <- read_csv("Data/Trait_key.csv")
coverisotope <- merge(isotopedata, traitkey)
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

#Count of sample size by species
countall <- isotopedata %>% group_by(Species) %>% count()
#Merging with means and averages
meanD13C <- isotopedata %>% group_by(Species) %>%
  summarise(mean_D13C = mean(delta, na.rm= TRUE),
            sd_D13C = sd(delta, na.rm=TRUE))
summaryD13C <- merge(meanD13C, countall)

#Need to import information on sun vs. shade ****** ############
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
leafdatamean <- leafdata %>% group_by(Species, Collection_site, Collection_plot, Individual) %>%
  mutate(mean_SLA = mean(SLA, na.rm= TRUE),
            sd_SLA = sd(SLA, na.rm=TRUE))
leafdatamean <- leafdatamean %>% group_by(Species, Collection_site, Collection_plot, Individual) %>%
  mutate(mean_LDMC = mean(LDMC, na.rm= TRUE),
         sd_LDMC = sd(LDMC, na.rm=TRUE))
#Reducing to one row
leafsimple <- leafdatamean %>% group_by(Species, Collection_site, Collection_plot, 
                      Individual, mean_SLA, sd_SLA, mean_LDMC, sd_LDMC) %>%
  filter(row_number() == 1)

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

#Calculating one value for each species, based on average of individual averages
speciestraits <- leafsimple %>% group_by(Species) %>% summarise(SLA = mean(mean_SLA),
                                                                sd_SLA = sd(mean_SLA),
                                                                LDMC = mean(mean_LDMC),
                                                                sd_LDMC = mean(mean_LDMC))
#Is survival rate related to SLA, LDMC or D13?
alltraits <- merge(speciestraits, summaryD13C)
survivaltraits <- merge(survivalsp, alltraits)
ggplot(survivaltraits, aes(x = log(SLA), y = survival_rate))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
ggplot(survivaltraits, aes(x = log(LDMC), y = survival_rate))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
ggplot(survivaltraits, aes(x = mean_D13C, y = survival_rate))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
#Seeing how leaf traits are correlated
ggplot(survivaltraits, aes(x = log(SLA), y = mean_D13C))+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme
#By sun/shade
#Importing info on whether collection sites were in sun or shade
traitkey <- read_csv("Data/Trait_key.csv")
leafsimple2 <- merge(leafsimple, traitkey)
covertraits <- leafsimple2 %>% group_by(Species, Cover) %>% summarise(SLA = mean(mean_SLA),
                                                                sd_SLA = sd(mean_SLA),
                                                                LDMC = mean(mean_LDMC),
                                                                sd_LDMC = mean(mean_LDMC))
### I think I need to do the below instead ---
test <- leafsimple2 %>% group_by(Species, Collection_site, Collection_plot, Individual) %>% 
  summarise(SLA = mean(mean_SLA),
            sd_SLA = sd(mean_SLA),
            LDMC = mean(mean_LDMC),
            sd_LDMC = mean(mean_LDMC))
#Merge back in Cover info
test2 <- merge(test, traitkey)
#survivaltraitscover <- merge(survivalcover, covertraits)
ggplot(test2, aes(x = Cover, y = log(SLA)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species))
ggplot(test2, aes(x = Cover, y = log(LDMC)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species))

#####Vero - can look within plots 
vero <- leafsimple %>% filter(Species == "VERO")
veroplot <- vero %>% group_by(Collection_site, Collection_plot) %>% summarise(SLA = mean(mean_SLA),
                                                             sd_SLA = sd(mean_SLA),
                                                             LDMC = mean(mean_LDMC),
                                                             sd_LDMC = mean(mean_LDMC),
                                                             area = mean(Leaf_area_cm2),
                                                             sd_area = sd(Leaf_area_cm2))
#Merge back with info on treatments
verotreatments <- merge(veroplot, traitkey)
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
#Need to rename Collection_site and Collection_plot to match Site and Plot
#below isn't renaming them, creating a new column but hey, it works!
verotest <- verotreatments %>% mutate (Site = Collection_site, Plot = Collection_plot)
verostplot <- merge(veroplot, verotest)
str(verostplot)
verostplot$Site <- as.factor(verostplot$Site)
ggplot(verostplot, aes(x = Site, y = survival_rate))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verostplot, aes(x = Cover, y = survival_rate))+
  geom_boxplot()+
  geom_jitter()+
  theme_classic()+
  my_theme
ggplot(verostplot, aes(x = log(LDMC), y = survival_rate))+
  geom_jitter()+
  theme_classic()+
  my_theme
#Plotting traits against canopy cover as continuous variable
verocover <- merge(verotest, canopydata)
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
veroisotope <- isotopedata %>% filter(Species == "VERO") %>% mutate(Site = Collection_site,
                                                                    Plot = Collection_plot)
verocoverisotope <- merge(veroisotope, canopydata)
ggplot(verocoverisotope, aes(x = cc_percentage, y = delta))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "lm")+
  theme_classic()+
  my_theme

