### Alexandra Catling
### Traits analysis of demographic responses
# WA Perenjori 2020 Watering Experiment
#Created 11/06/2023, updated July 2023

## Have trait data for vero at the plot level and other species for sun/shade
# except hygl which was pooled sun and shade
# Interested in whether species-level traits explain demographic responses
# so only want species-level averages of data from sunny sites?

#### Loading packages and data ####
#Data imported from data preparation sheet
source("data_preparation.R")
#vitaldata has all datasets combined 
# germination, survival, seed production, neighbour info, traits, quad factors
# one row per subplot with seeds sown - NAs are very important since, e.g. survival info is only on germinated subplots

#Functions file
source("functions.R")

#Packages
library(ggplot2)
library(tidyverse)

###### Importing trait data - for data preparation ####
#Trait key says which collection sites align with which experimental sites (sun/shade)
traitkey <- read_csv("Data/Trait_key.csv")
## Isotope data
isotopedataraw <- read_csv("Data/leaf_isotope_data.csv")
IDkey <- read_csv("Data/leaf_isotope_ID_key.csv")
isotopedata <- cbind(isotopedataraw, IDkey)
#Transformation that John used to calculate delta from raw values (x):
# (-8-x)/(1+x/1000). This should probably be in functions
isotopedata <- isotopedata %>% mutate(delta = (-8-D13C)/(1+D13C/1000))
#join with trait key data
isotopedata <- left_join(isotopedata, traitkey)

#Count of sample size by species
countall <- isotopedata %>% group_by(Species) %>% count()
#Merging with means and averages
meanD13C <- isotopedata %>% group_by(Species) %>%
  summarise(mean_D13C = mean(delta, na.rm= TRUE),
            sd_D13C = sd(delta, na.rm=TRUE))
speciesD13C <- merge(meanD13C, countall)

## Leaf traits - area, SLA, LDMC
leafdataraw <- read_csv("Data/wa_2020_leaf_traits.csv")
#Calculating LDMC (dry mass/fresh mass)
leafdata <- leafdataraw %>% filter(!is.na(Dry_mass_mg)) %>% mutate(LDMC = Dry_mass_mg/Fresh_mass_mg)
#Calculating SLA
leafdata <- leafdata %>% filter(!is.na(Leaf_area_cm2)) %>% mutate(SLA = Leaf_area_cm2/Dry_mass_mg)

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

#Calculating one value for each species, based on average of individual averages
speciestraits <- leafsimple %>% group_by(Species) %>% summarise(SLA = mean(mean_SLA),
                                                                sd_SLA = sd(mean_SLA),
                                                                LDMC = mean(mean_LDMC),
                                                                sd_LDMC = mean(mean_LDMC))
#Merging with isotope species summaries
allspeciestraits <- merge(speciestraits, speciesD13C)

### Traits by sun/shade (note that hygl was pooled)
# #This gives a summary value for sun and shade traits per species
# leafsimple2 <- left_join(leafsimple, traitkey)
# covertraits <- leafsimple2 %>% group_by(Species, Cover) %>% summarise(SLA = mean(mean_SLA),
#                                                                       sd_SLA = sd(mean_SLA),
#                                                                       LDMC = mean(mean_LDMC),
#                                                                       sd_LDMC = mean(mean_LDMC))
### This gives info per individual on sun vs shade
# leafsimplecover <- leafsimple2 %>% group_by(Species, Collection_site, Collection_plot, Individual) %>% 
#   summarise(SLA = mean(mean_SLA),
#             sd_SLA = sd(mean_SLA),
#             LDMC = mean(mean_LDMC),
#             sd_LDMC = mean(mean_LDMC))
# #Merge back in Cover info
# leafsimplecover <- merge(leafsimplecover, traitkey)

##### Plot-level vero data (not using, just species-level) ####
# veroplot <- leafsimple %>% filter(Species == "VERO") %>% group_by(Collection_site, Collection_plot) %>% 
#   summarise(SLA = mean(mean_SLA),
#             sd_SLA = sd(mean_SLA),
#             LDMC = mean(mean_LDMC),
#             sd_LDMC = mean(mean_LDMC),
#             area = mean(Leaf_area_cm2),
#             sd_area = sd(Leaf_area_cm2))
# #Plotting traits against canopy cover as continuous variable
# veroplot <- veroplot %>% mutate(Site = Collection_site, Plot = Collection_plot)
# veroplot <- veroplot %>% ungroup() %>% select(!(c(Collection_site, Collection_plot)))
# #Add isotope data
# veroisotope <- isotopedata %>% filter(Species == "VERO") %>% mutate(Site = Collection_site,
#                                                                     Plot = Collection_plot)
# veroplot <- left_join(veroplot, verocoverisotope) %>% select(Site, Plot, Species, SLA, sd_SLA, LDMC, sd_LDMC, area, sd_area, delta)

##### Merging trait data with demographic data #####
allspeciestraits <- allspeciestraits %>% select(Species, SLA, sd_SLA, LDMC, sd_LDMC, mean_D13C, sd_D13C)
vitaldata <- left_join(vitaldata, allspeciestraits)

#And merging with popdata
popdata <- left_join(popdata, allspeciestraits)

#### Check distributions? up to here* ####
#Distributions of traits and PC1, have to use datanotscaled to log
#SLA better logged?
#LDMC still strange logged
#WUE/D13C better logged but still left-skewed

#with(vitaldata, pairs(~SLA + LDMC + mean_D13C + PC1, diag.panel = panel.hist))
#with(vitaldata, pairs(~log(SLA+2) + LDMC + log(mean_D13C+2) + PC1, diag.panel = panel.hist))
#Logging and standardising variables
vitaldata$log_SLA <-log(vitaldata$SLA)
vitaldata$log_D13C <-log(vitaldata$mean_D13C)

vitaldata$std_log_SLA <- scale(vitaldata$log_SLA, center = TRUE, scale = TRUE)[,1]
vitaldata$std_LDMC <- scale(vitaldata$LDMC, center = TRUE, scale = TRUE)[,1]
vitaldata$std_log_D13C <- scale(vitaldata$log_D13C, center = TRUE, scale = TRUE)[,1]

