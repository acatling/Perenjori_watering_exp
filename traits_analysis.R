### Alexandra Catling
### Traits analysis of demographic responses
# WA Perenjori 2020 Watering Experiment
#Created 11/06/2023, updated July 2023

## Have trait data for vero at the plot level and other species for sun/shade
# except hygl which was pooled sun and shade
# Interested in whether species-level traits explain demographic responses
# pooling sun and shade traits

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
library(DHARMa)
library(glmmTMB)

#library(see) # problems installing this to use check_model

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
  summarise(D13C = mean(delta, na.rm= TRUE),
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

### Importing John's data on species seed mass and maximum height (and SLA) ###
johntraitsraw <- read_csv("Data_and_code_from_others/Dwyer_WA_species_list_with_basic_traits.csv")
##select my species
johntraits <- johntraitsraw[c(28,47,52,159,144,14,15,124),]
#assign species names
johntraits$Species <- 'placeholder'
johntraits$Species <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "TRCY", "TROR", "VERO")
johntraits <- johntraits %>% select(Species, john_sla = SLA, max_height, seed_mass)

#Merge with my trait dataframe
allspeciestraits <- left_join(allspeciestraits, johntraits)

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
allspeciestraits <- allspeciestraits %>% select(Species, SLA, sd_SLA, LDMC, sd_LDMC, D13C, sd_D13C, max_height, seed_mass)
vitaldata <- left_join(vitaldata, allspeciestraits)

#Making vital data no pole
vitaldatanopole <- vitaldata %>% filter(!(Species== "POLE"))

#And mergingtraits  with popdata
popdata <- left_join(popdata, allspeciestraits)

## And removing POLE, not enough data for population growth rates
popdatanopole <- popdata %>% filter(!(Species== "POLE"))

#### Checking distributions and standardising variables ####
#seed mass could be square rooted?
vitaldatanopole$std_SLA <- scale(vitaldatanopole$SLA, center = TRUE, scale = TRUE)[,1]
vitaldatanopole$std_LDMC <- scale(vitaldatanopole$LDMC, center = TRUE, scale = TRUE)[,1]
vitaldatanopole$std_D13C <- scale(vitaldatanopole$D13C, center = TRUE, scale = TRUE)[,1]
vitaldatanopole$std_height <- scale(vitaldatanopole$max_height, center = TRUE, scale = TRUE)[,1]
vitaldatanopole$std_seedmass <- scale(vitaldatanopole$seed_mass, center = TRUE, scale = TRUE)[,1]

popdatanopole$std_SLA <- scale(popdatanopole$SLA, center = TRUE, scale = TRUE)[,1]
popdatanopole$std_LDMC <- scale(popdatanopole$LDMC, center = TRUE, scale = TRUE)[,1]
popdatanopole$std_D13C <- scale(popdatanopole$D13C, center = TRUE, scale = TRUE)[,1]
popdatanopole$std_height <- scale(popdatanopole$max_height, center = TRUE, scale = TRUE)[,1]
popdatanopole$std_seedmass <- scale(popdatanopole$seed_mass, center = TRUE, scale = TRUE)[,1]

stdspeciestraitspop <- popdatanopole %>% select(Species, std_SLA, std_LDMC, std_D13C, std_height, std_seedmass) %>%
  group_by(Species) %>% filter(row_number() == 1) %>% select(-(Neighbours01))

stdspeciestraitsvital <- vitaldatanopole %>% select(Species, std_SLA, std_LDMC, std_D13C, std_height, std_seedmass) %>%
  group_by(Species) %>% filter(row_number() == 1)

#### Plot traits by species ####
ggplot(allspeciestraits, aes(x= Species, y = SLA))+
  geom_col()+
  theme_classic()
ggplot(stdspeciestraits, aes(x= Species, y = std_SLA))+
  geom_col()+
  theme_classic()
ggplot(allspeciestraits, aes(x= Species, y = LDMC))+
  geom_col()+
  theme_classic()
ggplot(stdspeciestraits, aes(x= Species, y = std_LDMC))+
  geom_col()+
  theme_classic()
ggplot(allspeciestraits, aes(x= Species, y = D13C))+
  geom_col()+
  theme_classic()
ggplot(stdspeciestraits, aes(x= Species, y = std_D13C))+
  geom_col()+
  theme_classic()
ggplot(allspeciestraits, aes(x= Species, y = max_height))+
  geom_col()+
  theme_classic()
ggplot(allspeciestraits, aes(x= Species, y = seed_mass))+
  geom_col()+
  theme_classic()

ggplot(vitaldatanopole, aes(x=Neighbours01, y=surv_to_produce_seeds))+
  geom_boxplot()+
  geom_jitter(alpha=0.1, height=0.05)+
  theme_classic()

#### 2023 - Traits modulating population growth rate responses ####
# ## Model that converges with all traits is full model with interactions without PC2
# 
# ### SLA
# #simple model with PC1 and SLA only
# popslapc1 <- lmer(log_lambda ~ std_PC1 + std_SLA + std_PC1:std_SLA + (1|Site/Plot) + (1+std_PC1|Species), popdata)
# popslapc1dharma <- simulateResiduals(popslapc1)
# plot(popslapc1dharma)
# #not normally distributed, lambda is left-skewed
# #check_model(popslapc1)
# summary(popslapc1)
# 
# #Full model minus interactions
# #Allowing SLA to modulate responses to Treatment, PC1 and neighbours with PC2 as a covariate
# popsla <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
#                      Treatment:std_SLA + std_PC1:std_SLA + Neighbours01:std_SLA +
#                      (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), popdatanopole)
# popsladharma <- simulateResiduals(popsla)
# plot(popsladharma)
# summary(popsla)
# # 
# # #Full model with interactions
# # #Allowing SLA to modulate responses to Treatment, PC1 and neighbours with PC2 as a covariate
# # popslabig <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
# #                   Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
# #                  Treatment:std_SLA + std_PC1:std_SLA + Neighbours01:std_SLA + 
# #                  (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), popdata)
# # popslabigdharma <- simulateResiduals(popslabig)
# # plot(popslabigdharma)
# # summary(popslabig)
# # 
# # ### LDMC
# # popLDMCpc1 <- lmer(log_lambda ~ std_PC1 + std_LDMC + std_PC1:std_LDMC + (1|Site/Plot) + (1+std_PC1|Species), popdata)
# # popLDMCpc1dharma <- simulateResiduals(popLDMCpc1)
# # plot(popLDMCpc1dharma)
# # #not normally distributed, lambda is left-skewed
# # #check_model(popLDMCpc1)
# # summary(popLDMCpc1)
# # 
# # #Full model minus interactions
# # #Allowing LDMC to modulate responses to Treatment, PC1 and neighbours with PC2 as a covariate
# # popLDMC <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
# #                   Treatment:std_LDMC + std_PC1:std_LDMC + Neighbours01:std_LDMC +
# #                   (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), popdata)
# # popLDMCdharma <- simulateResiduals(popLDMC)
# # plot(popLDMCdharma)
# # summary(popLDMC)
# 
# #Full model with interactions
# popLDMCbig <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
#                      Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
#                      Treatment:std_LDMC + std_PC1:std_LDMC + Neighbours01:std_LDMC +
#                      (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), popdatanopole)
# #model fails to converge
# popLDMCbigdharma <- simulateResiduals(popLDMCbig)
# plot(popLDMCbigdharma)
# summary(popLDMCbig)
# 
# # popLDMCbig2 <- glmmTMB(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
# #                      Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
# #                      Treatment:std_LDMC + std_PC1:std_LDMC + Neighbours01:std_LDMC + 
# #                      (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")), popdata)
# # #control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
# # #control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
# # 
# # ### D13C
# # popD13Cpc1 <- lmer(log_lambda ~ std_PC1 + std_D13C + std_PC1:std_D13C + (1|Site/Plot) + (1+std_PC1|Species), popdata)
# # popD13Cpc1dharma <- simulateResiduals(popD13Cpc1)
# # plot(popD13Cpc1dharma)
# # #not normally distributed, lambda is left-skewed
# # #check_model(popD13Cpc1)
# # summary(popD13Cpc1)
# # 
# # #Full model minus interactions
# # #Allowing D13C to modulate responses to Treatment, PC1 and neighbours with PC2 as a covariate
# # #model DOES converge if I remove PC2 from the random effect. 
# # popD13C <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
# #                   Treatment:std_D13C + std_PC1:std_D13C + Neighbours01:std_D13C + 
# #                   (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), popdata)
# # #model does not converge
# # popD13Cdharma <- simulateResiduals(popD13C)
# # plot(popD13Cdharma)
# # summary(popD13C)
# # 
# # popD13C2 <- glmmTMB(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
# #                   Treatment:std_D13C + std_PC1:std_D13C + Neighbours01:std_D13C + 
# #                   (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")), popdata)
# # #optimisers don't help
# # 
# #Full model with interactions
# popD13Cbig <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
#                      Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
#                      Treatment:std_D13C + std_PC1:std_D13C + Neighbours01:std_D13C +
#                      (1|Site/Plot) + (Treatment + std_PC1 + std_PC2 + Neighbours01|Species), popdatanopole)
# #model does not converge
# popD13Cbigdharma <- simulateResiduals(popD13Cbig)
# plot(popD13Cbigdharma)
# summary(popD13Cbig)

#### Full models with interactions WITHOUT PC2 ####
slamodel <- lmer(log_lambda ~ Treatment + std_PC1 + Neighbours01 +
                      Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
                      Treatment:std_SLA + std_PC1:std_SLA + Neighbours01:std_SLA + 
                      (1|Site/Plot) + (Treatment + std_PC1 + Neighbours01|Species), popdatanopole)
slamodeldharma <- simulateResiduals(slamodel)
plot(slamodeldharma)
#happy!
summary(slamodel)

### Full models with interactions WITHOUT PC2
ldmcmodel <- lmer(log_lambda ~ Treatment + std_PC1 + Neighbours01 +
                       Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
                       Treatment:std_LDMC + std_PC1:std_LDMC + Neighbours01:std_LDMC + 
                       (1|Site/Plot) + (Treatment + std_PC1 + Neighbours01|Species), popdatanopole)
ldmcmodeldharma <- simulateResiduals(ldmcmodel)
#not great
plot(ldmcmodeldharma)
summary(ldmcmodel)

### Full models with interactions WITHOUT PC2
d13cmodel <- lmer(log_lambda ~ Treatment + std_PC1 + Neighbours01 +
                       Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 +
                       Treatment:std_D13C + std_PC1:std_D13C + Neighbours01:std_D13C + 
                       (1|Site/Plot) + (Treatment + std_PC1 + Neighbours01|Species), popdatanopole)
d13cmodeldharma <- simulateResiduals(d13cmodel)
plot(d13cmodeldharma)
#not great
summary(d13cmodel)

#### 21/06/23 Working with John - simplified models #####
#Took out watering treatment - otherwise asking too much of the model
#with only 8 species/data points
slamodel2 <- lmer(log_lambda ~ std_PC1 + Neighbours01 + std_SLA +
                   std_PC1:Neighbours01 + std_PC1:std_SLA + Neighbours01:std_SLA + 
                   (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
slamodel2dharma <- simulateResiduals(slamodel2)
plot(slamodel2dharma)
summary(slamodel2)

#Does SLA explain the three way interaction?
slamodel2 <- lmer(log_lambda ~ std_PC1*Neighbours01*std_SLA +
                    (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
slamodel2dharma <- simulateResiduals(slamodel2)
plot(slamodel2dharma)
summary(slamodel2)

#ldmc
ldmcmodel2 <- lmer(log_lambda ~ std_PC1 + Neighbours01 + std_LDMC + 
                     std_PC1:Neighbours01 + std_PC1:std_LDMC + Neighbours01:std_LDMC + 
                     (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
ldmcmodel2dharma <- simulateResiduals(ldmcmodel2)
plot(ldmcmodel2dharma)
summary(ldmcmodel2)

#Does ldmc explain the three way interaction?
ldmcmodel3 <- lmer(log_lambda ~ std_PC1*Neighbours01*std_LDMC +
                     (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
ldmcmodel3dharma <- simulateResiduals(ldmcmodel3)
plot(ldmcmodel3dharma)
summary(ldmcmodel3)

#D13C
d13cmodel2 <- lmer(log_lambda ~ std_PC1 + Neighbours01 + std_D13C + 
                     std_PC1:Neighbours01 + std_PC1:std_D13C + Neighbours01:std_D13C + 
                     (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
d13cmodel2dharma <- simulateResiduals(d13cmodel2)
plot(d13cmodel2dharma)
summary(d13cmodel2)

#Does d13c explain the three way interaction?
d13cmodel2 <- lmer(log_lambda ~ std_PC1*Neighbours01*std_D13C +
                     (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
d13cmodel2dharma <- simulateResiduals(d13cmodel2)
plot(d13cmodel2dharma)
summary(d13cmodel2)

#seedmass
seedmassmodel2 <- lmer(log_lambda ~ std_PC1 + Neighbours01 + std_seedmass + 
                         std_PC1:Neighbours01 + std_PC1:std_seedmass + Neighbours01:std_seedmass + 
                         (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
seedmassmodel2dharma <- simulateResiduals(seedmassmodel2)
plot(seedmassmodel2dharma)
summary(seedmassmodel2)

#Does seedmass explain the three way interaction?
seedmassmodel2 <- lmer(log_lambda ~ std_PC1*Neighbours01*std_seedmass +
                         (1|Site/Plot) + (std_PC1 + Neighbours01|Species), popdatanopole)
seedmassmodel2dharma <- simulateResiduals(seedmassmodel2)
plot(seedmassmodel2dharma)
summary(seedmassmodel2)

#responses to the experiment (things we manipulated) - 
#watering and neighbour removal treatments
#ALWAYS NEED TRAIT AS MAIN EFFECT
slawater <- lmer(log_lambda ~ Treatment + Neighbours01 + std_SLA + 
                   Treatment:Neighbours01 + Treatment:std_SLA + Neighbours01:std_SLA + 
                   (1|Site/Plot) + (Treatment + Neighbours01|Species), popdatanopole)
slawaterdharma <- simulateResiduals(slawater)
plot(slawaterdharma)
summary(slawater)

#this is exactly the same - just different notation
#building a model with fixed effects of each variable and two-way interactions
#between all of them
# slawater <- lmer(log_lambda ~ (Treatment + Neighbours01 + std_SLA)^2  +
#                    (1|Site/Plot) + (Treatment + Neighbours01|Species), popdatanopole)

#seed mass
#larger seeds - thicker seedling stems - faster shoot elongation rates, larger cotyledons
#John's seeds - good year, many batches

#### Simplified population growth models - one trait, one treatment ####
### PC1
#not sla, ldmc, d13c, seed mass or max height
#very simple model asking does SLA modulate the response of PC1 to the environment?
# similarly, does your response to PC1 depend on your sla?
slamodel <- lmer(log_lambda ~ std_PC1*std_SLA + (1|Site/Plot) + (std_PC1|Species), popdatanopole)
summary(slamodel)
#ldmc
ldmcmodel <- lmer(log_lambda ~ std_PC1*std_LDMC + (1|Site/Plot) + (std_PC1|Species), popdatanopole)
summary(ldmcmodel)
#D13C
d13cmodel <- lmer(log_lambda ~ std_PC1*std_D13C + (1|Site/Plot) + (std_PC1|Species), popdatanopole)
summary(d13cmodel)
#max height
heightmodel <- lmer(log_lambda ~ std_PC1*std_height + (1|Site/Plot) + (std_PC1|Species), popdatanopole)
summary(heightmodel)
#seed mass
seedmassmodel <- lmer(log_lambda ~ std_PC1*std_seedmass + (1|Site/Plot) + (std_PC1|Species), popdatanopole)
summary(seedmassmodel)

###Treatment
#not sla, ldmc, d13c, height or seed mass.
slamodel <- lmer(log_lambda ~ Treatment*std_SLA + (1|Site/Plot) + (Treatment|Species), popdatanopole)
summary(slamodel)
#ldmc
ldmcmodel <- lmer(log_lambda ~ Treatment*std_LDMC + (1|Site/Plot) + (Treatment|Species), popdatanopole)
summary(ldmcmodel)
#D13C
d13cmodel <- lmer(log_lambda ~ Treatment*std_D13C + (1|Site/Plot) + (Treatment|Species), popdatanopole)
summary(d13cmodel)
#max height
heightmodel <- lmer(log_lambda ~ Treatment*std_height + (1|Site/Plot) + (Treatment|Species), popdatanopole)
summary(heightmodel)
#seed mass
seedmassmodel <- lmer(log_lambda ~ Treatment*std_seedmass + (1|Site/Plot) + (Treatment|Species), popdatanopole)
summary(seedmassmodel)

###Neighbour removal
#not sla, ldmc, d13c, max height or seed mass.
slamodel <- lmer(log_lambda ~ Neighbours01*std_SLA + (1|Site/Plot) + (Neighbours01|Species), popdatanopole)
summary(slamodel)
#ldmc
ldmcmodel <- lmer(log_lambda ~ Neighbours01*std_LDMC + (1|Site/Plot) + (Neighbours01|Species), popdatanopole)
summary(ldmcmodel)
#D13C
d13cmodel <- lmer(log_lambda ~ Neighbours01*std_D13C + (1|Site/Plot) + (Neighbours01|Species), popdatanopole)
summary(d13cmodel)
#max height
heightmodel <- lmer(log_lambda ~ Neighbours01*std_height + (1|Site/Plot) + (Neighbours01|Species), popdatanopole)
summary(heightmodel)
#seed mass
seedmassmodel <- lmer(log_lambda ~ Neighbours01*std_seedmass + (1|Site/Plot) + (Neighbours01|Species), popdatanopole)
summary(seedmassmodel)

##### Models for each trait for germination, survival and seed production ####
##SURVIVAL
#survival - neighbours
#tolerance of neighbours in terms of survival depends on SLA, LDMC and seed mass
# abundance neighbours - yes SLA, LDMC, no seed mass (marginal)
survslamodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_SLA + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                          family = binomial, vitaldatanopole)
survsladharma <- simulateResiduals(survslamodel)
plot(survsladharma)
summary(survslamodel)
r.squaredGLMM(survslamodel)

survslamodel <- glmer(surv_to_produce_seeds ~ Neighbours01*std_SLA + (1|Site/Plot) + (Neighbours01|Species), 
                      family = binomial, vitaldatanopole)
survsladharma <- simulateResiduals(survslamodel)
plot(survsladharma)
summary(survslamodel)
r.squaredGLMM(survslamodel)

survLDMCmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_LDMC + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                       family = binomial, vitaldatanopole)
survLDMCdharma <- simulateResiduals(survLDMCmodel)
plot(survLDMCdharma)
summary(survLDMCmodel)
r.squaredGLMM(survLDMCmodel)

survLDMCmodel <- glmer(surv_to_produce_seeds ~ Neighbours01*std_LDMC + (1|Site/Plot) + (Neighbours01|Species), 
                       family = binomial, vitaldatanopole)
survLDMCdharma <- simulateResiduals(survLDMCmodel)
plot(survLDMCdharma)
summary(survLDMCmodel)
r.squaredGLMM(survLDMCmodel)

survD13Cmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_D13C + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                       family = binomial, vitaldatanopole)
survD13Cdharma <- simulateResiduals(survD13Cmodel)
plot(survD13Cdharma)
summary(survD13Cmodel)
r.squaredGLMM(survD13Cmodel)

survD13Cmodel <- glmer(surv_to_produce_seeds ~ Neighbours01*std_D13C + (1|Site/Plot) + (Neighbours01|Species), 
                       family = binomial, vitaldatanopole)
survD13Cdharma <- simulateResiduals(survD13Cmodel)
plot(survD13Cdharma)
summary(survD13Cmodel)

survheightmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_height + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                         family = binomial, vitaldatanopole)
survheightdharma <- simulateResiduals(survheightmodel)
plot(survheightdharma)
summary(survheightmodel)
survheightmodel <- glmer(surv_to_produce_seeds ~ Neighbours01*std_height + (1|Site/Plot) + (Neighbours01|Species), 
                       family = binomial, vitaldatanopole)
survheightdharma <- simulateResiduals(survheightmodel)
plot(survheightdharma)
summary(survheightmodel)

survseedmassmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_seedmass + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                           family = binomial, vitaldatanopole)
survseedmassdharma <- simulateResiduals(survseedmassmodel)
plot(survseedmassdharma)
summary(survseedmassmodel)
r.squaredGLMM(survseedmassmodel)

survseedmassmodel <- glmer(surv_to_produce_seeds ~ Neighbours01*std_seedmass + (1|Site/Plot) + (Neighbours01|Species), 
                           family = binomial, vitaldatanopole)
survseedmassdharma <- simulateResiduals(survseedmassmodel)
plot(survseedmassdharma)
summary(survseedmassmodel)
r.squaredGLMM(survseedmassmodel)

### survival - watering treatment
#response to watering not dependent on sla, ldmc, D13C, height or seed mass
survslamodel <- glmer(surv_to_produce_seeds ~ Treatment*std_SLA + (1|Site/Plot) + (Treatment|Species), 
                      family = binomial, vitaldatanopole)
survsladharma <- simulateResiduals(survslamodel)
plot(survsladharma)
summary(survslamodel)

survLDMCmodel <- glmer(surv_to_produce_seeds ~ Treatment*std_LDMC + (1|Site/Plot) + (Treatment|Species), 
                       family = binomial, vitaldatanopole)
survLDMCdharma <- simulateResiduals(survLDMCmodel)
plot(survLDMCdharma)
summary(survLDMCmodel)

survD13Cmodel <- glmer(surv_to_produce_seeds ~ Treatment*std_D13C + (1|Site/Plot) + (Treatment|Species), 
                       family = binomial, vitaldatanopole)
survD13Cdharma <- simulateResiduals(survD13Cmodel)
plot(survD13Cdharma)
summary(survD13Cmodel)

survheightmodel <- glmer(surv_to_produce_seeds ~ Treatment*std_height + (1|Site/Plot) + (Treatment|Species), 
                       family = binomial, vitaldatanopole)
summary(survheightmodel)

survseedmassmodel <- glmer(surv_to_produce_seeds ~ Treatment*std_seedmass + (1|Site/Plot) + (Treatment|Species), 
                         family = binomial, vitaldatanopole)
summary(survseedmassmodel)

### survival - PC1
#no response to PC1 explained by sla, ldmc, D13C, height or seed mass
survslamodel <- glmer(surv_to_produce_seeds ~ std_PC1*std_SLA + (1|Site/Plot) + (std_PC1|Species), 
                      family = binomial, vitaldatanopole)
survsladharma <- simulateResiduals(survslamodel)
plot(survsladharma)
summary(survslamodel)

survLDMCmodel <- glmer(surv_to_produce_seeds ~ std_PC1*std_LDMC + (1|Site/Plot) + (std_PC1|Species), 
                       family = binomial, vitaldatanopole)
survLDMCdharma <- simulateResiduals(survLDMCmodel)
plot(survLDMCdharma)
summary(survLDMCmodel)

survD13Cmodel <- glmer(surv_to_produce_seeds ~ std_PC1*std_D13C + (1|Site/Plot) + (std_PC1|Species), 
                       family = binomial, vitaldatanopole)
survD13Cdharma <- simulateResiduals(survD13Cmodel)
plot(survD13Cdharma)
summary(survD13Cmodel)

survheightmodel <- glmer(surv_to_produce_seeds ~ std_PC1*std_height + (1|Site/Plot) + (std_PC1|Species), 
                       family = binomial, vitaldatanopole)
summary(survheightmodel)

survseedmassmodel <- glmer(surv_to_produce_seeds ~ std_PC1*std_seedmass + (1|Site/Plot) + (std_PC1|Species), 
                         family = binomial, vitaldatanopole)
summary(survseedmassmodel)

##FECUNDITY
#### seed production - neighbours
#running into convergence issues with these models
#tolerance of neighbours in terms of seed production depends on SLA only
#continous neighbours - not SLA
seedslamodel <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund*std_SLA + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                        family = nbinom2, vitaldatanopole)
summary(seedslamodel)
seedslamodel <- glmmTMB(No_viable_seeds_grouped ~ Neighbours01*std_SLA + (1|Site/Plot) + (Neighbours01|Species), 
                        family = nbinom2, vitaldatanopole)
# seedslamodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01*std_SLA + (1|Site/Plot) + (Neighbours01|Species), 
#                         vitaldatanopole)
seedsladharma <- simulateResiduals(seedslamodel)
plot(seedsladharma)
summary(seedslamodel)

seedLDMCmodel <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund*std_LDMC + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                        family = nbinom2, vitaldatanopole)
summary(seedLDMCmodel)
seedLDMCmodel <- glmmTMB(No_viable_seeds_grouped ~ Neighbours01*std_LDMC + (1|Site/Plot) + (Neighbours01|Species), 
                         family = nbinom2, vitaldatanopole)
seedLDMCdharma <- simulateResiduals(seedLDMCmodel)
plot(seedLDMCdharma)
summary(seedLDMCmodel)

seedD13Cmodel <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund*std_D13C + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                         family = nbinom2, vitaldatanopole)
summary(seedD13Cmodel)
seedD13Cmodel <- glmmTMB(No_viable_seeds_grouped ~ Neighbours01*std_D13C + (1|Site/Plot) + (Neighbours01|Species), 
                         family = nbinom2, vitaldatanopole)
seedD13Cdharma <- simulateResiduals(seedD13Cmodel)
plot(seedD13Cdharma)
summary(seedD13Cmodel)

seedheightmodel <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund*std_height + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                         family = nbinom2, vitaldatanopole)
summary(seedheightmodel)
seedheightmodel <- glmmTMB(No_viable_seeds_grouped ~ Neighbours01*std_height + (1|Site/Plot) + (Neighbours01|Species), 
                         family = nbinom2, vitaldatanopole)
summary(seedheightmodel)

seedseedmassmodel <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund*std_seedmass + (1|Site/Plot) + (std_logp1_totalabund|Species), 
                           family = nbinom2, vitaldatanopole)
summary(seedseedmassmodel)
seedseedmassmodel <- glmmTMB(No_viable_seeds_grouped ~ Neighbours01*std_seedmass + (1|Site/Plot) + (Neighbours01|Species), 
                           family = nbinom2, vitaldatanopole)
summary(seedseedmassmodel)

### seed production - watering treatment
#response to watering not dependent on sla, ldmc, D13C, seed height or seed mass
seedslamodel <- glmmTMB(No_viable_seeds_grouped ~ Treatment*std_SLA + (1|Site/Plot) + (Treatment|Species), 
                        family = nbinom2, vitaldatanopole)
seedsladharma <- simulateResiduals(seedslamodel)
plot(seedsladharma)
summary(seedslamodel)

seedLDMCmodel <- glmmTMB(No_viable_seeds_grouped ~ Treatment*std_LDMC + (1|Site/Plot) + (Treatment|Species), 
                         family = nbinom2, vitaldatanopole)
seedLDMCdharma <- simulateResiduals(seedLDMCmodel)
plot(seedLDMCdharma)
summary(seedLDMCmodel)

seedD13Cmodel <- glmmTMB(No_viable_seeds_grouped ~ Treatment*std_D13C + (1|Site/Plot) + (Treatment|Species), 
                         family = nbinom2, vitaldatanopole)
seedD13Cdharma <- simulateResiduals(seedD13Cmodel)
plot(seedD13Cdharma)
summary(seedD13Cmodel)

seedheightmodel <- glmmTMB(No_viable_seeds_grouped ~ Treatment*std_height + (1|Site/Plot) + (Treatment|Species), 
                         family = nbinom2, vitaldatanopole)
summary(seedheightmodel)

seedseedmassmodel <- glmmTMB(No_viable_seeds_grouped ~ Treatment*std_seedmass + (1|Site/Plot) + (Treatment|Species), 
                           family = nbinom2, vitaldatanopole)
summary(seedseedmassmodel)

### seed production - PC1
#no response to PC1 explained by sla, ldmc, D13C, height or seed mass
seedslamodel <- glmmTMB(No_viable_seeds_grouped ~ std_PC1*std_SLA + (1|Site/Plot) + (std_PC1|Species), 
                        family = nbinom2, vitaldatanopole)
seedsladharma <- simulateResiduals(seedslamodel)
plot(seedsladharma)
summary(seedslamodel)

seedLDMCmodel <- glmmTMB(No_viable_seeds_grouped ~ std_PC1*std_LDMC + (1|Site/Plot) + (std_PC1|Species), 
                         family = nbinom2, vitaldatanopole)
seedLDMCdharma <- simulateResiduals(seedLDMCmodel)
plot(seedLDMCdharma)
summary(seedLDMCmodel)

seedD13Cmodel <- glmmTMB(No_viable_seeds_grouped ~ std_PC1*std_D13C + (1|Site/Plot) + (std_PC1|Species), 
                         family = nbinom2, vitaldatanopole)
seedD13Cdharma <- simulateResiduals(seedD13Cmodel)
plot(seedD13Cdharma)
summary(seedD13Cmodel)

seedheightmodel <- glmmTMB(No_viable_seeds_grouped ~ std_PC1*std_height + (1|Site/Plot) + (std_PC1|Species), 
                         family = nbinom2, vitaldatanopole)
summary(seedheightmodel)

seedseedmassmodel <- glmmTMB(No_viable_seeds_grouped ~ std_PC1*std_seedmass + (1|Site/Plot) + (std_PC1|Species), 
                           family = nbinom2, vitaldatanopole)
summary(seedseedmassmodel)

###### Plotting #####

#traits_for_model needs to be STANDARDISED SLA
#ARCA: -0.27 #HYGL -1.4
test <- popdatanopole %>% select(Species, std_SLA, std_LDMC, std_D13C) %>%
  group_by(Species) %>% filter(row_number() == 1) %>% select(-(Neighbours01))

traits_for_model <- speciestraits[match(rownames(coef(slamodel)$Species), speciestraits$Species),c(1:2)]
varfix <- vcov(slamodel)[2,2]
#vcov has four rows - fixed effects - should be grabbing PC1:PC1?? check notes
re <- ranef(slamodel,condVar=TRUE)
#now sure what this below line is doing - not sure if grabbing the right cell
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(slamodel)$Species$std_PC1 + coef(slamodel)$Species$`std_PC1:std_SLA`*traits_for_model$SLA)

## Always need to multiple by STANDARDISED SLA value and plot against 
# STANDARDISED SLA value
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdatanopole$std_SLA))
with(slopedata, plot(slope_PC1 ~ SLA, ylab= "Slope of pop growth-PC1 relationship", xlab="SLA (standardised)", ylim=c(-0.6,0.6), pch = 19, col="grey70", cex.lab = 1, cex.axis = 2, cex = 2))
with(slopedata, arrows(SLA, slope_PC1+slope.SEs, SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 3, cex = 3))
abline(h=0, lty=3, lwd=3)
abline(a=-0.15980, b=-0.13497, lwd=2)

#### 28/06/23 attempt LDMC slope surv~NA
survLDMCmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_LDMC + (1|Site/Plot) + (std_logp1_totalabund|Species),
                                              family = binomial, vitaldatanopole)
traits_for_model <- stdspeciestraitsvital[match(rownames(coef(survLDMCmodel)$Species), stdspeciestraitsvital$Species),c(1,3)]
#Need standardised LDMC values
#not sure why below line isn't working
#slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(survLDMCmodel)$Species$std_logp1_totalabund + coef(survLDMCmodel)$Species$`std_logp1_totalabund:std_LDMC`*traits_for_model$std_LDMC)
slopedata <- traits_for_model
slopedata$slope_NA <- coef(survLDMCmodel)$Species$std_logp1_totalabund + coef(survLDMCmodel)$Species$`std_logp1_totalabund:std_LDMC`*traits_for_model$std_LDMC
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(vitaldatanopole$std_LDMC))
with(slopedata, plot(slope_NA ~ std_LDMC, ylab= "Slope of survival-neighbour abundance relationship", xlab="LDMC (standardised)", pch = 19, col="grey70", cex.lab = 1, cex.axis = 2, cex = 2))
#with(slopedata, arrows(std_LDMC, slope_NA+slope.SEs, std_LDMC, slope_NA-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 3, cex = 3))
#abline a (intercept) is the estimate of std_NA
#abline b (slope) is estimate of std_NA:std_LDMC 
summary(survLDMCmodel)
abline(a=0, b=0, lwd=2)
abline(a=-0.33149, b=0.15776, lwd=2)


# data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdatanopole$std_SLA))
# pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(slamodel)[c(6,9)]
# sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(slamodel)[c(6,9), c(6,9)], as.matrix(data.for.sla.slope.regression))))
# with(slopedata, plot(slope_PC1 ~ SLA, ylab= "Slope of pop growth-PC1 relationship", xlab="SLA", ylim=c(-0.6,0.6), pch = 19, col="grey70", cex.lab = 1, cex.axis = 2, cex = 3))
# with(slopedata, arrows(SLA, slope_PC1+slope.SEs, SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
# plot.CI.func(x.for.plot=seq.func(popdatanopole$std_SLA), pred=pred.sla.slope.regression,
#              upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
#              lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
#              env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
# abline(h=0, lty=3, lwd=3, cex = 3)

#### 2023 - plotting slope of pop-PC1 and pop-NA relationships against traits ####
#John: this bit grabs the random slopes for each species as well as their associated SEs
#NOT SURE THESE ARE THE RIGHT VALUES?** I changed this to [8,8] to be SLA column and row but
#not sure that the random effects part is correct. What does [2,2,] pull out??
#I think [2,2,] pulls out the random slope values for PC1 which should be correct as [2,2,] since it is 
# the first random slope fitted in the model and [1,1,] is the random intercept?

#Calculating the slope values for the response of fecundity to PC1 as modulated by SLA
# slopes = the slope responses of PC2 for each species + the interaction of PC2:SLA*each species' SLA value
#Unique value for PC2 slopes per species (what's left over), one value for PC2:SLA for all species, unique values of SLA per species
#This changes the position of each species: PC2:SLA - which is one it's one value
# This is all because the fixed effects have explained some of the response already

### Plotting response of fecundity to PC1 as a function of SLA
#c(6,9) is pulling out std_PC1 and std_PC1*std_log_SLA

#### 2023
# Why are all of these PC1:SLA values the same? All of the interactions
# have the same values for all species -- because they haven't been multiplied 
#by their invididual trait values yet! modulated by SLA
# Each species has different intercepts but the same slopes? No, different slopes
#coef(slamodel)$Species[,14]

## We need std_SLA rather than just SLA to match
#Does it make sense to have sla standardised in the model?
# Need to plot against not standardised sla...

#All have the same amount of error?
# This plot doesn't line up with summary of the model though, which suggests
# there is almost a significant relationship in PC1:SLA 

summary(slamodel)
traits_for_model <- speciestraits[match(rownames(coef(slamodel)$Species), speciestraits$Species),c(1:2)]
#below is standardised
# traits_for_model <- popdatanopole %>% select(Species, std_SLA) %>% 
#   group_by(Species) %>% filter(row_number() == 1) %>% select(Species, std_SLA)
varfix <- vcov(slamodel)[2,2]
re <- ranef(slamodel,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(slamodel)$Species[,4] + coef(slamodel)$Species[,14]*traits_for_model$SLA)
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdatanopole$std_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(slamodel)[c(6,9)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(slamodel)[c(6,9), c(6,9)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ SLA, ylab= "Slope of pop growth-PC1 relationship", xlab="SLA", ylim=c(-0.6,0.6), pch = 19, col="grey70", cex.lab = 1, cex.axis = 2, cex = 3))
with(slopedata, arrows(SLA, slope_PC1+slope.SEs, SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
plot.CI.func(x.for.plot=seq.func(popdatanopole$std_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
abline(h=0, lty=3, lwd=3, cex = 3)

### LDMC
summary(ldmcmodel)
traits_for_model <- speciestraits[match(rownames(coef(ldmcmodel)$Species), speciestraits$Species),c(1:4)]
traits_for_model <- traits_for_model %>% select(Species, LDMC)
varfix <- vcov(ldmcmodel)[2,2]
re <- ranef(ldmcmodel,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(ldmcmodel)$Species[,4] + coef(ldmcmodel)$Species[,14]*traits_for_model$LDMC)
data.for.ldmc.slope.regression<-cbind(int=1, x=seq.func(popdatanopole$std_LDMC))
pred.ldmc.slope.regression<- data.for.ldmc.slope.regression%*%fixef(ldmcmodel)[c(6,9)]
ldmc.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.ldmc.slope.regression) %*% tcrossprod(vcov(ldmcmodel)[c(6,9), c(6,9)], as.matrix(data.for.ldmc.slope.regression))))
with(slopedata, plot(slope_PC1 ~ LDMC, ylab= "Slope of pop growth-PC1 relationship", xlab="ldmc", ylim=c(-0.6,0.6), pch = 19, col="grey70", cex.lab = 1, cex.axis = 2, cex = 3))
with(slopedata, arrows(LDMC, slope_PC1+slope.SEs, LDMC, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
plot.CI.func(x.for.plot=seq.func(popdatanopole$std_LDMC), pred=pred.ldmc.slope.regression,
             upper=pred.ldmc.slope.regression + 1.96*ldmc.slope.regression.SEs,
             lower=pred.ldmc.slope.regression - 1.96*ldmc.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
abline(h=0, lty=3, lwd=3, cex = 3)

#### Plotting survival probability ~ neighbours by SLA, LDMC and seed mass ####

#col=ifelse(vitaldatanopole$std_SLA<0, alpha("#CC79A7", 0.1), alpha("#0072B2", 0.1))

dev.off()
pdf("Output/Figures/surv~nbh_by_trait.pdf", width=16, height=7)
par(mfrow=c(1,2), oma = c(6, 6, 1, 6), mar =c(2,2,1,1))

#sla
x_to_plot_high <- seq.func(vitaldatanopole$std_logp1_totalabund[vitaldatanopole$std_SLA>0])
x_to_plot <- seq.func(vitaldatanopole$std_logp1_totalabund)
x_to_plot_low <- seq.func(vitaldatanopole$std_logp1_totalabund[vitaldatanopole$std_SLA<0])
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,3), pch=19, col=alpha("grey60", 0.1), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, vitaldatanopole)
survslamodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_SLA + (1|Site/Plot) + (std_logp1_totalabund|Species),
                      family = binomial, vitaldatanopole)
model <- survslamodel
#ambient - black
# plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
# plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#high sla - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 2, x_to_plot_high*2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_high, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#low sla - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, -1.4, x_to_plot_low*-1.4), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_low, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#LDMC
x_to_plot_high <- seq.func(vitaldatanopole$std_logp1_totalabund[vitaldatanopole$std_LDMC>0])
x_to_plot <- seq.func(vitaldatanopole$std_logp1_totalabund)
x_to_plot_low <- seq.func(vitaldatanopole$std_logp1_totalabund[vitaldatanopole$std_LDMC<0])
plot(surv_to_produce_seeds~ std_logp1_totalabund, xlim=c(-1,3), pch=19, col=alpha("grey60", 0.1), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, vitaldatanopole)
survLDMCmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_LDMC + (1|Site/Plot) + (std_logp1_totalabund|Species),
                       family = binomial, vitaldatanopole)
model <- survLDMCmodel
#ambient - black
#plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
#plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#high LDMC - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1.7, x_to_plot_high*1.7), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_high, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#low LDMC - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, -1.1, x_to_plot_low*-1.1), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_low, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

dev.off()

ggplot(vitaldatanopole, aes(x = Neighbours01, y = surv_to_produce_seeds))+
  geom_boxplot()+
  #geom_jitter()+
  theme_classic()

#############
#seedmass
x_to_plot_high <- seq.func(vitaldatanopole$std_logp1_totalabund[vitaldatanopole$std_seedmass>0])
x_to_plot <- seq.func(vitaldatanopole$std_logp1_totalabund)
x_to_plot_low <- seq.func(vitaldatanopole$std_logp1_totalabund[vitaldatanopole$std_seedmass<0])
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,3), pch=19, col=alpha("grey60", 0.1), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, vitaldatanopole)
survseedmassmodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund*std_seedmass + (1|Site/Plot) + (std_logp1_totalabund|Species),
                           family = binomial, vitaldatanopole)
model <- survseedmassmodel
#ambient - black
#plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
#plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#high seedmass - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1.4, x_to_plot_high*1.4), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_high, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#low seedmass - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, -1.6, x_to_plot_low*-1.6), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_low, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

dev.off()

### Plot SLA seed production ####
seedslamodel <- glmmTMB(No_viable_seeds_grouped ~ Neighbours01*std_SLA + (1|Site/Plot) + (std_PC1|Species), 
                        family = nbinom2, vitaldatanopole)
seedsladharma <- simulateResiduals(seedslamodel)
plot(seedsladharma)
summary(seedslamodel)


##Neighbours 01 survival
#### Old code below here - traits explaining demographic responses ####
#Plotting PC1 simply
species_slopes<-coef(lambdaD13C)$Species[,4] + coef(lambdaD13C)$Species[,9]*speciesD13C$std_log_D13C
plot(species_slopes ~ speciesD13C$std_log_D13C)
curve(cbind(1,x)%*%fixef(lambdaD13C)[c(4,9)], add=T)
abline(h=0, lty=3)

# #Trait data needed for plotting
# #Creating a dataframe with SLA data - need the as.data.frame part!
# speciessla <- seedmodeldata %>% select("Species", "std_log_SLA") %>% group_by(Species) %>% filter(row_number() == 1)
# speciessla <- as.data.frame(speciessla)
# 
# #Creating a dataframe with D13C data
# speciesD13C <- seedmodeldata %>% select("Species", "std_log_D13C") %>% group_by(Species) %>% filter(row_number() == 1)
# speciesD13C <- as.data.frame(speciesD13C)

# ### John's work 07/11/2021 plotting simple model sla modulating response of fecundity to PC2 ###
# seedsla2 <- glmer.nb(seeds_percent ~ std_PC2 + std_log_SLA + std_PC2:std_log_SLA +
#                        (1|Site/Plot) + (1+std_PC2|Species), family = nbinom2, seedmodeldata)
# summary(seedsla2)
# #My interpretation: the slopes are the slope responses of PC2 for each species plus
# #the interaction of PC2:SLA*each species' SLA value.
# #Unique value for PC2 slopes, one value for PC2:SLA, unique values of SLA
# species_slopes<-coef(seedsla2)$Species[,2] + coef(seedsla2)$Species[,4]*speciessla$std_log_SLA
# 
# plot(species_slopes ~ speciessla$std_log_SLA)
# curve(cbind(1,x)%*%fixef(seedsla2)[c(2,4)], add=T)
# abline(h=0, lty=3)
# 
# ##Does SLA modulate the responses of survival to PC1, PC2, neighbour abundance or treatment?
# ### Survival models ###
# #SLA
# survsla <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_SLA + 
#                    std_PC1:std_log_SLA + std_logp1_totalabund:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
#                    (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, vitaldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# survsladharma <- simulateResiduals(survsla)
# plot(survsladharma)
# #good
# summary(survsla)
# #No interactions significant
# 
# #LDMC
# survldmc <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_LDMC + 
#                     std_PC1:std_LDMC + std_logp1_totalabund:std_LDMC + std_PC2:std_LDMC + Treatment:std_LDMC +
#                     (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, vitaldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# survldmcdharma <- simulateResiduals(survldmc)
# #good but outlier test significant?
# plot(survldmcdharma)
# summary(survldmc)
# #No interactions significant
# 
# #D13C
# survD13C <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_D13C + 
#                     std_PC1:std_log_D13C + std_logp1_totalabund:std_log_D13C + std_PC2:std_log_D13C + Treatment:std_log_D13C +
#                     (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, vitaldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# survD13Cdharma <- simulateResiduals(survD13C)
# plot(survD13Cdharma)
# #good
# summary(survD13C)
# #No interactions significant
# 
# ### Fecundity models ###
# #SLA
# #control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
# #Converges with glmer.nb instead of glmmTMB!
# seedsla <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_SLA +
#                       std_PC1:std_log_SLA + std_logp1_totalabund:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
#                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, seedmodeldata)
# ## Saving the model to load later
# save(seedsla, file = "Output/Models/seedslamodel.RData")
# seedsladharma <- simulateResiduals(seedsla)
# plot(seedsladharma)
# #okay?
# summary(seedsla)
# #PC1 and PC2 interactions significant
# 
# #LDMC
# #Will only converge with optimiser
# seedldmc <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + std_LDMC + 
#                        std_PC1:std_LDMC + std_logp1_totalabund:std_LDMC + std_PC2:std_LDMC + Treatment:std_LDMC +
#                        (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, 
#                      seedmodeldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# ## Saving the model to load later
# save(seedldmc, file = "Output/Models/seedldmcmodel.RData")
# seedldmcdharma <- simulateResiduals(seedldmc)
# plot(seedldmcdharma)
# #good
# summary(seedldmc)
# #no significant interactions
# 
# #D13C
# seedD13C <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + std_log_D13C + 
#                        std_PC1:std_log_D13C + std_logp1_totalabund:std_log_D13C + std_PC2:std_log_D13C + Treatment:std_log_D13C +
#                        (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, seedmodeldata,
#                      control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# ## Saving the model to load later
# save(seedD13C, file = "Output/Models/seedD13Cmodel.RData")
# seedD13Cdharma <- simulateResiduals(seedD13C)
# plot(seedD13Cdharma)
# #Okay, KS test significant
# summary(seedD13C)
# #No significant interactions
# 
# ### Lambda models ###
# #Trying simple model first - SLA modulating reponse to PC1
# lambdaslasimple <- lmer(log_lambda ~ std_PC1 + std_SLA + std_PC1:std_SLA + (1|Site/Plot) + (std_PC1|Species), popdata)
# lambdaslasimpledharma <- simulateResiduals(lambdaslasimple)
# plot(lambdaslasimpledharma)
# #Okay, KS test significant
# summary(lambdaslasimple)
# 
# #Simple model with PC1 and Treatment
# lambdaslasimple2 <- lmer(log_lambda ~ std_PC1 + Treatment + std_SLA + std_PC1:std_SLA + Treatment:std_SLA + (1|Site/Plot) + (Treatment + std_PC1|Species), popdata)
# lambdaslasimpledharma <- simulateResiduals(lambdaslasimple2)
# plot(lambdaslasimpledharma)
# #Okay, KS test significant
# summary(lambdaslasimple2)
# 
# #Simple model just treatment
# lambdaslasimple4 <- lmer(log_lambda ~ Treatment + std_SLA + Treatment:std_SLA + (1|Site/Plot) + (Treatment|Species), popdata)
# lambdaslasimpledharma <- simulateResiduals(lambdaslasimple4)
# plot(lambdaslasimpledharma)
# #Poor, KS test significant
# summary(lambdaslasimple4)
# 
# #Simpler model with PC1, PC2 and neighbours
# lambdasla <- lmer(log_lambda_p1 ~ std_PC1 + std_PC2 + neighbours01 + std_log_SLA + 
#                     std_PC1:std_log_SLA + std_PC2:std_log_SLA + neighbours01:std_log_SLA + (1|Site/Plot) + (neighbours01 + std_PC1 + std_PC2|Species), popdata)
# lambdaslasimpledharma <- simulateResiduals(lambdasla)
# plot(lambdaslasimpledharma)
# #Okay, KS test significant
# summary(lambdasla)
# #No significant interactions
# 
# ##LDMC - converges!
# lambdaldmc <- lmer(log_lambda_p1 ~ neighbours01 + Treatment + std_PC1 + std_PC2 + std_LDMC + 
#                      neighbours01:std_LDMC + Treatment:std_LDMC + std_PC1:std_LDMC + std_PC2:std_LDMC +
#                      (1|Site/Plot) + (neighbours01 + Treatment + std_PC1 + std_PC2|Species), popdata)
# lambdaldmcdharma <- simulateResiduals(lambdaldmc)
# plot(lambdaldmcdharma)
# #not great residuals, KS test significant
# summary(lambdaldmc)
# #No significant interactions
# 
# #D13C - converges!
# lambdaD13C <- lmer(log_lambda_p1 ~ neighbours01 + Treatment + std_PC1 + std_PC2 + std_log_D13C + 
#                      neighbours01:std_log_D13C + Treatment:std_log_D13C + std_PC1:std_log_D13C + std_PC2:std_log_D13C +
#                      (1|Site/Plot) + (neighbours01 + Treatment + std_PC1 + std_PC2|Species), popdata)
# lambdaD13Cdharma <- simulateResiduals(lambdaD13C)
# plot(lambdaD13Cdharma)
# #not great residuals, KS test significant
# summary(lambdaD13C)
# #PC1 and TreatmentWet significant
# r.squaredLR(lambdaD13C)
#### Plotting significant trait interactions ####
#using John's Catford script, help in lab retreat and my amendments
######## Fecundity SLA PC1 PC2
# Loading in the model from saved output
load("Output/Models/seedslamodel.RData")
#Matching SLA values to Species in the model
traits_for_model<-speciessla[match(rownames(coef(seedsla)$Species), speciessla$Species),c(1:2)]
#John: this bit grabs the random slopes for each species as well as their associated SEs
#NOT SURE THESE ARE THE RIGHT VALUES?** I changed this to [8,8] to be SLA column and row but
#not sure that the random effects part is correct. What does [2,2,] pull out??
#I think [2,2,] pulls out the random slope values for PC1 which should be correct as [2,2,] since it is 
# the first random slope fitted in the model and [1,1,] is the random intercept?
#Not sure why there are six columns and rows though, I would expect 5
varfix <- vcov(seedsla)[8,8]
re <- ranef(seedsla,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
#Calculating the slope values for the response of fecundity to PC1 as modulated by SLA
# slopes = the slope responses of PC2 for each species + the interaction of PC2:SLA*each species' SLA value
#Unique value for PC2 slopes per species (what's left over), one value for PC2:SLA for all species, unique values of SLA per species
#This changes the position of each species: PC2:SLA - which is one it's one value
# This is all because the fixed effects have explained some of the response already
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(seedsla)$Species[,6] + coef(seedsla)$Species[,9]*traits_for_model$std_log_SLA,
                                         slope_PC2 = coef(seedsla)$Species[,7] + coef(seedsla)$Species[,11]*traits_for_model$std_log_SLA)

### Plotting response of fecundity to PC1 as a function of SLA
#c(6,9) is pulling out std_PC1 and std_PC1*std_log_SLA
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(6,9)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(6,9), c(6,9)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ std_log_SLA, ylab= "Slope of fecundity - PC1 relationship", xlab="log(SLA) (standardised)", ylim=c(-1.1,.3), pch = 19, col="grey70"))
with(slopedata, arrows(std_log_SLA, slope_PC1+slope.SEs, std_log_SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
abline(h=0, lty=3, lwd=2)

#Calculating R squared
pc1sla <- lm(slope_PC1 ~ std_log_SLA, slopedata)
summary(pc1sla)

### Plotting response of fecundity to PC2 as a function of SLA
#c(7,11) is pulling out std_PC2 and std_PC2*std_log_SLA

data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(7,11)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(7,11), c(7,11)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC2 ~ std_log_SLA, ylab= "Slope of fecundity - PC2 relationship", xlab="log(SLA) (standardised)", ylim = c(-0.5, 0.4), pch = 19, col="grey70"))
with(slopedata, arrows(std_log_SLA, slope_PC2+slope.SEs, std_log_SLA, slope_PC2-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
abline(h=0, lty=3, lwd=2)

#Calculating R squared
pc2sla <- lm(slope_PC2 ~ std_log_SLA, slopedata)
summary(pc2sla)

#As a check, run a model without interaction and these should be the slope values we want ultimately
# after adding the things, this:
## slopes = the slope responses of PC2 for each species + the interaction of PC2:SLA*each species' SLA value
#Very similar - not exactly the same!

## NEED TO UPDATE MY MODIFIERS OF SLOPES FROM MY BIG MODEL 
# everything that is interacting with SLA! 
#Check this?*

######## Plotting lambda Treatment PC1
summary(lambdaD13C)
traits_for_model<-speciesD13C[match(rownames(coef(lambdaD13C)$Species), speciesD13C$Species),c(1:2)]
varfix <- vcov(lambdaD13C)[2,2]
re <- ranef(lambdaD13C,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(lambdaD13C)$Species[,5] + coef(lambdaD13C)$Species[,11]*traits_for_model$std_log_D13C)

data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdata$std_log_D13C))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(lambdaD13C)[c(5,11)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(lambdaD13C)[c(5,11), c(5,11)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ std_log_D13C, ylab= "Slope of lambda - PC1 relationship", ylim = c(-0.65, 0.45), xlab="log(D13C) (standardised)", pch = 19, col="grey70"))
with(slopedata, arrows(std_log_D13C, slope_PC1+slope.SEs, std_log_D13C, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
plot.CI.func(x.for.plot=seq.func(popdata$std_log_D13C), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
abline(h=0, lty=3, lwd=2)

#Calculating R squared
pc1d13c <- lm(slope_PC1 ~ std_log_D13C, slopedata)
summary(pc1d13c)

### Plotting significant TreatmentWet:std_D13C lambda response
#Just wet
traits_for_model<-speciesD13C[match(rownames(coef(lambdaD13C)$Species), speciesD13C$Species),c(1:2)]
varfix <- vcov(lambdaD13C)[2,2]
re <- ranef(lambdaD13C,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1_dry = coef(lambdaD13C)$Species[,3] + coef(lambdaD13C)$Species[,9]*traits_for_model$std_log_D13C,
                                         slope_PC1_wet = coef(lambdaD13C)$Species[,4] + coef(lambdaD13C)$Species[,10]*traits_for_model$std_log_D13C)
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdata$std_log_D13C))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(lambdaD13C)[c(4,10)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(lambdaD13C)[c(4,10), c(4,10)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1_wet ~ std_log_D13C, ylab= "Slope of lambda - TreatmentWet relationship", ylim = c(-1.0, 1.4), xlab="log(D13C) (standardised)", pch = 19, col="grey70"))
with(slopedata, arrows(std_log_D13C, slope_PC1_wet+slope.SEs, std_log_D13C, slope_PC1_wet-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
plot.CI.func(x.for.plot=seq.func(popdata$std_log_D13C), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
abline(h=0, lty=3, lwd=2)
#Calculating R squared
lambawetD13C <- lm(slope_PC1_wet ~ std_log_D13C, slopedata)
summary(lambawetD13C)

# y axis slope of response
# x axis D13C
# two lines - wet and dry (not sure how to add ambient)
summary(lambdaD13C)
traits_for_model<-speciesD13C[match(rownames(coef(lambdaD13C)$Species), speciesD13C$Species),c(1:2)]
varfix <- vcov(lambdaD13C)[2,2]
re <- ranef(lambdaD13C,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1_dry = coef(lambdaD13C)$Species[,3] + coef(lambdaD13C)$Species[,9]*traits_for_model$std_log_D13C,
                                         slope_PC1_wet = coef(lambdaD13C)$Species[,4] + coef(lambdaD13C)$Species[,10]*traits_for_model$std_log_D13C)

data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdata$std_log_D13C))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(lambdaD13C)[c(4,10)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(lambdaD13C)[c(4,10), c(4,10)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1_wet ~ std_log_D13C, ylab= "Slope of lambda - Treatment relationship", ylim = c(-1.0, 1.4), xlab="log(D13C) (standardised)", pch = 19, col="lightblue"))
with(slopedata, arrows(std_log_D13C, slope_PC1_wet+slope.SEs, std_log_D13C, slope_PC1_wet-slope.SEs, length = 0, angle = 30, code = 2, col = "lightblue", lwd = 1))
plot.CI.func(x.for.plot=seq.func(popdata$std_log_D13C), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='lightblue', env.trans=40, line.colour='blue', line.weight=3, line.type=1)
text(x = 2.0, y = 1.2, "*", cex = 2, col = "blue")
abline(h=0, lty=3, lwd=2)
with(slopedata, points(slope_PC1_dry ~ std_log_D13C, pch = 19, col="pink"))
with(slopedata, arrows(std_log_D13C, slope_PC1_dry+slope.SEs, std_log_D13C, slope_PC1_dry-slope.SEs, length = 0, angle = 30, code = 2, col = "pink", lwd = 1))
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdata$std_log_D13C))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(lambdaD13C)[c(3,9)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(lambdaD13C)[c(3,9), c(3,9)], as.matrix(data.for.sla.slope.regression))))
plot.CI.func(x.for.plot=seq.func(popdata$std_log_D13C), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='pink', env.trans=40, line.colour='red', line.weight=3, line.type=1)

#### Panel four traits ##########
###Making a panel with four figures representing significant trait relationships
#A - fecundity-PC1 modulated by SLA; C - fecundity-PC2 modulated by SLA
# B - lambda-PC1 modulated by D13C; D - lambda-TreatmentWet modulated by D13C

dev.off()
pdf("Output/Figures/panel_traits.pdf", width=21, height=21)
par(mfrow=c(2,2), mar =c(8,8,1,1), oma = c(1, 1, 1, 1))
#A
traits_for_model<-speciessla[match(rownames(coef(seedsla)$Species), speciessla$Species),c(1:2)]
varfix <- vcov(seedsla)[2,2]
re <- ranef(seedsla,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(seedsla)$Species[,6] + coef(seedsla)$Species[,9]*traits_for_model$std_log_SLA,
                                         slope_PC2 = coef(seedsla)$Species[,7] + coef(seedsla)$Species[,11]*traits_for_model$std_log_SLA)
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(6,9)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(6,9), c(6,9)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ std_log_SLA, ylab= "Slope of fecundity - PC1 relationship", xlab=NA, ylim=c(-1.1,.3), pch = 19, col="grey70", cex.lab = 3, cex.axis = 2, cex = 3))
with(slopedata, arrows(std_log_SLA, slope_PC1+slope.SEs, std_log_SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
abline(h=0, lty=3, lwd=3, cex = 3)
#B
traits_for_model<-speciesD13C[match(rownames(coef(lambdaD13C)$Species), speciesD13C$Species),c(1:2)]
varfix <- vcov(lambdaD13C)[2,2]
re <- ranef(lambdaD13C,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(lambdaD13C)$Species[,5] + coef(lambdaD13C)$Species[,11]*traits_for_model$std_log_D13C)
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdata$std_log_D13C))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(lambdaD13C)[c(5,11)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(lambdaD13C)[c(5,11), c(5,11)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ std_log_D13C, ylab= "Slope of lambda - PC1 relationship", ylim = c(-0.65, 0.45), xlab=NA, pch = 19, col="grey70", cex.lab = 3, cex.axis = 2, cex = 3))
with(slopedata, arrows(std_log_D13C, slope_PC1+slope.SEs, std_log_D13C, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
plot.CI.func(x.for.plot=seq.func(popdata$std_log_D13C), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
abline(h=0, lty=3, lwd=3, cex = 3)
#C
traits_for_model<-speciessla[match(rownames(coef(seedsla)$Species), speciessla$Species),c(1:2)]
varfix <- vcov(seedsla)[2,2]
re <- ranef(seedsla,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(seedsla)$Species[,6] + coef(seedsla)$Species[,9]*traits_for_model$std_log_SLA,
                                         slope_PC2 = coef(seedsla)$Species[,7] + coef(seedsla)$Species[,11]*traits_for_model$std_log_SLA)
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(7,11)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(7,11), c(7,11)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC2 ~ std_log_SLA, ylab= "Slope of fecundity - PC2 relationship", xlab="log(SLA) (standardised)", ylim = c(-0.5, 0.4), cex.lab = 3, cex.axis = 2, pch = 19, col="grey70", cex = 3))
with(slopedata, arrows(std_log_SLA, slope_PC2+slope.SEs, std_log_SLA, slope_PC2-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
abline(h=0, lty=3, lwd=3, cex = 3)
#D
traits_for_model<-speciesD13C[match(rownames(coef(lambdaD13C)$Species), speciesD13C$Species),c(1:2)]
varfix <- vcov(lambdaD13C)[2,2]
re <- ranef(lambdaD13C,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
slopedata <- traits_for_model %>% mutate(slope_PC1_dry = coef(lambdaD13C)$Species[,3] + coef(lambdaD13C)$Species[,9]*traits_for_model$std_log_D13C,
                                         slope_PC1_wet = coef(lambdaD13C)$Species[,4] + coef(lambdaD13C)$Species[,10]*traits_for_model$std_log_D13C)
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(popdata$std_log_D13C))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(lambdaD13C)[c(4,10)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(lambdaD13C)[c(4,10), c(4,10)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1_wet ~ std_log_D13C, ylab= "Slope of lambda - TreatmentWet relationship", ylim = c(-1.0, 1.4), xlab="log(D13C) (standardised)", cex.lab = 3, cex.axis = 2, pch = 19, col="grey70", cex = 3))
with(slopedata, arrows(std_log_D13C, slope_PC1_wet+slope.SEs, std_log_D13C, slope_PC1_wet-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1, cex = 3))
plot.CI.func(x.for.plot=seq.func(popdata$std_log_D13C), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=4, line.type=1)
abline(h=0, lty=3, lwd=3, cex = 3)

dev.off()

