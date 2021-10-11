### Answering my specific questions
# WA Perenjori Watering Experiment
#October 2021
library(tidyverse)  
library(lme4)
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(emmeans)

#Data imported from data preparation sheet
source("data_preparation.R")
#dataall has everything (germination, survival, seed production and neighbour info) combined
#datanotscaled is same dataset but without standardised predictors
#Functions file
source("R_functions/functions.R")

#### Question 1 ####
##What is the relative importance of abiotic and biotic factors for
## germination, survival and fecundity?
#Additive model: response ~ total_abundance + water + PC1 + PC2 + PC3 + RE
### Relative importance abiotic and biotic - Survival to produce seeds ###
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survival <- glmer(ProducedSeeds ~ log(Total_abundance+1) + Treatment + PC1 + PC2 + PC3 + 
                  (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]))
  print(summary(survival))
}
#Assigning model names by species
survivalmodels <- list()
for (i in 1:length(specieslist)){
  survivalmodels[[i]] <- glmer(ProducedSeeds ~ log(Total_abundance+1) + Treatment + PC1 + PC2 + PC3 + 
                                 (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]))
}
plot_models(survivalmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

### Relative importance abiotic and biotic - Fecundity ##
#Loops aren't working because of convergence issues with ARCA
fecunditymodels <- list()
for (i in 1:length(specieslist)){
  print(specieslist[i])
  fecundity <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + Treatment + PC1 + PC2 + PC3 +
                      (1|Site/Plot), data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(fecundity))
}
# #Assigning model names by species
# for (i in 1:length(specieslist)){
#   nam <- paste0("seedmod", specieslist[i])
#   assign(nam, glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + Treatment + PC1 + PC2 + PC3 +
#                       (1|Site/Plot), data = filter(seedmodeldata, Species == specieslist[i])))
# }
# arcaseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedarca)
# summary(arcaseedmod)  
# hyglseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedhygl)
# summary(hyglseedmod)
# laroseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedlaro)
# summary(laroseedmod)
# peaiseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedpeai)
# summary(peaiseedmod)
# pldeseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedplde)
# summary(pldeseedmod)
# poleseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedpole)
# summary(poleseedmod)
# trcyseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedtrcy)
# summary(trcyseedmod)
#####Why isn't this converging?!
trorseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedtror)
# summary(trorseedmod)
# veroseedmod <- glmer.nb(No_viable_seeds_grouped ~ log(Total_abundance+1) + 
#                           Treatment + PC1 + PC2 + PC3 + (1|Site/Plot), seedvero)
# summary(veroseedmod)
seedmodels <- list()
seedmodels[[1]] <- seedmodARCA
seedmodels[[2]] <- seedmodHYGL
seedmodels[[3]] <- seedmodLARO
seedmodels[[4]] <- seedmodPEAI
seedmodels[[5]] <- seedmodPLDE
seedmodels[[6]] <- seedmodPOLE
seedmodels[[7]] <- seedmodTRCY
seedmodels[[8]] <- seedmodTROR
seedmodels[[9]] <- seedmodVERO
tab_model(seedmodels, transform = NULL)
plot_models(seedmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

##### Question 2 ####
## Do species-level leaf traits explain responses to abiotic environment?
#Distributions of traits and PC1, have to use datanotscaled to log
#SLA better logged?
#LDMC still strange logged
#WUE/D13C better logged but still left-skewed
#PC1 - right-skewed, logging doesn't help
with(datanotscaled, pairs(~SLA + LDMC + mean_D13C + PC1, diag.panel = panel.hist))
with(datanotscaled, pairs(~log(SLA) + LDMC + log(mean_D13C) + PC1, diag.panel = panel.hist))

with(dataall, pairs(~SLA + LDMC + mean_D13C + PC1, diag.panel = panel.hist))
with(dataall, pairs(~log(SLA+2) + LDMC + log(mean_D13C+2) + PC1, diag.panel = panel.hist))

#Is this a legitimate solution?! Trying to have a fixed log, log+2 of scaled data

##Starting with simple example. Does SLA explain the relationship between
# survival and PC1?
### Simple survival model, example
survpc1sla <- glmer(ProducedSeeds ~ PC1 + log(SLA+2) + PC1:log(SLA+2) + (1|Site/Plot) + 
                      (PC1|Species), family = binomial, dataall)
summary(survpc1sla)
survpc1ldmc <- glmer(ProducedSeeds ~ PC1 + LDMC + PC1:LDMC + (1|Site/Plot) + 
                      (PC1|Species), family = binomial, dataall)
summary(survpc1ldmc)
survpc1wue <- glmer(ProducedSeeds ~ PC1 + mean_D13C + PC1:mean_D13C + (1|Site/Plot) + 
                       (PC1|Species), family = binomial, dataall)
summary(survpc1wue)

#Full survival model
survtraitsmod <- glmer(ProducedSeeds ~ Treatment + PC1 + SLA + Treatment:SLA + PC1:SLA + 
                                  LDMC + Treatment:LDMC + PC1:LDMC + 
                                  mean_D13C + Treatment:mean_D13C + PC1:mean_D13C + 
                                  (1|Site/Plot) + (Treatment + PC1|Species), 
                                  family = binomial, dataall)
#Failed to converge
summary(survtraitsmod)

### Simple fecundity model, example
seedpc1sla <- glmer.nb(No_viable_seeds_grouped ~ Treatment + PC1 + SLA +
                         Treatment:SLA + PC1:SLA + (1|Site/Plot) + 
                         (Treatment + PC1 |Species), dataall)
summary(seedpc1sla)

#Full fecundity model
seedtraitsmod <- glmer.nb(No_viable_seeds_grouped ~ Treatment + PC1 + 
                            log(SLA+2) + Treatment:log(SLA+2) + PC1:log(SLA+2) + 
                            LDMC + Treatment:LDMC + PC1:LDMC + 
                            log(mean_D13C+2) + Treatment:log(mean_D13C+2) + PC1:log(mean_D13C+2) + 
                            (1|Site/Plot) + (Treatment + PC1|Species), dataall)
summary(seedtraitsmod)

#Need to add biotic elements into this model? Total_abundance





