### Answering my specific questions
# WA Perenjori Watering Experiment
#October 2021
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(MuMIn)
library(DHARMa)
library(glmmTMB)

#Functions file
source("R_functions/functions.R")
#Data imported from data preparation sheet
source("data_preparation.R")
#dataall has everything (germination, survival, seed production and neighbour info) combined

#### Question 1 ####
##What is the relative importance of abiotic and biotic factors for survival and fecundity?
#Additive model: response ~ total_abundance + water + PC1 + PC2 + PC3 + RE
### Survival to produce seeds ####
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survival <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + 
                  (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]))
  print(summary(survival))
}
###Putting models into a list for coef plotting
survivalmodels <- list()
for (i in 1:length(specieslist)){
  survivalmodels[[i]] <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + Dodder01 +
                                 (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]))
}
#Making coefficient plot for survival as a function of environment
plot_models(survivalmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
### Assigning model names by species to plot
for (i in 1:length(specieslist)){
  nam <- paste0("survmod", specieslist[i])
  assign(nam, glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + Dodder01 +
                      (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i])))
}

#Want to plot survival in response to neighbour abundance
#ARCA survival model, model created in loop above
arcasurvdharma <- simulateResiduals(survmodARCA)
plot(arcasurvdharma)
summary(survmodARCA)
x_to_plot<-seq.func(dataall$std_logp1_totalabund)
with(arcadata, plot(ProducedSeeds ~ std_logp1_totalabund))
arcapreddata <- with(arcadata, data.frame(1, x_to_plot, 0, 0, 0, 0, 0, 0))
arcapred <- glmm.predict(mod = survmodARCA, newdat = arcapreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

### Relative importance abiotic and biotic - Fecundity #####
fecunditymodels <- list()
for (i in 1:length(specieslist)){
  print(specieslist[i])
  fecundity <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + Dodder01 +
                      (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(fecundity))
}

#Assigning model names by species
for (i in 1:length(specieslist)){
  nam <- paste0("seedmod", specieslist[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + 
                      (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i])))
}
#Putting models in a list for tab_model and coef plots
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

#### Example plotting model output - this works!
#Plotting number of viable seeds against PC1, holding total abundance at mean (0, standardised)
arcaseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 +
                          (1|Site/Plot), family = nbinom2, data = seedarca)
arcaseedmod1test <- simulateResiduals(arcaseedmod1)
plot(arcaseedmod1test)
summary(arcaseedmod1)
x_to_plot<-seq.func(seedarca$std_logp1_totalabund)
with(seedarca, plot(No_viable_seeds_grouped ~ seedarca$std_logp1_totalabund))
arcapreddata <- with(seedarca, data.frame(1, x_to_plot, 0, 0, 0))
arcapred <- glmm.predict(mod = arcaseedmod1, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

## Trying above but with the full model
arcaseedmod2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + Dodder01 +
                          (1|Site/Plot), family = nbinom2, data = seedarca)
arcaseedmod2resid <- simulateResiduals(arcaseedmod2)
plot(arcaseedmod2resid)
#This is telling me it's unhappy with the model!?
arcaseedmod3 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 +
                          (1|Site/Plot), family = nbinom2, data = seedarca)
arcaseedmod3resid <- simulateResiduals(arcaseedmod3)
plot(arcaseedmod3resid)
# Much happier with this

# Plotting all species' responses to PC1 which is coming up as sig. for 4 species 
#with TA + Treatment + PC1/2/3 + Dodder
# But ARCA isn't happy with that (bad residuals)
arcaseedmod2resid <- simulateResiduals(arcaseedmod2)
plot(arcaseedmod2resid)
seedmodHYGLresid <- simulateResiduals(seedmodHYGL)
plot(seedmodHYGLresid)
#HYGL fine
seedmodLAROresid <- simulateResiduals(seedmodLARO)
plot(seedmodLAROresid)
#LARO fine
seedmodPOLEresid <- simulateResiduals(seedmodPOLE)
plot(seedmodPOLEresid)
#POLE fine
seedmodPEAIresid <- simulateResiduals(seedmodPEAI)
plot(seedmodPEAIresid)
#PEAI unhappy
seedmodPLDEresid <- simulateResiduals(seedmodPLDE)
plot(seedmodPLDEresid)
#PLDE unhappy
seedmodTRCYresid <- simulateResiduals(seedmodTRCY)
plot(seedmodTRCYresid)
#TRCY fine
seedmodTRORresid <- simulateResiduals(seedmodTROR)
plot(seedmodTRORresid)
#TROR unhappy
seedmodVEROresid <- simulateResiduals(seedmodVERO)
plot(seedmodVEROresid)
#VERO unhappy

#What if I run model without Dodder01?
# Viableseeds ~ TA + Treatment + PC1/2/3
#PC1 signif 3 times.
seedmodARCAresid <- simulateResiduals(seedmodARCA)
plot(seedmodARCAresid)
#ARCA fine
seedmodHYGLresid <- simulateResiduals(seedmodHYGL)
plot(seedmodHYGLresid)
#HYGL fine
seedmodLAROresid <- simulateResiduals(seedmodLARO)
plot(seedmodLAROresid)
#LARO fine
seedmodPOLEresid <- simulateResiduals(seedmodPOLE)
plot(seedmodPOLEresid)
#POLE very unhappy
seedmodPEAIresid <- simulateResiduals(seedmodPEAI)
plot(seedmodPEAIresid)
#PEAI fine
seedmodPLDEresid <- simulateResiduals(seedmodPLDE)
plot(seedmodPLDEresid)
#PLDE unhappy
seedmodTRCYresid <- simulateResiduals(seedmodTRCY)
plot(seedmodTRCYresid)
#TRCY fine
seedmodTRORresid <- simulateResiduals(seedmodTROR)
plot(seedmodTRORresid)
#TROR fine
seedmodVEROresid <- simulateResiduals(seedmodVERO)
plot(seedmodVEROresid)
#VERO unhappy

# Running a model with just viable seeds ~ TA + Treatment + PC1 + PC2
laro1 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 +
                          (1|Site/Plot), family = nbinom2, data = seedlaro)
laro1resids <- simulateResiduals(laro1)
plot(laro1resids)
#Significant with PC1 are: LARO, PEAI and TRCY
#ARCA fine
#HYGL fine
#LARO fine
#PEAI fine
#PLDE unhappy
#POLE very unhappy
#TRCY fine
#TROR fine
#VERO unhappy

par(mfrow=c(1,1))
summary(mod)
x_to_plot<-seq.func(seedlaro$std_PC1)
with(seedlaro, plot(No_viable_seeds_grouped ~ seedlaro$std_PC1))
arcapreddata <- with(seedlaro, data.frame(1, 0, 0, 0, x_to_plot, 0))
arcapred <- glmm.predict(mod = mod, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Know that I can plot these in a panel, but not sure which models
# to plot since some species have poor residuals


##### Question 2 ####
## Do species-level leaf traits explain responses to abiotic environment?
##Starting with simple example. Does SLA explain the relationship between
# survival and PC1?
### Simple survival model, example
survpc1sla <- glmer(ProducedSeeds ~ std_PC1 + std_log_SLA + std_PC1:std_log_SLA + (1|Site/Plot) + 
                      (std_PC1|Species), family = binomial, dataall)
summary(survpc1sla)
survpc1ldmc <- glmer(ProducedSeeds ~ std_PC1 + std_LDMC + std_PC1:std_LDMC + (1|Site/Plot) + 
                      (std_PC1|Species), family = binomial, dataall)
summary(survpc1ldmc)
survpc1wue <- glmer(ProducedSeeds ~ std_PC1 + std_log_D13C + std_PC1:std_log_D13C + (1|Site/Plot) + 
                       (std_PC1|Species), family = binomial, dataall)
summary(survpc1wue)
r.squaredGLMM(survpc1wue)
#WUE by itself seems to explain survival.
with(dataall, plot(ProducedSeeds ~ std_log_D13C))
x_to_plot<-seq.func(dataall$std_log_D13C)
wuepreddata <- with(dataall, data.frame(1, 0, x_to_plot, 0*x_to_plot))
wuepred <- glmm.predict(mod = survpc1wue, newdat = wuepreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = wuepred$y, upper = wuepred$upper, lower = wuepred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Full survival model
survtraitsmod <- glmer(ProducedSeeds ~ Treatment + std_PC1 + std_log_SLA + Treatment:std_log_SLA + std_PC1:std_log_SLA + 
                                  std_LDMC + Treatment:std_LDMC + std_PC1:std_LDMC + 
                                  std_log_D13C + Treatment:std_log_D13C + std_PC1:std_log_D13C + 
                                  (1|Site/Plot) + (Treatment + std_PC1|Species), 
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

#### Calculating species richness ####
# Column for the number of neighbouring species that each focal has

#sp_richness = count of number of rows? but dodder. tally of neighbour sp?
#How to consider dodder?
#vegan - specnumber
library(vegan)

#Not working yet
#test <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep)
#  diversity(seedabundance$Neighbour_sp, index = "shannon")
#specnumber(seedabundance$Neighbour_sp)
#Let's calculate Shannon's Index rather than just species number
#Neighbour count is the number of neighbours in that row.
#Neighbour_sp is the species.
#Get formula for Shannon's index and manually calculate it?
  #will have to multiply neighbour by neighbour_sp




