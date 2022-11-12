### Answering my specific questions
# WA Perenjori 2020 Watering Experiment
#Updated July 2022

#### Loading packages and data ####
library(ggplot2)
library(lmerTest)
library(sjPlot)
library(MuMIn)
library(DHARMa)
library(glmmTMB)
library(gridExtra)
library(ggExtra)
library(cowplot)
library(car)
library(kableExtra)

#Functions file
source("R_functions/functions.R")
#ggplot here for ease of plotting. Use theme_classic()+ my_theme
my_theme <- theme(axis.title.x = element_text(size = 14, face = 'bold'),
                  axis.title.y = element_text(size = 14, face = 'bold'),
                  axis.text = element_text(size = 14),
                  strip.text.x = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.title = element_blank())
#Data imported from data preparation sheet
source("data_preparation.R")
#vitaldata has all datasets combined 
# germination, survival, seed production, neighbour info, traits, abiotic factors
# one row per subplot with seeds sown - NAs are very important since, e.g. survival info is only on germinated subplots
#1614 subplots with germination data (1171 germinated)
#1139* subplots with survival data (1171 have surv_to_produce_seeds data, 1139 from survtomerge, mortalitydataraw 1164)
#653 subplots with seed production data
#1100 subplots with neighbourhood surveys done (534 have neighbours)

## Note that germination analysis is separate, in germination_analysis.R
## Note that neighbour abundance and PC2 are correlated, so will not model them together

#### Is neighbour abundance correlated with abiotic environmental factors or diversity? ####
##PC1
pc1plot1 <- ggplot(vitaldata, aes(x = std_PC1, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()
pc1plot2 <- ggplot(datanonly, aes(x = std_PC1, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("")+
  theme_classic()
#Testing marginal histogram plot or marginal density plot
ggMarginal(pc1plot2)
ggMarginal(pc1plot2, type = "histogram")
##PC2
pc2plot1 <- ggplot(vitaldata, aes(x = std_PC2, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()
pc2plot2 <- ggplot(datanonly, aes(x = std_PC2, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("")+
  theme_classic()

##Water availability
waterplot1 <- ggplot(vitaldata, aes(x = Treatment, y = log(Total_abundance+1)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()
waterplot2 <- ggplot(datanonly, aes(x = Treatment, y = log(Total_abundance+1)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  ylab("")+
  theme_classic()
pdf("Output/Figures/n_abundance_corr.pdf")
  #grid.arrange(pc1plot1, pc1plot2, pc2plot1, pc2plot2, waterplot1, waterplot2, ncol = 2)
#hjust adjusts position of label in cowplot
#Interestingly different position for E and F compared to rest
plot_grid(pc1plot1, pc1plot2, pc2plot1, pc2plot2, waterplot1, waterplot2, 
          ncol = 2, labels = c("A", "B", "C", "D", "E", "F"), hjust = -3.5)+
          theme(plot.margin = unit(c(5,1,1,1), "points"))
dev.off()

##Plotting only this without thinned neighbourhoods, and not fitting linear relationship
pc1plot2b <- ggplot(datanonly, aes(x = std_PC1, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth()+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()+
  my_theme
pc2plot2b <- ggplot(datanonly, aes(x = std_PC2, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth()+
  ylab("")+
  theme_classic()+
  my_theme
pdf("Output/Figures/n_abund_corr_smooth.pdf")
plot_grid(pc1plot2b, pc2plot2b, labels = c("A", "B"))+
  theme(plot.margin = unit(c(5,1,1,1), "points"))
dev.off()

#Models
modelabundpc1 <- lm(log(Total_abundance+1) ~ std_PC1, datanonly)
summary(modelabundpc1)
modelabundpc2 <- lm(log(Total_abundance+1) ~ std_PC2, datanonly)
summary(modelabundpc2)
modelabundpc2b <- lm(log(Total_abundance+1) ~ std_PC2, vitaldata)
summary(modelabundpc2b)
modelabundpc2c <- lm(log(Total_abundance+1) ~ std_PC2 + std_PC1, vitaldata)
summary(modelabundpc2c)
vif(modelabundpc2c)
modelabundpc2d <- lm(No_viable_seeds_grouped ~ std_logp1_totalabund + std_PC1 + std_PC2, vitaldata)
summary(modelabundpc2d)
vif(modelabundpc2d)
modelabundwater <- aov(log(Total_abundance+1) ~ Treatment, datanonly)
summary(modelabundwater)
TukeyHSD(modelabundwater)

### More plots
##Plotting PC1 vs abund for watering treatments
waterpc1plot <- ggplot(datanonly, aes(x = std_PC1, y = logp1_totalabund))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  ylab("log(neighbour abundance + 1)")+
  theme_classic()+
  my_theme+
  facet_wrap(~Treatment)
##Plotting PC2 vs abund for watering treatments
waterpc2plot <- ggplot(datanonly, aes(x = std_PC2, y = logp1_totalabund))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  ylab("log(neighbour abundance + 1)")+
  theme_classic()+
  my_theme+
  facet_wrap(~Treatment)
pdf("Output/Figures/PCs_abund_by_watering.pdf")
plot_grid(waterpc1plot, waterpc2plot, ncol = 1, labels = c("A", "B")) + 
  theme(plot.margin = unit(c(1,50,1,50), "points"))
dev.off()

### How is neighbour abundance distributed across species across abiotic env? ####
sppc1abund <- ggplot(vitaldata, aes(x = std_PC1, y = logp1_totalabund))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()+
  facet_wrap(~Species)
sppc2abund <- ggplot(vitaldata, aes(x = std_PC2, y = logp1_totalabund))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()+
  facet_wrap(~Species)
# Most species don't have data for high neighbour abundance at low values of PC2
# POLE doesn't have much coverage for low values of PC2 or PC1 (lesser extent)
# All in all it's great coverage!
#Try cowplot for plotting
#The order of edges for plot.margin is unit(c(top, right, bottom, left), units)
pdf("Output/Figures/nbh_coverage.pdf")
plot_grid(sppc1abund, sppc2abund, ncol = 1, labels = c("A", "B")) + 
  theme(plot.margin = unit(c(1,50,1,50), "points"))
dev.off()

#### Does the diversity of neighbours influence vital rates? ####
# Can only look at this in subplots that have neighbours, otherwise data are zero conflated with respect to SDI and neighbour abundance
#datanonly is vitaldata filtered to Total_abundance > 0 

### Is neighbour abundance correlated with diversity?
# Simpson's diversity index, Shannon Diversity Index and species richness
ggplot(datanonly, aes(x = logp1_totalabund, y = log(sp_richness)))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

modelabundrichness <- lm(logp1_totalabund ~ log(sp_richness), datanonly)
summary(modelabundrichness)

ggplot(datanonly, aes(x = logp1_totalabund, y = SDI))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

modelabundrichness <- lm(logp1_totalabund ~ SDI, datanonly)
summary(modelabundrichness)

ggplot(datanonly, aes(x = logp1_totalabund, y = shannon))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

modelabundrichness <- lm(logp1_totalabund ~ shannon, datanonly)
summary(modelabundrichness)

ggplot(datanonly, aes(x = SDI, y = log(sp_richness)))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

modelrichnessSDI <- lm(SDI ~ sp_richness, datanonly)
summary(modelrichnessSDI)

ggplot(datanonly, aes(x = SDI, y = shannon))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

modelrichnessSDI <- lm(SDI ~ shannon, datanonly)
summary(modelrichnessSDI)

### Testing variance inflation factors with species richness 
arcatest <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + sp_richness + (1|Site/Plot),
                  family = binomial, arcadata)
summary(arcatest)
vif(arcatest)
#sp richness > 3!

arcatest2 <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                   family = binomial, arcadata)
summary(arcatest2)
vif(arcatest2)

arcatest4 <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + shannon + (1|Site/Plot),
                   family = binomial, arcadata)
summary(arcatest4)
vif(arcatest4)
#shannon and SDI 7

#only using shannon! as a measure of diversity

#Create a model per species
# Survival
#without POLE
for (i in 1:length(specieslist.nop)){
  print(specieslist.nop[i])
  survival <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + shannon + (1|Site/Plot),
                    family = binomial, data = filter(datanonly, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(survival))
}

#TRCY and TROR significant positive response to SDI

### Quantifying variance in vital rates explained by abiotic only vs biotic only ####
#### r squared table ####
### Survival ###
#run models from 'AIC Model Selection Abiotic vs Biotic' section below
#abioticmodsurv
#bioticmodsurv
#abioticmodseed
#bioticmodseed

#Extracting values for theoretical marginal R squared
summary(abioticmodsurvARCA)

arcaabiotic <- r.squaredGLMM(abioticmodsurvARCA)[1,1]
hyglabiotic <- r.squaredGLMM(abioticmodsurvHYGL)[1,1]
laroabiotic <- r.squaredGLMM(abioticmodsurvLARO)[1,1]
peaiabiotic <- r.squaredGLMM(abioticmodsurvPEAI)[1,1]
pldeabiotic <- r.squaredGLMM(abioticmodsurvPLDE)[1,1]
trcyabiotic <- r.squaredGLMM(abioticmodsurvTRCY)[1,1]
trorabiotic <- r.squaredGLMM(abioticmodsurvTROR)[1,1]
veroabiotic <- r.squaredGLMM(abioticmodsurvVERO)[1,1]

### Biotic only - additive
# Need optimiser to converge
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survivalbiotic <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Dodder01 + shannon + (1|Site/Plot), 
                          family = binomial, data = filter(vitaldata, Species == specieslist[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(survivalbiotic))
}

for (i in 1:length(specieslist)){
  nam <- paste0("bioticmod", specieslist[i])
  assign(nam, glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Dodder01 + shannon + (1|Site/Plot), 
                    family = binomial, data = filter(vitaldata, Species == specieslist[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}

summary(bioticmodARCA)
arcabiotic <- r.squaredGLMM(bioticmodARCA)[1,1]
summary(bioticmodHYGL)
hyglbiotic <- r.squaredGLMM(bioticmodHYGL)[1,1]
summary(bioticmodLARO)
larobiotic <- r.squaredGLMM(bioticmodLARO)[1,1]
summary(bioticmodPEAI)
peaibiotic <- r.squaredGLMM(bioticmodPEAI)[1,1]
summary(bioticmodPLDE)
pldebiotic <- r.squaredGLMM(bioticmodPLDE)[1,1]
summary(bioticmodPOLE)
polebiotic <- r.squaredGLMM(bioticmodPOLE)[1,1]
summary(bioticmodTRCY)
trcybiotic <- r.squaredGLMM(bioticmodTRCY)[1,1]
summary(bioticmodTROR)
trorbiotic <- r.squaredGLMM(bioticmodTROR)[1,1]
summary(bioticmodVERO)
verobiotic <- r.squaredGLMM(bioticmodVERO)[1,1]

#Create table to put these values in
#create matrix with 2 columns filled with random value, 1
rsquaredtable <- matrix(rep(1, times=2), ncol=2, nrow = 9)
#define column names and row names of matrix
colnames(rsquaredtable) <- c('R squared biotic model', 'R squared abiotic model')
rownames(rsquaredtable) <- c('ARCA', 'HYGL', 'LARO', 'PEAI', 'POLE', 'PLDE', 'TRCY', 'TROR', 'VERO')
#convert matrix to table 
rsquaredtable <- as.data.frame(rsquaredtable)
#view table 
rsquaredtable
##Replace values with appropriate ones
#Biotic model
rsquaredtable[1,1] <- arcabiotic
rsquaredtable[2,1] <- hyglbiotic
rsquaredtable[3,1] <- larobiotic
rsquaredtable[4,1] <- peaibiotic
rsquaredtable[5,1] <- pldebiotic
rsquaredtable[6,1] <- polebiotic
rsquaredtable[7,1] <- trcybiotic
rsquaredtable[8,1] <- trorbiotic
rsquaredtable[9,1] <- verobiotic
#Abiotic model
rsquaredtable[1,2] <- arcaabiotic
rsquaredtable[2,2] <- hyglabiotic
rsquaredtable[3,2] <- laroabiotic
rsquaredtable[4,2] <- peaiabiotic
rsquaredtable[5,2] <- pldeabiotic
rsquaredtable[6,2] <- poleabiotic
rsquaredtable[7,2] <- trcyabiotic
rsquaredtable[8,2] <- trorabiotic
rsquaredtable[9,2] <- veroabiotic

#Export table as a csv file
write.csv(rsquaredtable,"Output/Tables/surv_rsquaredtable.csv")

####
### Fecundity ###
### Abiotic only
for (i in 1:length(specieslist)){
  print(specieslist[i])
  abioticseed <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                         family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(abioticseed))
}

for (i in 1:length(specieslist)){
  nam <- paste0("abioticmodseed", specieslist[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i])))
}
#Extracting values for theoretical marginal R squared
summary(abioticmodseedARCA)
arcaabioticseed <- r.squaredGLMM(abioticmodseedARCA)[1,1]
summary(abioticmodseedHYGL)
hyglabioticseed <- r.squaredGLMM(abioticmodseedHYGL)[1,1]
summary(abioticmodseedLARO)
laroabioticseed <- r.squaredGLMM(abioticmodseedLARO)[1,1]
summary(abioticmodseedPEAI)
peaiabioticseed <- r.squaredGLMM(abioticmodseedPEAI)[1,1]
summary(abioticmodseedPLDE)
pldeabioticseed <- r.squaredGLMM(abioticmodseedPLDE)[1,1]
summary(abioticmodseedPOLE)
poleabioticseed <- r.squaredGLMM(abioticmodseedPOLE)[1,1]
summary(abioticmodseedTRCY)
trcyabioticseed <- r.squaredGLMM(abioticmodseedTRCY)[1,1]
summary(abioticmodseedTROR)
trorabioticseed <- r.squaredGLMM(abioticmodseedTROR)[1,1]
summary(abioticmodseedVERO)
veroabioticseed <- r.squaredGLMM(abioticmodseedVERO)[1,1]

### Biotic only
for (i in 1:length(specieslist)){
  print(specieslist[i])
  bioticseed <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Dodder01 + shannon + (1|Site/Plot), 
                        family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(bioticseed))
}

for (i in 1:length(specieslist)){
  nam <- paste0("bioticmodseed", specieslist[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Dodder01 + shannon + (1|Site/Plot), 
                      family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i])))
}

summary(bioticmodseedARCA)
arcabioticseed <- r.squaredGLMM(bioticmodseedARCA)[1,1]
summary(bioticmodseedHYGL)
hyglbioticseed <- r.squaredGLMM(bioticmodseedHYGL)[1,1]
summary(bioticmodseedLARO)
larobioticseed <- r.squaredGLMM(bioticmodseedLARO)[1,1]
summary(bioticmodseedPEAI)
peaibioticseed <- r.squaredGLMM(bioticmodseedPEAI)[1,1]
summary(bioticmodseedPLDE)
pldebioticseed <- r.squaredGLMM(bioticmodseedPLDE)[1,1]
summary(bioticmodseedPOLE)
polebioticseed <- r.squaredGLMM(bioticmodseedPOLE)[1,1]
summary(bioticmodseedTRCY)
trcybioticseed <- r.squaredGLMM(bioticmodseedTRCY)[1,1]
summary(bioticmodseedTROR)
trorbioticseed <- r.squaredGLMM(bioticmodseedTROR)[1,1]
summary(bioticmodseedVERO)
verobioticseed <- r.squaredGLMM(bioticmodseedVERO)[1,1]

#Create table to put these values in
#create matrix with 2 columns filled with random value, 1
rsquaredtableseed <- matrix(rep(1, times=2), ncol=2, nrow = 9)
#define column names and row names of matrix
colnames(rsquaredtableseed) <- c('R squared biotic model', 'R squared abiotic model')
rownames(rsquaredtableseed) <- c('ARCA', 'HYGL', 'LARO', 'PEAI', 'POLE', 'PLDE', 'TRCY', 'TROR', 'VERO')
#convert matrix to table 
rsquaredtableseed <- as.data.frame(rsquaredtableseed)
#view table 
rsquaredtableseed
##Replace values with appropriate ones
#Biotic model
rsquaredtableseed[1,1] <- arcabioticseed
rsquaredtableseed[2,1] <- hyglbioticseed
rsquaredtableseed[3,1] <- larobioticseed
rsquaredtableseed[4,1] <- peaibioticseed
rsquaredtableseed[5,1] <- pldebioticseed
rsquaredtableseed[6,1] <- polebioticseed
rsquaredtableseed[7,1] <- trcybioticseed
rsquaredtableseed[8,1] <- trorbioticseed
rsquaredtableseed[9,1] <- verobioticseed
#Abiotic model
rsquaredtableseed[1,2] <- arcaabioticseed
rsquaredtableseed[2,2] <- hyglabioticseed
rsquaredtableseed[3,2] <- laroabioticseed
rsquaredtableseed[4,2] <- peaiabioticseed
rsquaredtableseed[5,2] <- pldeabioticseed
rsquaredtableseed[6,2] <- poleabioticseed
rsquaredtableseed[7,2] <- trcyabioticseed
rsquaredtableseed[8,2] <- trorabioticseed
rsquaredtableseed[9,2] <- veroabioticseed

#Export table as a csv file
write.csv(rsquaredtableseed,"Output/Tables/seed_rsquaredtable.csv")

###### Population growth rate, lambda ###
### Abiotic only
for (i in 1:length(specieslist)){
  print(specieslist[i])
  abioticlambda <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                        data = filter(popdata, Species == specieslist[i]))
  print(summary(abioticlambda))
}

for (i in 1:length(specieslist)){
  nam <- paste0("abioticmodlambda", specieslist[i])
  assign(nam, glmmTMB(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                      data = filter(popdata, Species == specieslist[i])))
}
#Extracting values for theoretical marginal R squared
summary(abioticmodlambdaARCA)
arcaabioticlambda <- r.squaredGLMM(abioticmodlambdaARCA)[1,1]
summary(abioticmodlambdaHYGL)
hyglabioticlambda <- r.squaredGLMM(abioticmodlambdaHYGL)[1,1]
summary(abioticmodlambdaLARO)
laroabioticlambda <- r.squaredGLMM(abioticmodlambdaLARO)[1,1]
summary(abioticmodlambdaPEAI)
peaiabioticlambda <- r.squaredGLMM(abioticmodlambdaPEAI)[1,1]
summary(abioticmodlambdaPLDE)
pldeabioticlambda <- r.squaredGLMM(abioticmodlambdaPLDE)[1,1]
summary(abioticmodlambdaPOLE)
poleabioticlambda <- r.squaredGLMM(abioticmodlambdaPOLE)[1,1]
summary(abioticmodlambdaTRCY)
trcyabioticlambda <- r.squaredGLMM(abioticmodlambdaTRCY)[1,1]
summary(abioticmodlambdaTROR)
trorabioticlambda <- r.squaredGLMM(abioticmodlambdaTROR)[1,1]
summary(abioticmodlambdaVERO)
veroabioticlambda <- r.squaredGLMM(abioticmodlambdaVERO)[1,1]

### Biotic only - hm, only have categorical neighbours.... 
#unless I could calculate high and low shannon? Consider this*
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survivalbioticlambda <- lmer(log_lambda_p1 ~ neighbours01 + (1|Site/Plot), 
                               data = filter(popdata, Species == specieslist[i]))
  print(summary(survivalbioticlambda))
}

for (i in 1:length(specieslist)){
  nam <- paste0("bioticmodlambda", specieslist[i])
  assign(nam, lmer(log_lambda_p1 ~ neighbours01 + (1|Site/Plot), 
                   data = filter(popdata, Species == specieslist[i])))
}

summary(bioticmodlambdaARCA)
arcabioticlambda <- r.squaredGLMM(bioticmodlambdaARCA)[1,1]
summary(bioticmodlambdaHYGL)
hyglbioticlambda <- r.squaredGLMM(bioticmodlambdaHYGL)[1,1]
summary(bioticmodlambdaLARO)
larobioticlambda <- r.squaredGLMM(bioticmodlambdaLARO)[1,1]
summary(bioticmodlambdaPEAI)
peaibioticlambda <- r.squaredGLMM(bioticmodlambdaPEAI)[1,1]
summary(bioticmodlambdaPLDE)
pldebioticlambda <- r.squaredGLMM(bioticmodlambdaPLDE)[1,1]
summary(bioticmodlambdaPOLE)
polebioticlambda <- r.squaredGLMM(bioticmodlambdaPOLE)[1,1]
summary(bioticmodlambdaTRCY)
trcybioticlambda <- r.squaredGLMM(bioticmodlambdaTRCY)[1,1]
summary(bioticmodlambdaTROR)
trorbioticlambda <- r.squaredGLMM(bioticmodlambdaTROR)[1,1]
summary(bioticmodlambdaVERO)
verobioticlambda <- r.squaredGLMM(bioticmodlambdaVERO)[1,1]

#Create table to put these values in
#create matrix with 2 columns filled with random value, 1
rsquaredtablelambda <- matrix(rep(1, times=2), ncol=2, nrow = 9)
#define column names and row names of matrix
colnames(rsquaredtablelambda) <- c('R squared biotic model', 'R squared abiotic model')
rownames(rsquaredtablelambda) <- c('ARCA', 'HYGL', 'LARO', 'PEAI', 'POLE', 'PLDE', 'TRCY', 'TROR', 'VERO')
#convert matrix to table 
rsquaredtablelambda <- as.data.frame(rsquaredtablelambda)
#view table 
rsquaredtablelambda
##Replace values with appropriate ones
#Biotic model
rsquaredtablelambda[1,1] <- arcabioticlambda
rsquaredtablelambda[2,1] <- hyglbioticlambda
rsquaredtablelambda[3,1] <- larobioticlambda
rsquaredtablelambda[4,1] <- peaibioticlambda
rsquaredtablelambda[5,1] <- pldebioticlambda
rsquaredtablelambda[6,1] <- polebioticlambda
rsquaredtablelambda[7,1] <- trcybioticlambda
rsquaredtablelambda[8,1] <- trorbioticlambda
rsquaredtablelambda[9,1] <- verobioticlambda
#Abiotic model
rsquaredtablelambda[1,2] <- arcaabioticlambda
rsquaredtablelambda[2,2] <- hyglabioticlambda
rsquaredtablelambda[3,2] <- laroabioticlambda
rsquaredtablelambda[4,2] <- peaiabioticlambda
rsquaredtablelambda[5,2] <- pldeabioticlambda
rsquaredtablelambda[6,2] <- poleabioticlambda
rsquaredtablelambda[7,2] <- trcyabioticlambda
rsquaredtablelambda[8,2] <- trorabioticlambda
rsquaredtablelambda[9,2] <- veroabioticlambda

#Export table as a csv file
write.csv(rsquaredtablelambda,"Output/Tables/lambda_rsquaredtable.csv")

### Merge them all into a nice table
### Binding survival, fecundity and lambda tables together ###
rsquaredtablemerged <- cbind(rsquaredtable, rsquaredtableseed, rsquaredtablelambda)
#Changing row and column names
rownames(rsquaredtablemerged) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
colnames(rsquaredtablemerged) <- c("biotic", "abiotic", "biotic", "abiotic", "biotic", "abiotic")

#Plotting with kableR
rsquaredtablemerged %>% kbl(caption = "<b>Supplementary X</b>. R squared values from abiotic and biotic factor only models for survival, fecundity and population growth rate.", digits = 2) %>%
  kable_classic(full_width = F, html_font = "Times") %>%
  column_spec(1, italic = T) %>%
  #row_spec(0, bold = T) %>%
  add_header_above(c("", "Survival" = 2, "Fecundity" = 2, "Lambda" = 2))

#Struggling to save this using save_kable and as_image() atm.
#Can copy it from the Viewer using copy to clipboard, maintain aspect ratio, first value 500

#### AIC Model Selection Abiotic vs Biotic ####
## Rerunning the models
#removing SDI/shannon
for (i in 1:length(specieslist.nop)){
  nam <- paste0("abioticmodsurv", specieslist.nop[i])
  assign(nam, glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                    family = binomial, data = filter(vitaldata, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("bioticmodsurv", specieslist.nop[i])
  assign(nam, glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                    family = binomial, data = filter(vitaldata, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("abioticmodseed", specieslist.nop[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = nbinom2, data = filter(seedmodeldata, Species == specieslist.nop[i])))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("bioticmodseed", specieslist.nop[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                      family = nbinom2, data = filter(seedmodeldata, Species == specieslist.nop[i])))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("abioticmodlambda", specieslist.nop[i])
  assign(nam, lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + (1|Site/Plot), 
                      data = filter(popdata, Species == specieslist.nop[i])))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("bioticmodlambda", specieslist.nop[i])
  assign(nam, lmer(log_lambda ~ Neighbours01 + (1|Site/Plot), 
                   data = filter(popdata, Species == specieslist.nop[i])))
}
## Abiotic and biotic together
for (i in 1:length(specieslist.nop)){
  nam <- paste0("bothmodsurv", specieslist.nop[i])
  assign(nam, glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                    family = binomial, data = filter(vitaldata, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("bothmodseed", specieslist.nop[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                      family = nbinom2, data = filter(seedmodeldata, Species == specieslist.nop[i])))
}
for (i in 1:length(specieslist)){
  nam <- paste0("bothmodlambda", specieslist[i])
  assign(nam, lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + (1|Site/Plot), 
                      data = filter(popdata, Species == specieslist[i])))
}
### Abiotic and biotic interacting
for (i in 1:length(specieslist.nop)){
  nam <- paste0("intmodsurv", specieslist.nop[i])
  assign(nam, glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                    family = binomial, data = filter(vitaldata, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}
for (i in 1:length(specieslist.nop)){
  nam <- paste0("intmodseed", specieslist.nop[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = nbinom2, data = filter(seedmodeldata, Species == specieslist.nop[i])))
}
for (i in 1:length(specieslist)){
  nam <- paste0("intmodlambda", specieslist[i])
  assign(nam, lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), 
                      data = filter(popdata, Species == specieslist[i])))
}

### Survival
a <- AIC(bioticmodsurvARCA, abioticmodsurvARCA, bothmodsurvARCA, intmodsurvARCA)
b <- AIC(bioticmodsurvHYGL, abioticmodsurvHYGL, bothmodsurvHYGL, intmodsurvHYGL)
c <- AIC(bioticmodsurvLARO, abioticmodsurvLARO, bothmodsurvLARO, intmodsurvLARO)
d <- AIC(bioticmodsurvPEAI, abioticmodsurvPEAI, bothmodsurvPEAI, intmodsurvPEAI)
e <- AIC(bioticmodsurvPLDE, abioticmodsurvPLDE, bothmodsurvPLDE, intmodsurvPLDE)
f <- AIC(bioticmodsurvTRCY, abioticmodsurvTRCY, bothmodsurvTRCY, intmodsurvTRCY)
g <- AIC(bioticmodsurvTROR, abioticmodsurvTROR, bothmodsurvTROR, intmodsurvTROR)
h <- AIC(bioticmodsurvVERO, abioticmodsurvVERO, bothmodsurvVERO, intmodsurvVERO)
### Fecundity
aseed <- AIC(bioticmodseedARCA, abioticmodseedARCA, bothmodseedARCA, intmodseedARCA)
bseed <- AIC(bioticmodseedHYGL, abioticmodseedHYGL, bothmodseedHYGL, intmodseedHYGL)
cseed <- AIC(bioticmodseedLARO, abioticmodseedLARO, bothmodseedLARO, intmodseedLARO)
dseed <- AIC(bioticmodseedPEAI, abioticmodseedPEAI, bothmodseedPEAI, intmodseedPEAI)
eseed <- AIC(bioticmodseedPLDE, abioticmodseedPLDE, bothmodseedPLDE, intmodseedPLDE)
fseed <- AIC(bioticmodseedTRCY, abioticmodseedTRCY, bothmodseedTRCY, intmodseedTRCY)
gseed <- AIC(bioticmodseedTROR, abioticmodseedTROR, bothmodseedTROR, intmodseedTROR)
hseed <- AIC(bioticmodseedVERO, abioticmodseedVERO, bothmodseedVERO, intmodseedVERO)
### Lambda
alambda <- AIC(bioticmodlambdaARCA, abioticmodlambdaARCA, bothmodlambdaARCA, intmodlambdaARCA)
blambda <- AIC(bioticmodlambdaHYGL, abioticmodlambdaHYGL, bothmodlambdaHYGL, intmodlambdaHYGL)
clambda <- AIC(bioticmodlambdaLARO, abioticmodlambdaLARO, bothmodlambdaLARO, intmodlambdaLARO)
dlambda <- AIC(bioticmodlambdaPEAI, abioticmodlambdaPEAI, bothmodlambdaPEAI, intmodlambdaPEAI)
elambda <- AIC(bioticmodlambdaPLDE, abioticmodlambdaPLDE, bothmodlambdaPLDE, intmodlambdaPLDE)
flambda <- AIC(bioticmodlambdaTRCY, abioticmodlambdaTRCY, bothmodlambdaTRCY, intmodlambdaTRCY)
glambda <- AIC(bioticmodlambdaTROR, abioticmodlambdaTROR, bothmodlambdaTROR, intmodlambdaTROR)
hlambda <- AIC(bioticmodlambdaVERO, abioticmodlambdaVERO, bothmodlambdaVERO, intmodlambdaVERO)

## Presenting them in tables
### Survival
surv_AIC <- cbind(a, b[,2], c[,2], d[,2], e[,2], f[,2], g[,2], h[,2])
colnames(surv_AIC) <- c('df', "Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
rownames(surv_AIC) <- c('Biotic', 'Abiotic', 'Both', 'Interaction')
surv_AIC <- surv_AIC %>% rownames_to_column("Model")
#Plotting with kableR
surv_AIC %>%
  kbl(align = 'lccccccccc', digits = 1) %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("Survival" = 10), align = "l", background = "lightgrey") %>%
  row_spec(0, italic = T)
  #column_spec(4, bold = min(surv_AIC[4]))
  #cell_spec(surv_AIC[1,5], bold = T)
#Can't figure out conditional bolding
#mutate(`Arctotheca calendula` = cell_spec(`Arcotheca calendula`, bold = ifelse(min(surv_AIC$`Arctotheca calendula`), TRUE, FALSE)))

### Fecundity
seed_AIC <- cbind(aseed, bseed[,2], cseed[,2], dseed[,2], eseed[,2], fseed[,2], gseed[,2], hseed[,2])
colnames(seed_AIC) <- c('df', "Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
rownames(seed_AIC) <- c('Biotic', 'Abiotic', 'Both', 'Interaction')
seed_AIC <- seed_AIC %>% rownames_to_column("Model")
seed_AIC %>%
  kbl(align = 'lccccccccc', digits = 1) %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("Seed production" = 10), align = "l", background = "lightgrey")%>%
  row_spec(0, italic = T)

### Lambda
lambda_AIC <- cbind(alambda, blambda[,2], clambda[,2], dlambda[,2], elambda[,2], flambda[,2], glambda[,2], hlambda[,2])
colnames(lambda_AIC) <- c('df', "Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
rownames(lambda_AIC) <- c('Biotic', 'Abiotic', 'Both', 'Interaction')
lambda_AIC <- lambda_AIC %>% rownames_to_column("Model")
lambda_AIC %>%
  kbl(align = 'lccccccccc', digits = 1) %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("Population growth" = 10), align = "l", background = "lightgrey")%>%
  row_spec(0, italic = T)
#resize 600

### Looking at distribution of data - linear or quadratic? ####
#Meeting with John 13/03
#Watering*N*Cover
## With quadratic
# surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2)
#x is std_PC1
#cbind // plogis // 1 - int, x for PC1, x^2
# (1, x, x^2)
# for loops
species.list.s<-list(arcadata, hygldata, larodata, peaidata, pldedata, poledata, trcydata, trordata, verodata)
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

## Germination ##
arcasimplemod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, arcadata)
summary(arcasimplemod)
plot(percent_germ ~ std_PC1, arcadata)
curve(exp(cbind(1, x, x^2)%*%fixef(arcasimplemod)), add = T)
title(main = "ARCA")

# Making a bunch of plots including quadratic term to compare line fits
### Survival, PC2
dev.off()
pdf("Output/Figures/germ_PC2_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC2, pch=19, col="grey60", ylab="Percentage of seeds that germinated", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(cbind(total_germ, total_no_germ)~std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

## Survival ##
#ARCA
arcasimplemod <- glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, arcadata)
summary(arcasimplemod)
plot(surv_to_produce_seeds ~ std_PC1, arcadata)
curve(exp(cbind(1, x, x^2)%*%fixef(arcasimplemod)), add = T)
title(main = "ARCA")

#Loop all species
for (i in 1:length(specieslist)){
  print(specieslist[i])
  mod <- glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                  family = binomial, data = filter(datanonly, Species == specieslist[i]))
  print(summary(mod))
  plot(surv_to_produce_seeds ~ std_PC1, data = filter(vitaldata, Species == specieslist[i]))
  curve(exp(cbind(1, x, x^2)%*%fixef(mod)), add = T)
  title(main = specieslist[i])
}
### Survival, PC1
dev.off()
pdf("Output/Figures/Surv_PC1_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC1, pch=19, col="grey60", ylab="Probability of survival", xlab="PC1", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

### Survival, PC2
dev.off()
pdf("Output/Figures/Surv_PC2_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC2, pch=19, col="grey60", ylab="Probability of survival", xlab="PC2", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(surv_to_produce_seeds ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC2)
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

### Survival, neighbour abundance
dev.off()
pdf("Output/Figures/Surv_NA_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_logp1_totalabund, pch=19, col="grey60", ylab="Probability of survival", xlab="Total neighbour abundance (log plus 1, standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(surv_to_produce_seeds ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_logp1_totalabund)
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()


## Seed production ##
##Seed production ~ PC1
#ARCA
arcasimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), arcadata)
summary(arcasimplemod)
plot(No_viable_seeds_grouped ~ std_PC1, arcadata)
curve(exp(cbind(1, x, x^2)%*%fixef(arcasimplemod)), add = T)
title(main = "ARCA")

#PEAI
peaisimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), peaidata)
summary(peaisimplemod)
plot(No_viable_seeds_grouped ~ std_PC1, peaidata)
curve(exp(cbind(1, x, x^2)%*%fixef(peaisimplemod)), add = T)
title(main = "PEAI")

# With a loop
#Need to use glmmTMB, getting model convergence issues using glmer.nb (2) and optimiser (1)
#Species list ARCA, HYGL, etc.
for (i in 1:length(specieslist)){
  print(specieslist[i])
  mod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                 data = filter(datanonly, Species == specieslist[i]))
  print(summary(mod))
  plot(No_viable_seeds_grouped ~ std_PC1, data = filter(vitaldata, Species == specieslist[i]))
  curve(exp(cbind(1, x, x^2)%*%fixef(mod)), add = T)
  title(main = specieslist[i])
}

#One didn't converge, figuring out which
arcasimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), arcadata)
hyglsimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), hygldata)
larosimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), larodata)
peaisimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), peaidata)
pldesimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), pldedata)
polesimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), poledata)
trcysimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), trcydata)
trorsimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), trordata)
#TROR not converging
trorsimplemod <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, trordata)
verosimplemod <- glmer.nb(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), verodata)
# VERO not converging
verosimplemod <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, trordata)
ggplot(trordata, aes(x = std_PC1, y = No_viable_seeds_grouped))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

### Seed production, PC1
dev.off()
pdf("Output/Figures/SP_PC1_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_PC1, pch=19, col="grey60", ylab="Number of viable seeds", xlab="PC1", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()
#Check what is happening with POLE, confidence interval is entire plot
polesimplemod <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, poledata)
summary(polesimplemod)
polesimplemoddharma <- simulateResiduals(polesimplemod)
plot(polesimplemoddharma)
#very bad residual plot

### Seed production, PC2
dev.off()
pdf("Output/Figures/SP_PC2_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_PC2, pch=19, col="grey60", ylab="Number of viable seeds", xlab="PC2", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmmTMB(No_viable_seeds_grouped ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC2)
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

### Seed production, neighbour abundance
dev.off()
pdf("Output/Figures/SP_NA_quadratic.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_logp1_totalabund, pch=19, col="grey60", ylab="Number of viable seeds", xlab="Total neighbour abundance (standardised, log plus 1)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_logp1_totalabund)
  preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()
#What is happening with pole?
polesimplemod2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = nbinom2, poledata)
summary(polesimplemod2)
polesimplemoddharma2 <- simulateResiduals(polesimplemod2)
plot(polesimplemoddharma2)
#bad dharma

### Checking model fits and vif ####
## Checking for multicollinearity in models using variance inflation factor
# VIF over 5 is a problem
for (i in 1:length(specieslist.nop)){
  nam <- paste0("SDIsurvivalmod", specieslist[i])
  assign(nam, glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + shannon + (1|Site/Plot),
                    family = binomial, data = filter(datanonly, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}

SDIsurvivalmodARCA <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                             family = binomial, filter(datanonly, Species == 'ARCA'))
arcaSDIsurvdharma <- simulateResiduals(SDIsurvivalmodARCA)
plot(arcaSDIsurvdharma)
summary(SDIsurvivalmodARCA)
vif(SDIsurvivalmodARCA)

SDIsurvivalmodHYGL <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'HYGL'))
#test <- datanonly %>% filter(Species == 'HYGL')
#5 B HYGL T 1 didn't germinate or have survival data but has 3 neighbours recorded, check this*
datanonly %>% filter(Species == 'HYGL') %>%
  ggplot(aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme
hyglSDIsurvdharma <- simulateResiduals(SDIsurvivalmodHYGL)
plot(hyglSDIsurvdharma)
summary(SDIsurvivalmodHYGL)
vif(SDIsurvivalmodHYGL)

SDIsurvivalmodLARO <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'LARO'))
laroSDIsurvdharma <- simulateResiduals(SDIsurvivalmodLARO)
plot(laroSDIsurvdharma)
summary(SDIsurvivalmodLARO)
vif(SDIsurvivalmodLARO)

#Only converges with optimiser
SDIsurvivalmodPEAI <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'PEAI'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
datanonly %>% filter(Species == 'PEAI') %>%
  ggplot(aes(x = Treatment, y = surv_to_produce_seeds))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme
PEAISDIsurvdharma <- simulateResiduals(SDIsurvivalmodPEAI)
plot(PEAISDIsurvdharma)
#dharma bit unhappy
summary(SDIsurvivalmodPEAI)
vif(SDIsurvivalmodPEAI)
#vif fine

SDIsurvivalmodPLDE <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'PLDE'))
PLDESDIsurvdharma <- simulateResiduals(SDIsurvivalmodPLDE)
plot(PLDESDIsurvdharma)
summary(SDIsurvivalmodPLDE)
vif(SDIsurvivalmodPLDE)
#I think this is okay - looking at GVIF^(1/2*Df)) PC1??

SDIsurvivalmodTRCY <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'TRCY'))
TRCYSDIsurvdharma <- simulateResiduals(SDIsurvivalmodTRCY)
plot(TRCYSDIsurvdharma)
summary(SDIsurvivalmodTRCY)
vif(SDIsurvivalmodTRCY)

SDIsurvivalmodTROR <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'TROR'))
TRORSDIsurvdharma <- simulateResiduals(SDIsurvivalmodTROR)
plot(TRORSDIsurvdharma)
summary(SDIsurvivalmodTROR)
vif(SDIsurvivalmodTROR)

SDIsurvivalmodVERO <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI + (1|Site/Plot),
                            family = binomial, filter(datanonly, Species == 'VERO'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
VEROSDIsurvdharma <- simulateResiduals(SDIsurvivalmodVERO)
plot(VEROSDIsurvdharma)
#dharma unhappy
datanonly %>% filter(Species == 'VERO') %>%
  ggplot(aes(x = std_PC1, y = surv_to_produce_seeds))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme
# All plants with low PC1 survived, everything else seems fine
summary(SDIsurvivalmodVERO)
vif(SDIsurvivalmodVERO)

#### Is competition stronger where species prefer (abiotic conditions)? ####
#Need a dataset for no neighbours only.
data_no_nbh_vital <- vitaldata %>% filter(Total_abundance == 0)
data_no_nbh_seed <- seedmodeldata %>% filter(Total_abundance == 0)
data_no_nbh_lambda <- popdata %>% filter(Neighbours01 == 'Neighbours0')

ggplot(data_no_nbh_lambda, aes(x = std_PC1, y = log_lambda))+
  geom_point(alpha = 0.4)+
  theme_classic()+
  geom_smooth(method='lm')+
  facet_wrap(~Species)

ggplot(popdata, aes(x = std_PC1, y = log_lambda, colour = Neighbours01))+
  geom_point(alpha = 0.4)+
  theme_classic()+
  geom_smooth(method='lm')+
  facet_wrap(~Species)

### PC1 lambda coloured by neighbours or no
#Can't have Site/Plot
#Lambda - vero example
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=ifelse(Neighbours01=='Neighbours1', alpha("forestgreen", 0.3), alpha("purple", 0.3)), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#neighbours
model2<-lmer(log_lambda~std_PC1 + (1|Site), filter(lambdavero, Neighbours01=='Neighbours1'))
preddata2 <- with(model2, data.frame(1, x_to_plot))
plotted.pred2 <- glmm.predict(mod = model2, newdat = preddata2, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred2$y, upper = plotted.pred2$upper, lower = plotted.pred2$lower, env.colour = "forestgreen", env.trans = 50, line.colour = "forestgreen", line.weight = 2, line.type = 1)
#no neighbours
model3<-lmer(log_lambda~std_PC1 + (1|Site), filter(lambdavero, Neighbours01=='Neighbours0'))
preddata3 <- with(model3, data.frame(1, x_to_plot))
plotted.pred3 <- glmm.predict(mod = model3, newdat = preddata3, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred3$y, upper = plotted.pred3$upper, lower = plotted.pred3$lower, env.colour = "purple", env.trans = 50, line.colour = "purple", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 6, col = "red")

#No pole
#arca - PC2
#hygl - wet
#peai - dry, wet, pc2
# no significance for others incl. vero. Interesting! Not actually preference for high PC1
for (i in 1:length(specieslist.nop)){
  print(specieslist.nop[i])
  no_nbh_mod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Treatment:std_PC1 + (1|Site), 
                       data = filter(data_no_nbh_lambda, Species == specieslist.nop[i]))
  print(summary(no_nbh_mod))
}

### How about for fecundity?
#plde and vero maybe evidence of biotic interactions dampening 
#performance in preferred environment
ggplot(vitaldata, aes(x = std_PC1, y = log(No_viable_seeds_grouped)+1, colour = Neighbours01))+
  geom_point(alpha = 0.4)+
  theme_classic()+
  geom_smooth(method='lm')+
  facet_wrap(~Species)

#peai pc1, vero pc1 - neither had interaction of PC1:neighbours so nope
for (i in 1:length(specieslist.nop)){
  print(specieslist.nop[i])
  no_nbh_mod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + Dodder01 + Treatment:std_PC1 + (1|Site), 
                     family = nbinom2, data = filter(data_no_nbh_seed, Species == specieslist.nop[i]))
  print(summary(no_nbh_mod))
}

## Survival
ggplot(vitaldata, aes(x = std_PC1, y = surv_to_produce_seeds, colour = Neighbours01))+
  geom_point(alpha = 0.4)+
  theme_classic()+
  geom_smooth(method='lm')+
  facet_wrap(~Species)

#none significant for PC1
for (i in 1:length(specieslist.nop)){
  print(specieslist.nop[i])
  no_nbh_mod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + Dodder01 + Treatment:std_PC1 + (1|Site), 
                        family = binomial, data = filter(data_no_nbh_vital, Species == specieslist.nop[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(no_nbh_mod))
}


### Need to update below here* ####
### Viable seed production ###
for (i in 1:length(specieslist)){
  print(specieslist[i])
  fecundity <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI +
                         (1|Site/Plot), family = nbinom2, data = filter(datanonly, Species == specieslist[i]))
  print(summary(fecundity))
}
#POLE significant negative response to SDI
# TRCY positive response to SDI

for (i in 1:length(specieslist)){
  nam <- paste0("SDIfecunditymod", specieslist[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + SDI +
                        (1|Site/Plot), family = nbinom2, data = filter(datanonly, Species == specieslist[i])))
}

arcaSDIfecunditydharma <- simulateResiduals(SDIfecunditymodARCA)
plot(arcaSDIfecunditydharma)
summary(SDIfecunditymodARCA)
vif(SDIfecunditymodARCA)
hyglSDIfecunditydharma <- simulateResiduals(SDIfecunditymodHYGL)
plot(hyglSDIfecunditydharma)
summary(SDIfecunditymodHYGL)
laroSDIfecunditydharma <- simulateResiduals(SDIfecunditymodLARO)
plot(laroSDIfecunditydharma)
summary(SDIfecunditymodLARO)
PEAISDIfecunditydharma <- simulateResiduals(SDIfecunditymodPEAI)
plot(PEAISDIfecunditydharma)
summary(SDIfecunditymodPEAI)
PLDESDIfecunditydharma <- simulateResiduals(SDIfecunditymodPLDE)
plot(PLDESDIfecunditydharma)
summary(SDIfecunditymodPLDE)
#PLDE slight issues
POLESDIfecunditydharma <- simulateResiduals(SDIfecunditymodPOLE)
plot(POLESDIfecunditydharma)
summary(SDIfecunditymodPOLE)
#POLE big issues
TRCYSDIfecunditydharma <- simulateResiduals(SDIfecunditymodTRCY)
plot(TRCYSDIfecunditydharma)
summary(SDIfecunditymodTRCY)
TRORSDIfecunditydharma <- simulateResiduals(SDIfecunditymodTROR)
plot(TRORSDIfecunditydharma)
summary(SDIfecunditymodTROR)
VEROSDIfecunditydharma <- simulateResiduals(SDIfecunditymodVERO)
plot(VEROSDIfecunditydharma)
summary(SDIfecunditymodVERO)

## Plots for TRCY
## Try plotting both of these from the models using lab retreat code
### Something funky happening here with y axis**
seedtrcy %>% filter(Total_abundance > 0) %>%
ggplot(aes(x = SDI, y = surv_to_produce_seeds))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()+
  my_theme

### Are PC1 and PC2 correlated? ####
#Only want one data point per plot
plotvar <- vitaldata %>% group_by(Site, Plot) %>% filter(row_number() == 1)
ggplot(plotvar, aes(x = std_PC1, y = std_PC2))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth()+
  theme_classic()
modpc1pc2 <- lm(std_PC2 ~ std_PC1, plotvar)
summary(modpc1pc2)
#### Question 1 relative importance A. and B. factors####
##What is the relative importance of abiotic and biotic factors for survival and fecundity?
#Additive model: response ~ total_abundance + water + PC1 + PC2 + PC3 + RE
### Q1 - Survival to produce seeds ####

# There are 129 subplots that do not have total_abundance information but did germinate
# And 82 subplots that weren't surveyed but germinated...?

specieslistnopole <- specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "TRCY", "TROR", "VERO")


for (i in 1:length(specieslistnopole)){
  print(specieslist[i])
  survival <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + (1|Site/Plot),
                  family = binomial, data = filter(vitaldata, Species == specieslistnopole[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(survival))
}

#### Putting models into a list for coef plotting ####
# GERMINATION
germmodels <- list()
for (i in 1:length(specieslistnopole)){
  germmodels[[i]] <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                               family = binomial, data = filter(vitaldata, Species == specieslistnopole[i]))
}
# , control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))

#Making coefficient plot for survival as a function of environment
plot_models(germmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

# SURVIVAL
survivalmodels <- list()
for (i in 1:length(specieslistnopole)){
  survivalmodels[[i]] <- glmer(surv_to_produce_seeds ~ std_PC1 + std_PC2 + Treatment + std_logp1_totalabund + Dodder01 +
                                 (1|Site/Plot), family = binomial, data = filter(vitaldata, Species == specieslistnopole[i]))
}
#Making coefficient plot for survival as a function of environment
plot_models(survivalmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
 scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
### Assigning model names by species to plot
for (i in 1:length(specieslistnopole)){
  nam <- paste0("survmod", specieslistnopole[i])
  assign(nam, glmer(surv_to_produce_seeds ~ std_PC1 + std_PC2 + Treatment + std_logp1_totalabund + Dodder01 +
                      (1|Site/Plot), family = "binomial", data = filter(vitaldata, Species == specieslistnopole[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
}

#ARCA survival model, model created in loop above
arcasurvdharma <- simulateResiduals(survmodARCA)
plot(arcasurvdharma)
summary(survmodARCA)
hyglsurvdharma <- simulateResiduals(survmodHYGL)
plot(hyglsurvdharma)
summary(survmodHYGL)
larosurvdharma <- simulateResiduals(survmodLARO)
plot(larosurvdharma)
summary(survmodLARO)
PEAIsurvdharma <- simulateResiduals(survmodPEAI)
plot(PEAIsurvdharma)
summary(survmodPEAI)
PLDEsurvdharma <- simulateResiduals(survmodPLDE)
plot(PLDEsurvdharma)
summary(survmodPOLE)
TRCYsurvdharma <- simulateResiduals(survmodTRCY)
plot(TRCYsurvdharma)
summary(survmodTRCY)
TRORsurvdharma <- simulateResiduals(survmodTROR)
plot(TRORsurvdharma)
summary(survmodTROR)
VEROsurvdharma <- simulateResiduals(survmodVERO)
plot(VEROsurvdharma)
summary(survmodVERO)

#Putting models in a list for tab_model and coef plots
survmodels <- list()
survmodels[[1]] <- survmodARCA
survmodels[[2]] <- survmodHYGL
survmodels[[3]] <- survmodLARO
survmodels[[4]] <- survmodPEAI
survmodels[[5]] <- survmodPLDE
survmodels[[6]] <- survmodTRCY
survmodels[[7]] <- survmodTROR
survmodels[[8]] <- survmodVERO
tab_model(survmodels, transform = NULL)
plot_models(survmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

# FECUNDITY
fecunditymodels <- list()
for (i in 1:length(specieslistnopole)){
  fecunditymodels[[i]] <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + std_PC2 + Treatment + std_logp1_totalabund + Dodder01 +
                                    (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslistnopole[i]))
}
#Making coefficient plot for survival as a function of environment
plot_models(fecunditymodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

# LAMBDA
lambdamodels <- list()
for (i in 1:length(specieslistnopole)){
  lambdamodels[[i]] <- lmer(log_lambda_p1 ~ std_PC1 + std_PC2 + Treatment + neighbours01 + (1|Site/Plot), 
                               data = filter(popdata, Species == specieslistnopole[i]))
}
#Making coefficient plot for survival as a function of environment
plot_models(lambdamodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

# #Plots of significant with neighbour abundance (ARCA and TROR) -- haven't updated
# par(mfrow=c(2,1))
# x_to_plot<-seq.func(vitaldata$std_logp1_totalabund)
# with(arcadata, plot(surv_to_produce_seeds ~ std_logp1_totalabund))
# arcapreddata <- with(arcadata, data.frame(1, x_to_plot, 0, 0, 0, 0, 0, 0))
# arcapred <- glmm.predict(mod = survmodARCA, newdat = arcapreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
# plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)
# with(trordata, plot(surv_to_produce_seeds ~ std_logp1_totalabund))
# trorpreddata <- with(trordata, data.frame(1, x_to_plot, 0, 0, 0, 0, 0))
# trorpred <- glmm.predict(mod = trorsurvmod2, newdat = trorpreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
# plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

### Plotting with all species together ####
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
dev.off()
pdf("Output/Figures/survival_PC1.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
species.list.s<-list(arcadata, hygldata, larodata, peaidata, pldedata, poledata, trcydata, trordata, verodata)
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC1, pch=19, col="grey60", ylab="Survival rate", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "PC1", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                 (1|Site/Plot), family = binomial, data = plotted.data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, 0, 0, 0, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

dev.off()
pdf("Output/Figures/survival_PC2.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
species.list.s<-list(arcadata, hygldata, larodata, peaidata, pldedata, poledata, trcydata, trordata, verodata)
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC2, pch=19, col="grey60", ylab="Survival rate", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "PC1", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                 (1|Site/Plot), family = binomial, data = plotted.data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  x_to_plot<-seq.func(plotted.data$std_PC2)
  preddata <- with(model, data.frame(1, 0, 0, 0, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

#Simpler version for ESA presentation:
dev.off()
pdf("Output/Figures/survival_PC1_ESA.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(2,2,6,0.5))
#Margins: bottom, left, top, right
par(pty="s")
species.list.s<-list(arcadata, hygldata, larodata, peaidata, pldedata, poledata, trcydata, trordata, verodata)
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC1, pch=19, col="grey60", ylab="", xlab="", cex = 3, cex.lab=3, cex.axis=3,tck=-0.02, alpha = 0.4, axes = 'FALSE', frame.plot = TRUE)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=4)
  model<-glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                 (1|Site/Plot), family = binomial, data = plotted.data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, 0, 0, 0, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, 
               env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 4, line.type = 1)
  Axis(side=1, labels=TRUE)
  Axis(side=2, labels=TRUE)
}
dev.off()

### Q1 - Viable seed production #####
# Should I use No_viable_seeds_grouped or seeds_percent???
#Using No_viable_seeds_grouped for individual species models
#and seeds_percent for trait models

# optimiser 
#control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))

specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  fecundity <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                      (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(fecundity))
}
#Assigning model names by species
for (i in 1:length(specieslist)){
  nam <- paste0("seedmod", specieslist[i])
  assign(nam, glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                      (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i])))
}

#Checking residuals
arcaseeddharma <- simulateResiduals(seedmodARCA)
plot(arcaseeddharma)
summary(seedmodARCA)
hyglseeddharma <- simulateResiduals(seedmodHYGL)
plot(hyglseeddharma)
summary(seedmodHYGL)
laroseeddharma <- simulateResiduals(seedmodLARO)
plot(laroseeddharma)
summary(seedmodLARO)
PEAIseeddharma <- simulateResiduals(seedmodPEAI)
plot(PEAIseeddharma)
summary(seedmodPEAI)
PLDEseeddharma <- simulateResiduals(seedmodPLDE)
plot(PLDEseeddharma)
summary(seedmodPLDE)
POLEseeddharma <- simulateResiduals(seedmodPOLE)
plot(POLEseeddharma)
summary(seedmodPOLE)
# POLE unhappy. Need to try simpler model 
#or don't run it for this vital rate. Just not enough data??
# Only one occurrence of Dodder. Not a lot of neighbour data, but nor does ARCA
testDispersion(POLEseeddharma)
#Underdispersed! 
testZeroInflation(POLEseeddharma)
#Not zero inflated
#Trying a Conway-Maxwell-Poisson distribution (good for underdispersion)
#Needs to be within glmmTMB to have random effects
library(mpcmp)
poleseedmod3 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                          (1|Site/Plot), family = compois, data = seedpole)
summary(poleseedmod3)
POLEseeddharma3 <- simulateResiduals(poleseedmod3)
plot(POLEseeddharma3)
#MUCH better residuals
#Testing it for plde - very similar result!! Same residuals
#pldeseedmod3 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
#                          (1|Site/Plot), family = compois, data = seedplde)
#summary(pldeseedmod3)
summary(seedmodPLDE)
poleseedmod2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                          (1|Site/Plot), family = nbinom2, data = seedpole)
summary(poleseedmod2)
with(seedpole, plot(No_viable_seeds_grouped ~ std_logp1_totalabund))
POLEseeddharma2 <- simulateResiduals(poleseedmod2)
plot(POLEseeddharma2)
#
TRCYseeddharma <- simulateResiduals(seedmodTRCY)
plot(TRCYseeddharma)
summary(seedmodTRCY)
# trcyseedmod2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
#                           (1|Site/Plot), family = nbinom2, data = seedtrcy)
#summary(trcyseedmod2)
TRORseeddharma <- simulateResiduals(seedmodTROR)
plot(TRORseeddharma)
summary(seedmodTROR)
VEROseeddharma <- simulateResiduals(seedmodVERO)
plot(VEROseeddharma)
summary(seedmodVERO)

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

#Trying poleseedmod3 plot, comparing to glmmTMB plot
summary(poleseedmod3)
x_to_plot<-seq.func(seedpole$std_PC1)
with(seedpole, plot(No_viable_seeds_grouped ~ seedpole$std_PC1))
polepreddata <- with(seedpole, data.frame(1, 0, 0, 0, x_to_plot, 0, 0))
polepred <- glmm.predict(mod = poleseedmod3, newdat = polepreddata, se.mult = 1.96, logit_link = FALSE, log_link = TRUE, glmmTMB = TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = polepred$y, upper = polepred$upper, lower = polepred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

######### Plotting everything as PDF - fecundity ###
###create a list of names for figure headings
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

dev.off()
pdf("Output/Figures/fecundity_PC1", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
species.list.f<-list(seedarca, seedhygl, seedlaro, seedpeai, seedplde, seedpole, seedtrcy, seedtror, seedvero)
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$seeds_percent~plotted.data$std_PC1, pch=19, col="grey60", ylab="Number seeds", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "PC1", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmmTMB(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                   (1|Site/Plot), family = nbinom2, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, 0, 0, 0, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

#Creating separate plots for TRCY (won't converge unless with percent seeds) and POLE
#arcaseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 +
#                          (1|Site/Plot), family = nbinom2, data = seedarca)
#arcaseedmod1test <- simulateResiduals(arcaseedmod1)
#plot(arcaseedmod1test)
#summary(arcaseedmod1)
#x_to_plot<-seq.func(seedarca$std_logp1_totalabund)
#with(seedarca, plot(No_viable_seeds_grouped ~ seedarca$std_logp1_totalabund))
#arcapreddata <- with(seedarca, data.frame(1, x_to_plot, 0, 0, 0))
#arcapred <- glmm.predict(mod = arcaseedmod1, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
#plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

#Simpler version for ESA presentation:
#Points for seeds_percent and no_viable seeds in same place!
dev.off()
pdf("Output/Figures/fecundity_PC1_ESA", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(2,2,6,0.5))
#Margins: bottom, left, top, right
par(pty="s")
species.list.f<-list(seedarca, seedhygl, seedlaro, seedpeai, seedplde, seedpole, seedtrcy, seedtror, seedvero)
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$seeds_percent~plotted.data$std_PC1, pch=19, col="grey60", ylab="", xlab="", cex.lab=3, cex.axis=2.00,tck=-0.02, axes = 'FALSE', frame.plot = TRUE)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=4)
  model<-glmmTMB(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                   (1|Site/Plot), family = nbinom2, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, 0, 0, 0, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, 
               env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 4, line.type = 1)
  }
dev.off()

##### Question 2 traits explaining responses ####
#Trait data needed for plotting
#Creating a dataframe with SLA data - need the as.data.frame part!
speciessla <- seedmodeldata %>% select("Species", "std_log_SLA") %>% group_by(Species) %>% filter(row_number() == 1)
speciessla <- as.data.frame(speciessla)

#Creating a dataframe with D13C data
speciesD13C <- seedmodeldata %>% select("Species", "std_log_D13C") %>% group_by(Species) %>% filter(row_number() == 1)
speciesD13C <- as.data.frame(speciesD13C)

### John's working 07/11/2021 plotting simple model sla modulating response of fecundity to PC2 ###
seedsla2 <- glmer.nb(seeds_percent ~ std_PC2 + std_log_SLA + std_PC2:std_log_SLA +
                       (1|Site/Plot) + (1+std_PC2|Species), family = nbinom2, seedmodeldata)
summary(seedsla2)
#My interpretation: the slopes are the slope responses of PC2 for each species plus
#the interaction of PC2:SLA*each species' SLA value.
#Unique value for PC2 slopes, one value for PC2:SLA, unique values of SLA
species_slopes<-coef(seedsla2)$Species[,2] + coef(seedsla2)$Species[,4]*speciessla$std_log_SLA

plot(species_slopes ~ speciessla$std_log_SLA)
curve(cbind(1,x)%*%fixef(seedsla2)[c(2,4)], add=T)
abline(h=0, lty=3)

##Does SLA modulate the responses of survival to PC1, PC2, neighbour abundance or treatment?
### Survival models ###
#SLA
survsla <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_SLA + 
                      std_PC1:std_log_SLA + std_logp1_totalabund:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
                      (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, vitaldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
survsladharma <- simulateResiduals(survsla)
plot(survsladharma)
#good
summary(survsla)
#No interactions significant

#LDMC
survldmc <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_LDMC + 
                       std_PC1:std_LDMC + std_logp1_totalabund:std_LDMC + std_PC2:std_LDMC + Treatment:std_LDMC +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, vitaldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
survldmcdharma <- simulateResiduals(survldmc)
#good but outlier test significant?
plot(survldmcdharma)
summary(survldmc)
#No interactions significant

#D13C
survD13C <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_D13C + 
                       std_PC1:std_log_D13C + std_logp1_totalabund:std_log_D13C + std_PC2:std_log_D13C + Treatment:std_log_D13C +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, vitaldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
survD13Cdharma <- simulateResiduals(survD13C)
plot(survD13Cdharma)
#good
summary(survD13C)
#No interactions significant

### Fecundity models ###
#SLA
#control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
#Converges with glmer.nb instead of glmmTMB!
seedsla <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_SLA +
                      std_PC1:std_log_SLA + std_logp1_totalabund:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
                      (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, seedmodeldata)
## Saving the model to load later
save(seedsla, file = "Output/Models/seedslamodel.RData")
seedsladharma <- simulateResiduals(seedsla)
plot(seedsladharma)
#okay?
summary(seedsla)
#PC1 and PC2 interactions significant

#LDMC
#Will only converge with optimiser
seedldmc <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + std_LDMC + 
                       std_PC1:std_LDMC + std_logp1_totalabund:std_LDMC + std_PC2:std_LDMC + Treatment:std_LDMC +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, 
                     seedmodeldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
## Saving the model to load later
save(seedldmc, file = "Output/Models/seedldmcmodel.RData")
seedldmcdharma <- simulateResiduals(seedldmc)
plot(seedldmcdharma)
#good
summary(seedldmc)
#no significant interactions

#D13C
seedD13C <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + std_log_D13C + 
                       std_PC1:std_log_D13C + std_logp1_totalabund:std_log_D13C + std_PC2:std_log_D13C + Treatment:std_log_D13C +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, seedmodeldata,
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
## Saving the model to load later
save(seedD13C, file = "Output/Models/seedD13Cmodel.RData")
seedD13Cdharma <- simulateResiduals(seedD13C)
plot(seedD13Cdharma)
#Okay, KS test significant
summary(seedD13C)
#No significant interactions

### Lambda models ###
#Trying simple model first - SLA modulating reponse to PC1
lambdaslasimple <- lmer(log_lambda_p1 ~ std_PC1 + std_log_SLA + std_PC1:std_log_SLA + (1|Site/Plot) + (std_PC1|Species), popdata)
lambdaslasimpledharma <- simulateResiduals(lambdaslasimple)
plot(lambdaslasimpledharma)
#Okay, KS test significant
summary(lambdaslasimple)

#Simple model with PC1 and Treatment --- won't converge, treatment seems to be the problem
lambdaslasimple2 <- lmer(log_lambda_p1 ~ std_PC1 + Treatment + std_log_SLA + std_PC1:std_log_SLA + Treatment:std_log_SLA + (1|Site/Plot) + (Treatment + std_PC1|Species), popdata)
lambdaslasimpledharma <- simulateResiduals(lambdaslasimple2)
plot(lambdaslasimpledharma)
#Okay, KS test significant
summary(lambdaslasimple2)

#Simple model just treatment
lambdaslasimple4 <- lmer(log_lambda_p1 ~ Treatment + std_log_SLA + Treatment:std_log_SLA + (1|Site/Plot) + (Treatment|Species), popdata)
lambdaslasimpledharma <- simulateResiduals(lambdaslasimple4)
plot(lambdaslasimpledharma)
#Okay, KS test significant
summary(lambdaslasimple4)
#Closest to significant but never is

##SLA - not converging with Treatment, so removing it
#lambdasla <- lmer(log_lambda_p1 ~ neighbours01 + Treatment + std_PC1 + std_PC2 + std_log_SLA +
#                        std_PC1:std_log_SLA + neighbours01:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
#                        (1|Site/Plot) + (std_PC1 + std_PC2 + neighbours01 + Treatment|Species), popdata)
#lambdasladharma <- simulateResiduals(lambdasla)
#plot(lambdasladharma)
#summary(lambdasla)

#Simpler model with PC1, PC2 and neighbours
lambdasla <- lmer(log_lambda_p1 ~ std_PC1 + std_PC2 + neighbours01 + std_log_SLA + 
                    std_PC1:std_log_SLA + std_PC2:std_log_SLA + neighbours01:std_log_SLA + (1|Site/Plot) + (neighbours01 + std_PC1 + std_PC2|Species), popdata)
lambdaslasimpledharma <- simulateResiduals(lambdasla)
plot(lambdaslasimpledharma)
#Okay, KS test significant
summary(lambdasla)
#No significant interactions

##LDMC - converges!
lambdaldmc <- lmer(log_lambda_p1 ~ neighbours01 + Treatment + std_PC1 + std_PC2 + std_LDMC + 
                     neighbours01:std_LDMC + Treatment:std_LDMC + std_PC1:std_LDMC + std_PC2:std_LDMC +
                     (1|Site/Plot) + (neighbours01 + Treatment + std_PC1 + std_PC2|Species), popdata)
lambdaldmcdharma <- simulateResiduals(lambdaldmc)
plot(lambdaldmcdharma)
#not great residuals, KS test significant
summary(lambdaldmc)
#No significant interactions

#D13C - converges!
lambdaD13C <- lmer(log_lambda_p1 ~ neighbours01 + Treatment + std_PC1 + std_PC2 + std_log_D13C + 
                         neighbours01:std_log_D13C + Treatment:std_log_D13C + std_PC1:std_log_D13C + std_PC2:std_log_D13C +
                         (1|Site/Plot) + (neighbours01 + Treatment + std_PC1 + std_PC2|Species), popdata)
lambdaD13Cdharma <- simulateResiduals(lambdaD13C)
plot(lambdaD13Cdharma)
#not great residuals, KS test significant
summary(lambdaD13C)
#PC1 and TreatmentWet significant
r.squaredLR(lambdaD13C)

#Plotting PC1 simply
species_slopes<-coef(lambdaD13C)$Species[,4] + coef(lambdaD13C)$Species[,9]*speciesD13C$std_log_D13C
plot(species_slopes ~ speciesD13C$std_log_D13C)
curve(cbind(1,x)%*%fixef(lambdaD13C)[c(4,9)], add=T)
abline(h=0, lty=3)

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


### John's code to calculate variance components in models ####
#March 2022
## Have not updated these summaries since updating data_preparation file

specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
##By species for germination
for (i in 1:length(specieslist)){
  print(specieslist[i])
  vca <- glmer(cbind(total_germ, total_no_germ) ~ 1 + (1|Site/Plot), 
            family = binomial, data = filter(vitaldata, Species == specieslist[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(vca))
  print(vca_func(vca))
}
#POLE 0
#pole had extremely low germination
polevarmod <- glmer(cbind(total_germ, total_no_germ) ~ 1 + (1|Site/Plot), 
                    family = binomial, data = poledata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
polevarmoddharma <- simulateResiduals(polevarmod)
plot(polevarmoddharma)
summary(polevarmod)
vca_func(polevarmod)

#ARCA 90% of variation in germination in among plots within blocks, 10% among blocks
#HYGL 52% among plots, 48% among blocks
#LARO 34% among plots, 66% among blocks
#PEAI 40% among plots, 60% among blocks
#PLDE 43% among plots, 57% among blocks
#POLE 65% among plots, 35% among blocks
#TRCY 69% among plots, 31% among blocks
#TROR 37% among plots, 63% among blocks
#VERO 46% among plots, 54% among blocks
# One failed to converge

##By species for survival
arcatestmod1 <- glmer(surv_to_produce_seeds ~ 1 + (1|Site/Plot), family = binomial, arcadata)
summary(arcatestmod1)
vca_func(arcatestmod1)
#Tells me that 55% of the variation in survival is among plots within blocks. 
# 45% is among blocks.
for (i in 1:length(specieslist)){
  print(specieslist[i])
  vca <- glmer(surv_to_produce_seeds ~ 1 + (1|Site/Plot), family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(vca))
  print(vca_func(vca))
}
#Something very wrong with HYGL, LARO, POLE and TRCY
#ARCA, PEAI 55% among plots, 45% among blocks
#PLDE 11% among plots, 89% among blocks
#TROR 21% among plots, 79% among blocks
#VERO 41% among plots, 59% among blocks

## By species for seed production
for (i in 1:length(specieslist)){
  print(specieslist[i])
  vca <- glmer.nb(seeds_percent ~ 1 + (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(vca))
  print(vca_func(vca))
}
#vca_func doesn't work for glmmTMB but it does with glmer.nb! :o
#Something very wrong with ARCA, HYGL, PEAI, PLDE, POLE, TROR
#LARO 63% variation in seed production among plots, 37% among blocks
#TRCY 7% among plots, 93% among blocks
#VERO 51% among plots, 49% among blocks
#Convergence issues

##### 20/04/22 meeting with John, quadratic responses to PC1/PC2/NA ####
#Plots of survival probability in response to neighbour abundance by watering treatment
# with a quadratic term and fitting lines and CIs from simple glm binomial model
ggplot(vitaldata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds, colour = Treatment))+
  geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ poly(x, 2))+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)

#Survival ~ NA by watering treatment, linear fits only
ggplot(vitaldata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds, colour = Treatment))+
  geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"))+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)

## Survival ~ PC1 by watering treatment, quadratic term
ggplot(vitaldata, aes(x = std_PC1, y = surv_to_produce_seeds, colour = Treatment))+
  geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ poly(x, 2))+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)
# Survival ~ PC1 by watering treatment, linear term
ggplot(vitaldata, aes(x = std_PC1, y = surv_to_produce_seeds, colour = Treatment))+
  geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"))+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)
## Survival ~ PC2 by watering treatment, quadratic
ggplot(vitaldata, aes(x = std_PC2, y = surv_to_produce_seeds, colour = Treatment))+
  geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ poly(x, 2))+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)
# Survival ~ PC2 by watering treatment, linear
ggplot(vitaldata, aes(x = std_PC2, y = surv_to_produce_seeds, colour = Treatment))+
  geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"))+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)






