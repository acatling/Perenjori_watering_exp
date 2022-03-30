### Answering my specific questions
# WA Perenjori 2020 Watering Experiment
#Updated 22/03/22

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
#dataall has everything (germination, survival, seed production and neighbour info) combined

## Note that germination analysis is separate, in germination_analysis.R

#### Is neighbour abundance correlated with abiotic environmental factors? ####
glimpse(dataall)
#Worried about neighbour 0 count being assigned to things that weren't surveyed
# did I survey the plots that never germinated? hm...
#dataall is only the subplots that germinated
#dataall has rows for subplots that germinated but died before surveying communities
#1728 subplots sown.
#Just need to check why some plots weren't surveyed - I thought I surveyed everything that germinated
# 53 rows
#test <- surveydataraw %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
#  filter(row_number() ==1)
# everything that died, regardless of when, I went back and surveyed
# did not survey subplots that never germinated
hist(dataall$Total_abundance)
hist(log(dataall$Total_abundance+1))
#dataset with neighbours only
datanonly <- dataall %>% filter(Total_abundance > 0)
hist(datanonly$Total_abundance)
hist(log(datanonly$Total_abundance+1))

#Trying to plot them side by side as pdf
# All of the below variables are standardised - I don't think that matters
#PC1
pc1plot1 <- ggplot(dataall, aes(x = std_PC1, y = log(Total_abundance+1)))+
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
#PC2
pc2plot1 <- ggplot(dataall, aes(x = std_PC2, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()
pc2plot2 <- ggplot(datanonly, aes(x = std_PC2, y = log(Total_abundance+1)))+
  geom_jitter(alpha = 0.1, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("")+
  theme_classic()
#Water availability
waterplot1 <- ggplot(dataall, aes(x = Treatment, y = log(Total_abundance+1)))+
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
modelabundwater <- aov(log(Total_abundance+1) ~ Treatment, datanonly)
summary(modelabundwater)
TukeyHSD(modelabundwater)

## More plots
#Rename Control watering treatment to Ambient
datanonly <- datanonly %>% mutate(Treatment = recode(Treatment, Control = 'Ambient'))
#Reorder watering treatments to  Dry, Ambient, Wet
datanonly$Treatment <- factor(datanonly$Treatment, level = c("Dry", "Ambient", "Wet"))

#Plotting PC1 vs abund for watering treatments
waterpc1plot <- ggplot(datanonly, aes(x = std_PC1, y = logp1_totalabund))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  geom_smooth()+
  ylab("log(neighbour abundance + 1)")+
  theme_classic()+
  my_theme+
  facet_wrap(~Treatment)
#Plotting PC2 vs abund for watering treatments
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

### Are PC1 and PC2 correlated? ####
#Only want one data point per plot
plotvar <- dataall %>% group_by(Site, Plot) %>% filter(row_number() == 1)
ggplot(plotvar, aes(x = std_PC1, y = std_PC2))+
  geom_jitter(width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  theme_classic()
modpc1pc2 <- lm(std_PC2 ~ std_PC1, plotvar)
summary(modpc1pc2)
### How is neighbour abundance distributed across species across abiotic env? ####
sppc1abund <- ggplot(dataall, aes(x = std_PC1, y = logp1_totalabund))+
  geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
  ylab("log(Neighbour abundance + 1)")+
  theme_classic()+
  facet_wrap(~Species)
sppc2abund <- ggplot(dataall, aes(x = std_PC2, y = logp1_totalabund))+
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

#### Question 1 relative importance A. and B. factors####
##What is the relative importance of abiotic and biotic factors for survival and fecundity?
#Additive model: response ~ total_abundance + water + PC1 + PC2 + PC3 + RE
### Q1 - Survival to produce seeds ####
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survival <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                  (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(survival))
}
# ###Putting models into a list for coef plotting
# survivalmodels <- list()
# for (i in 1:length(specieslist)){
#   survivalmodels[[i]] <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_PC3 + Dodder01 +
#                                  (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]))
# }
# #Making coefficient plot for survival as a function of environment
# plot_models(survivalmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
#             dot.size = 2, line.size = 1)+
#   ylab("Estimate")+
#   theme_classic()+
#  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
### Assigning model names by species to plot
#Running it with this bobqya optimiser so that vero converges
for (i in 1:length(specieslist)){
  nam <- paste0("survmod", specieslist[i])
  assign(nam, glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                      (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))))
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
summary(survmodPLDE)
POLEsurvdharma <- simulateResiduals(survmodPOLE)
plot(POLEsurvdharma)
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
survmodels[[6]] <- survmodPOLE
survmodels[[7]] <- survmodTRCY
survmodels[[8]] <- survmodTROR
survmodels[[9]] <- survmodVERO
tab_model(survmodels, transform = NULL)
plot_models(survmodels, transform = NULL, vline.color = "grey", legend.title = "Species",
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

# #Plots of significant with neighbour abundance (ARCA and TROR)
# par(mfrow=c(2,1))
# x_to_plot<-seq.func(dataall$std_logp1_totalabund)
# with(arcadata, plot(ProducedSeeds ~ std_logp1_totalabund))
# arcapreddata <- with(arcadata, data.frame(1, x_to_plot, 0, 0, 0, 0, 0, 0))
# arcapred <- glmm.predict(mod = survmodARCA, newdat = arcapreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
# plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)
# with(trordata, plot(ProducedSeeds ~ std_logp1_totalabund))
# trorpreddata <- with(trordata, data.frame(1, x_to_plot, 0, 0, 0, 0, 0))
# trorpred <- glmm.predict(mod = trorsurvmod2, newdat = trorpreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
# plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

### Plotting with all species together
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
dev.off()
pdf("Output/Figures/survival_PC1.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
species.list.s<-list(arcadata, hygldata, larodata, peaidata, pldedata, poledata, trcydata, trordata, verodata)
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$ProducedSeeds~plotted.data$std_PC1, pch=19, col="grey60", ylab="Survival rate", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "PC1", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                 (1|Site/Plot), family = binomial, data = plotted.data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  x_to_plot<-seq.func(plotted.data$std_PC1)
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
  plot(plotted.data$ProducedSeeds~plotted.data$std_PC1, pch=19, col="grey60", ylab="", xlab="", cex = 3, cex.lab=3, cex.axis=3,tck=-0.02, alpha = 0.4, axes = 'FALSE', frame.plot = TRUE)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=4)
  model<-glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
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

### Surv - Quantify variance in vital rates explained by abiotic only vs biotic only ####
## Another way of this answering this question, instead of just looking at # of significant interactions
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
#Abiotic only
#Need optimiser for convergence
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survivalabiotic <- glmer(ProducedSeeds ~ Treatment + std_PC1 + std_PC2 +
                      (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
  print(summary(survivalabiotic))
}


##################################################################
### Q1 - Viable seed production #####
#Using No_viable_seeds_grouped for individual species models
#and seeds_percent for trait models
#trcy doesn't converge with No_viable_seeds_grouped
#Same results both ways
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  fecundity <- glmmTMB(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                      (1|Site/Plot), family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(fecundity))
}
#Assigning model names by species
for (i in 1:length(specieslist)){
  nam <- paste0("seedmod", specieslist[i])
  assign(nam, glmmTMB(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
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
# POLE unhappy. Need to try simpler model 
#or don't run it for this vital rate. Just not enough data??
# Only one occurrence of Dodder. Not a lot of neighbour data, but nor does ARCA
POLEseeddharma <- simulateResiduals(seedmodPOLE)
plot(POLEseeddharma)
summary(seedmodPOLE)
testDispersion(POLEseeddharma)
#Underdispersed! 
testZeroInflation(POLEseeddharma)
#Not zero inflated
#Trying a Conway-Maxwell-Poisson distribution (good for underdispersion)
#Needs to be within glmmTMB to have random effects
library(mpcmp)
poleseedmod3 <- glmmTMB(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                          (1|Site/Plot), family = compois, data = seedpole)
summary(poleseedmod3)
POLEseeddharma3 <- simulateResiduals(poleseedmod3)
plot(POLEseeddharma3)
#MUCH better residuals
#Testing it for plde - very similar result!! Same residuals
#pldeseedmod3 <- glmmTMB(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
#                          (1|Site/Plot), family = compois, data = seedplde)
#summary(pldeseedmod3)
summary(seedmodPLDE)
poleseedmod2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                          (1|Site/Plot), family = nbinom2, data = seedpole)
summary(poleseedmod2)
with(seedpole, plot(seeds_percent ~ std_logp1_totalabund))
POLEseeddharma2 <- simulateResiduals(poleseedmod2)
plot(POLEseeddharma2)
#
TRCYseeddharma <- simulateResiduals(seedmodTRCY)
plot(TRCYseeddharma)
summary(seedmodTRCY)
#Doesn't converge with no_viable_seeds_grouped
#Seeing if it converges with an optimiser - no.
#control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
trcyseedmod2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 +
                          (1|Site/Plot), family = nbinom2, data = seedtrcy)
summary(trcyseedmod2)
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
species.list.s<-list(seedarca, seedhygl, seedlaro, seedpeai, seedplde, seedpole, seedtrcy, seedtror, seedvero)
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
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
species.list.s<-list(seedarca, seedhygl, seedlaro, seedpeai, seedplde, seedpole, seedtrcy, seedtror, seedvero)
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
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
## Do species-level leaf traits explain responses to environment?
##Starting with simple example. Does SLA modulate the respones of survival to PC1?
### Survival models
survpc1sla <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_SLA + 
                      std_PC1:std_log_SLA + std_logp1_totalabund:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
                      (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, dataall, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(survpc1sla)
survpc1ldmc <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_LDMC + 
                       std_PC1:std_LDMC + std_logp1_totalabund:std_LDMC + std_PC2:std_LDMC + Treatment:std_LDMC +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, dataall, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(survpc1ldmc)
survpc1D13C <- glmer(ProducedSeeds ~ std_logp1_totalabund + Treatment + Dodder01 + std_PC1 + std_PC2 + std_log_D13C + 
                       std_PC1:std_log_D13C + std_logp1_totalabund:std_log_D13C + std_PC2:std_log_D13C + Treatment:std_log_D13C +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = binomial, dataall, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(survpc1D13C)

with(dataall, plot(ProducedSeeds ~ std_log_D13C))
x_to_plot<-seq.func(dataall$std_log_D13C)
wuepreddata <- with(dataall, data.frame(1, 0, x_to_plot, 0*x_to_plot))
wuepred <- glmm.predict(mod = survpc1wue, newdat = wuepreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = wuepred$y, upper = wuepred$upper, lower = wuepred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

######### Fecundity traits ####
### Fecundity models

#Will only converge with optimiser
seedldmc <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + std_LDMC + 
                       std_PC1:std_LDMC + std_logp1_totalabund:std_LDMC + std_PC2:std_LDMC + Treatment:std_LDMC +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, 
                     seedmodeldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(seedldmc)
seedldmcdharma <- simulateResiduals(seedldmc)
plot(seedldmcdharma)

seedD13C <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + Dodder01 + std_log_D13C + 
                       std_PC1:std_log_D13C + std_logp1_totalabund:std_log_D13C + std_PC2:std_log_D13C + Treatment:std_log_D13C +
                       (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, seedmodeldata,
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(seedD13C)
seedD13Cdharma <- simulateResiduals(seedD13Cdharma)
plot(seedD13Cdharma)

######## SLA analysis and plotting below! #########
## Saving the model to load later
save(seedsla, file = "Output/Models/seedslamodel.RData")
# Loading in the model from saved output
load("Output/Models/seedslamodel.RData")

#Trying optimiser to aid convergence
#control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))
#Trying without Treatment and Dodder. Didn't work but did converge with glmer.nb instead of glmmTMB!
seedsla <- glmer.nb(seeds_percent ~ std_logp1_totalabund + Treatment + std_PC1 + std_PC2 + std_log_SLA + Dodder01 +
                        std_PC1:std_log_SLA + std_logp1_totalabund:std_log_SLA + std_PC2:std_log_SLA + Treatment:std_log_SLA +
                        (1|Site/Plot) + (std_PC1 + std_PC2 + std_logp1_totalabund + Treatment|Species), family = nbinom2, seedmodeldata)
summary(seedsla)
seedsladharma <- simulateResiduals(seedsla)
plot(seedsladharma)

#### Need to check if below here was updated after lab retreat, not sure! ###

### Plotting this with John's Catford script and my amendments
### This works!!
#Creating a dataframe with SLA data - need the as.data.frame part!
speciessla <- seedmodeldata %>% select("Species", "std_log_SLA") %>% filter(row_number() == 1)
speciessla <- as.data.frame(speciessla)
#Matching SLA values to Species in the model
traits_for_model<-speciessla[match(rownames(coef(seedsla)$Species), speciessla$Species),c(1:2)]
#John: this bit grabs the random slopes for each species as well as their associated SEs
#NOT SURE THESE ARE THE RIGHT VALUES?** [2,2?]
varfix <- vcov(seedsla)[2,2]
re <- ranef(seedsla,condVar=TRUE)
varcm <- attr(re$Species,"postVar")[2,2,]
vartot <- varfix+varcm
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
#Calculating the slope values for the response of fecundity to PC1 as modulated by SLA
# slopes = the slope responses of PC2 for each species + the interaction of PC2:SLA*each species' SLA value
#Unique value for PC2 slopes per species (what's left over), one value for PC2:SLA for all species, unique values of SLA per species
#This changes the position of each species: PC2:SLA - which is one it's one value
# This is all because the fixed effects have explained some of the response already
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(seedsla)$Species[,5] + coef(seedsla)$Species[,9]*traits_for_model$std_log_SLA,
                                         slope_PC2 = coef(seedsla)$Species[,6] + coef(seedsla)$Species[,11]*traits_for_model$std_log_SLA)
### Plotting response of fecundity to PC1 as a function of SLA
#c(5,9) is pulling out std_PC1 and std_PC1*std_log_SLA
#WHY DO WE CALL DATA FROM OUR BIG DATASET IN THIS FIRST LINE EXACTLY? To calculate the CIs or error around the points?
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(5,9)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(5,9), c(5,9)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ std_log_SLA, ylab= "Slope of fecundity - PC1 relationship", xlab="log(SLA) (standardised)", ylim=c(-1,.3), pch = 19, col="grey70"))
with(slopedata, arrows(std_log_SLA, slope_PC1+slope.SEs, std_log_SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
abline(h=0, lty=3, lwd=2)

### Plotting response of fecundity to PC2 as a function of SLA
#c(6,11) is pulling out std_PC2 and std_PC2*std_log_SLA
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(6,11)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(6,11), c(6,11)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC2 ~ std_log_SLA, ylab= "Slope of fecundity - PC2 relationship", xlab="log(SLA) (standardised)", ylim=c(-0.6, 0.5), pch = 19, col="grey70"))
with(slopedata, arrows(std_log_SLA, slope_PC2+slope.SEs, std_log_SLA, slope_PC2-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
abline(h=0, lty=3, lwd=2)

#### This plots just the regression line, without CIs - for PC1-SLA
#curve(cbind(1,x)%*%fixef(seedsla)[c(5,9)], add=T, lwd=3)
### Not necessary, but can label the rownames the species names
#rownames(speciessla) <- speciessla$Species
### pch = 19 for plotting makes the data points filled in

##THis pulls out the fixed effect values in a table
#test <- summary(seedsla)$coef

#As a bit of a check, run a model without interaction and these should be the slope values we want ultimately
# after adding the things, this:
## slopes = the slope responses of PC2 for each species + the interaction of PC2:SLA*each species' SLA value
#Very similar - not exactly the same!

## NEED TO UPDATE MY MODIFIERS OF SLOPES FROM MY BIG MODEL 
# everything that is interacting with SLA!

########### Haven't updated anything below here ###########
### Various working out below (what works is above):
#Creating a dataframe with SLA data, my code;
speciessla <- seedmodeldata %>% select("Species", "std_log_SLA") %>% filter(row_number() == 1)
speciessla <- as.data.frame(speciessla)
rownames(speciessla) <- speciessla$Species
#speciessla <- speciessla %>% select("std_log_SLA")

### John's working 07/11/2021
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
### End John's

##John's code:
traits_for_model<-speciessla[match(rownames(coef(seedsla)$Species), speciessla$Species),c(1:2)]
#GETTING STUCK IN BELOW LINE. Not sure I'm calling the right coef columns, tried many and none work
traits_for_model$slopes<- rowSums(cbind(1, traits_for_model) * as.data.frame(coef(seedsla)$Species[,c(7, 9, 11)]))
#Below is my interpretation of what that line does
traits_for_model <- traits_for_model %>% mutate(slopes = coef(seedsla)$Species[,c(7, 9, 11)])
#Back to John's code
traits_for_model$slope.SEs<-slope.SEs<-sqrt(vartot)
traits_for_model$Species<-rownames(coef(seedsla)$Species)

## Updating this different way 08/11
slopedata <- traits_for_model %>% mutate(slope_PC1 = coef(seedsla)$Species[,5] + coef(seedsla)$Species[,9]*traits_for_model$std_log_SLA,
                                         slope_PC2 = coef(seedsla)$Species[,6] + coef(seedsla)$Species[,11]*traits_for_model$std_log_SLA)

#Trying it a different way
#Extracting slopes for each species
#Note that the slopes are called std_PC1 and std_PC2 (confusing oops), renaming
beta <- coef(seedsla)$Species
signifbeta <- beta %>% select("std_PC1", "std_PC2")
signifbeta <- signifbeta %>% rownames_to_column("Species")
signifbeta <- signifbeta %>% rename(slope_PC1 = std_PC1,
                              slope_PC2 = std_PC2)
#Merging datasets
slopespsla <- merge(speciessla, signifbeta)
#Creating a column for SE of slopes
slopespsla <- slopespsla %>% mutate(slope.SEs = sqrt(vartot))

### SLA - PC1 slope response
#Revelation! c(5,9) is pulling out std_PC1 and std_PC1*std_log_SLA
#But shouldn't it be std_log_SLA and std_log_SLA*std_PC1...? Nah, not what John did
#Not sure this is working yet.
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(5,9)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(5,9), c(5,9)], as.matrix(data.for.sla.slope.regression))))
with(slopedata, plot(slope_PC1 ~ std_log_SLA, ylab= "Slope of fecundity - PC1 relationship", xlab="log(SLA) (standardised)", ylim=c(-1,.3)))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
with(slopedata, points(std_log_SLA, slope_PC1, col="grey70", pch=19))
with(slopedata, arrows(std_log_SLA, slope_PC1+slope.SEs, std_log_SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
#curve(cbind(1,x)%*%fixef(seedsla)[c(9,11)], add=T, lwd=3)
abline(h=0, lty=3, lwd=2)

####Testing a few things with Isis to see if the code is working
# with(slopespsla, plot(slope_PC1 ~ std_log_SLA))
# test <- slopespsla %>% select(Species, slope_PC1, slope_PC2)
# test2 <- merge(seedmodeldata, test)
# x_to_plot<-seq.func(test2$std_log_SLA)
# with(test2, plot(slope_PC1 ~ test2$std_log_SLA))
# arcapreddata <- with(test2, data.frame(1, 0, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot, 0*x_to_plot, 0*x_to_plot, 0*x_to_plot, 0*x_to_plot))
#GETTING AN ERROR
# arcapred <- glmm.predict(mod = seedsla, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
# plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

### cath ###
dev.off()
plot(slopespsla$std_log_SLA, slopespsla$slope_PC2)
arrows(slopespsla$std_log_SLA, slopespsla$slope_PC2+slopespsla$slope.SEs, slopespsla$std_log_SLA, slopespsla$slope_PC2-slopespsla$slope.SEs, code=3, angle=90, length=0)
m <- lm(slopespsla$slope_PC2~slopespsla$std_log_SLA) # not good practice 
abline(m)

set_theme(theme_classic())
plot_model(seedsla, type = "int", mdrt.values = "minmax", show.data = T, show.legend = F)

dev.off()
pdf("Output/slapc-0511.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
#### End Cath ###


#Testing another way; aw man the order of 5/9 or 9/5 matters too...
# data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
# pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(9, 5)] # why do you need to do a matrix multiplication? - cath 
# sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(9, 5), c(9, 5)], as.matrix(data.for.sla.slope.regression))))
# with(slopespsla, plot(slope_PC1 ~ std_log_SLA, ylab= "Slope of fecundity - PC1 relationship", xlab="log(SLA) (standardised)", ylim=c(-1,.3)))
# plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
#              upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
#              lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
#              env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
# with(slopespsla, points(std_log_SLA, slope_PC1, col="grey70", pch=19))
# with(slopespsla, arrows(std_log_SLA, slope_PC1+slope.SEs, std_log_SLA, slope_PC1-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
# #curve(cbind(1,x)%*%fixef(seedsla)[c(9,11)], add=T, lwd=3)
#abline(h=0, lty=3, lwd=2)

###### SLA - PC2 slope response
#6,11 and 7,11 give the same thing?
data.for.sla.slope.regression<-cbind(int=1, x=seq.func(seedmodeldata$std_log_SLA))
pred.sla.slope.regression<- data.for.sla.slope.regression%*%fixef(seedsla)[c(6, 11)]
sla.slope.regression.SEs<-sqrt(diag(as.matrix(data.for.sla.slope.regression) %*% tcrossprod(vcov(seedsla)[c(6, 11), c(6, 11)], as.matrix(data.for.sla.slope.regression))))
with(slopestest, plot(slope_PC2 ~ std_log_SLA, ylab= "Slope of fecundity - PC2 relationship", xlab="log(SLA) (standardised)", ylim=c(-0.6, 0.5)))
plot.CI.func(x.for.plot=seq.func(seedmodeldata$std_log_SLA), pred=pred.sla.slope.regression,
             upper=pred.sla.slope.regression + 1.96*sla.slope.regression.SEs,
             lower=pred.sla.slope.regression - 1.96*sla.slope.regression.SEs,
             env.colour='grey30', env.trans=40, line.colour='black', line.weight=3, line.type=1)
with(slopestest, points(std_log_SLA, slope_PC2, col="grey70", pch=19))
with(slopestest, arrows(std_log_SLA, slope_PC2+slope.SEs, std_log_SLA, slope_PC2-slope.SEs, length = 0, angle = 30, code = 2, col = "grey70", lwd = 1))
abline(h=0, lty=3, lwd=2)

#### Calculating species richness and diversity ####
#Not including dodder because we don't have an abundance
# count for it so can't calculate diversity inc dodder
#Species richness
diversitydata <- seedabundancerows %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  add_tally(Neighbour_count > 0, name = "sp_richness")
#SDI = Simpsons diversity index
diversitydata <- diversitydata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
    mutate(SDI_step = Neighbour_count/Total_abundance,
           SDI_step2 = SDI_step^2,
          SDI = 1-(sum(SDI_step2)))
#NAs are where there are no neighbours (or just dodder). Changing SDI to 0
diversitydata <- within(diversitydata, SDI[is.na(SDI)] <-  '0')
diversitydata <- diversitydata %>% filter(row_number() == 1)

### John's code to calculate variance components in models ####
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
##By species for germination
#Is having (1|id) as a random effect legitimate? Can't remember
#I have this for my other germination models
for (i in 1:length(specieslist)){
  print(specieslist[i])
  vca <- glmer(cbind(total_germ, total_no_germ) ~ 1 + (1|Site/Plot), 
            family = binomial, data = filter(germdata, Species == specieslist[i]))
  print(summary(vca))
  print(vca_func(vca))
}
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
arcatestmod1 <- glmer(ProducedSeeds ~ 1 + (1|Site/Plot), family = binomial, arcadata)
summary(arcatestmod1)
vca_func(arcatestmod1)
#Tells me that 55% of the variation in survival is among plots within blocks. 
# 45% is among blocks.
for (i in 1:length(specieslist)){
  print(specieslist[i])
  vca <- glmer(ProducedSeeds ~ 1 + (1|Site/Plot), family = binomial, data = filter(dataall, Species == specieslist[i]))
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


