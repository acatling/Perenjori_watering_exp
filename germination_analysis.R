######## Analysing germination data
## WA Perenjori 2020 Watering Experiment
# Script updated 22/03/22

#### Loading packages and data ####
library(tidyverse)
library(sjPlot)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)
library(DHARMa)

source("data_preparation.R")
#germdata is the tibble I want to use!

# Note: germination rate is accounting for February germination as the total number germinated / total number seeds
# may want to change this later

#Did any germinate in February and then not later in the year?
test <- germdata %>% filter(February_germination > 0 & Number_germinated == 0)
#Yes, 16 of them! Wow. Maybe I remove these plots? Not sure, keeping for now
#germinationdata <- germinationdata %>% filter(!(February_germination > 0 & Number_germinated == 0))

#### Exploratory plots ####
ggplot(germdata, aes(x = Site, y = Germination_rate))+
  geom_boxplot()+
  geom_point(alpha = 0.1)+
  theme_classic()+
  ylab("Germination rate")+
  facet_wrap(~Species)
ggplot(germdata, aes(x = PC1, y = Germination_rate))+
  geom_point(alpha = 0.1)+
  geom_smooth(method="lm")+
  theme_classic()+
  ylab("Germination rate")+
  facet_wrap(~Species)
ggplot(germdata, aes(x = PC2, y = Germination_rate))+
  geom_point(alpha = 0.1)+
  geom_smooth(method="lm")+
  theme_classic()+
  ylab("Germination rate")+
  facet_wrap(~Species)
ggplot(germdata, aes(x = PC3, y = Germination_rate))+
  geom_point(alpha = 0.1)+
  geom_smooth(method="lm")+
  theme_classic()+
  ylab("Germination rate")+
  facet_wrap(~Species)

#Make a new column for Germinated (1 or 0)
#germinationdata <- germinationdata %>% 
#  mutate(Germinated01 = case_when(Germination_rate > "0" ~ "1",
#                      Germination_rate == "0" ~ "0"))

#Trying to add an individual-level random effect to fix overdispersion
#Creating a column for individual row ID
#Omg this fixed all of the qqplot residuals!!!
germdata <- germdata %>% unite("id", Site:Column, remove = "false")

####Creating a dataset and model for each species ####
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  nam <- paste0("germ", specieslist[i])
  assign(nam, germdata %>% filter(Species == specieslist[i]))
}

#Modelling using cbind with negative binomial model
for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot) + (1|id), 
           family = binomial, data = filter(germdata, Species == specieslist[i]))))
}
## Naming models so that I can plot coefs from them
for (i in 1:length(specieslist)){
  nam <- paste0("germmod", specieslist[i])
  assign(nam, glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot) + (1|id), 
                    family = binomial, data = filter(germdata, Species == specieslist[i])))
}

#Note that this works
#Trials = cbind(germdata$total_germ, germdata$total_no_germ)
#test <- glmer(Trials ~ std_PC1 + std_PC2 + std_PC3 + (1|Site/Plot), data = germdata, family = binomial)
# Not sure if it's appropriate?!

###Checking residuals
#Problems with them - need to fix
arcagermdharma <- simulateResiduals(germmodARCA)
plot(arcagermdharma)
hyglgermdharma <- simulateResiduals(germmodHYGL)
plot(hyglgermdharma)
larogermdharma <- simulateResiduals(germmodLARO)
plot(larogermdharma)
summary(germmodLARO)
PEAIgermdharma <- simulateResiduals(germmodPEAI)
plot(PEAIgermdharma)
summary(germmodPEAI)
PLDEgermdharma <- simulateResiduals(germmodPLDE)
plot(PLDEgermdharma)
#Extremely unhappy
summary(germmodPLDE)
POLEgermdharma <- simulateResiduals(germmodPOLE)
plot(POLEgermdharma)
summary(germmodPOLE)
TRCYgermdharma <- simulateResiduals(germmodTRCY)
plot(TRCYgermdharma)
summary(germmodTRCY)
TRORgermdharma <- simulateResiduals(germmodTROR)
plot(TRORgermdharma)
summary(germmodTROR)
VEROgermdharma <- simulateResiduals(germmodVERO)
plot(VEROgermdharma)
summary(germmodVERO)

## Trying to figure out what's wrong
testDispersion(VEROgermdharma)
#Overdispersed - everything except vero, tror
par(mfrow = c(1,2))
plot(germPEAI$PC1, germPEAI$percent_germ, xlab = "Envrionmental Predictor", ylab = "Response")
hist(germPEAI$percent_germ, xlab = "Response", main = "")
#Some of them are zero-inflated?
testZeroInflation(PLDEgermdharma)

############# Plotting germination panels/figures ######
#Something wrong with this, 22/03/22
#ARCA
#x_to_plot<-seq.func(germARCA$std_PC1)
#with(germARCA, plot(percent_germ ~ std_PC1))
#arcapreddata <- with(germmodARCA, data.frame(1, x_to_plot, 0, 0))
#arcapred <- glmm.predict(mod = germmodARCA, newdat = arcapreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
#plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)

species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
dev.off()
pdf("Output/Figures/germ_PC1.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
#Margins: bottom, left, top, right
#A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region
par(pty="s")
species.list<-list(germARCA, germHYGL, germLARO, germPEAI, germPLDE, germPOLE, germTRCY, germTROR, germVERO)
for(i in 1:length(species.list)){
  plotted.data<-as.data.frame(species.list[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC1, pch=19, col="grey60", ylab="Germination rate", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "PC1", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(cbind(total_germ, total_no_germ)~std_PC1+std_PC2+std_PC3+(1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC1) 
  preddata <- with(model, data.frame(1, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
   }
dev.off()

#Simpler version for ESA presentation:
dev.off()
pdf("Output/Figures/germ_PC1_ESA.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(2,2,6,0.5))
#Margins: bottom, left, top, right
par(pty="s")
species.list<-list(germARCA, germHYGL, germLARO, germPEAI, germPLDE, germPOLE, germTRCY, germTROR, germVERO)
for(i in 1:length(species.list)){
  plotted.data<-as.data.frame(species.list[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC1, pch=19, col="grey60", ylab="", xlab="", cex = 2, cex.lab=3, cex.axis=3,tck=-0.02, axes = 'FALSE', frame.plot = TRUE,ylim=c(0,1.0))
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=4)
  model<-glmer(cbind(total_germ, total_no_germ)~std_PC1+std_PC2+std_PC3+(1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC1)
  preddata <- with(model, data.frame(1, x_to_plot, 0, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, 
               env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 4, line.type = 1)
  Axis(side=1, labels=FALSE)
  Axis(side=2, labels=FALSE)
}
dev.off()

## And for PC2
dev.off()
pdf("Output/Figures/germ_PC2.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
#Margins: bottom, left, top, right
#A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region
par(pty="s")
species.list<-list(germARCA, germHYGL, germLARO, germPEAI, germPLDE, germPOLE, germTRCY, germTROR, germVERO)
for(i in 1:length(species.list)){
  plotted.data<-as.data.frame(species.list[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC2, pch=19, col="grey60", ylab="Germination rate", xlab="", cex.lab=3, cex.axis=3,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=2)
  title(xlab = "PC2", cex.lab=3)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=3)
  model<-glmer(cbind(total_germ, total_no_germ)~std_PC1+std_PC2+std_PC3+(1|Site/Plot), family = binomial, plotted.data)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  preddata <- with(model, data.frame(1, 0, x_to_plot, 0))
  plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
  plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
}
dev.off()

########### Haven't updated anything below here ##########

### Modelling germination rate by species ###
# for (i in 1:length(specieslist)){
#   print(specieslist[i])
#   print(
#     summary(
#     lmer(log_germ ~ std_PC1 + std_PC2 + std_PC3 + (1|Site/Plot), 
#           data = filter(germdata, Species == specieslist[i]))))
# }

tab_model(germmodelARCA, germmodelHYGL, germmodelLARO, germmodelPLDE, germmodelPOLE, germmodelPEAI, germmodelTRCY, germmodelTROR, germmodelVERO, transform = NULL)

##Plotting coefficients
allmodels <- list()
allmodels[[1]] <- germmodelARCA
allmodels[[2]] <- germmodelHYGL
allmodels[[3]] <- germmodelLARO
allmodels[[4]] <- germmodelPLDE
allmodels[[5]] <- germmodelPOLE
allmodels[[6]] <- germmodelPEAI
allmodels[[7]] <- germmodelTRCY
allmodels[[8]] <- germmodelTROR
allmodels[[9]] <- germmodelVERO

plot_models(allmodels, transform = NULL, vline.color = "grey")+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
#Can't make x axis limits any smaller using sJplot (can use a work around with plot_model in ggplot if I want to)


### OLD WORK BELOW ###
###### Need to calculate germination rates per subplot (number germinated/total number of seeds)
# Model by germination rate
# 
# for (i in 1:length(specieslist)){
#   print(specieslist[i])
#   print(
#     summary(
#       lmer(Germination_rate ~ log(cc_percentage+1) + (1|Site/Plot),
#            data = filter(germinationdata, Species == specieslist[i]))))
# }
# #All have a significant negative coefficient for canopy closure except ARCA
# 
# for (i in 1:length(specieslist)){
#   nam <- paste0("germmodelc", specieslist[i])
#   assign(nam, lmer(Germination_rate ~ log(cc_percentage+1) + (1|Site/Plot),  
#                     data = filter(germinationdata, Species == specieslist[i])))
# }
# ####Calculating germination rates
# #Import Sun vs. Shade info
# treatmentdata <- read_csv("Data/treatments_meta.csv")
# germinationdata <- merge(germinationdata, treatments)
# 
# germinationcounts <- germinationdata %>% group_by(Species, Cover, Germinated01) %>% tally()
# #I already have germination rates calculated... so I just need to average these?
# germinationsp <- germinationdata %>% group_by(Species) %>% summarise(mean_germ_rates = mean(Germination_rate),
#                                                             sd_germ_rates = sd(Germination_rate))
# germinationcoversp <- germinationdata %>% group_by(Species, Cover) %>% summarise(mean_germ_rates = mean(Germination_rate),
#                                                                           sd_germ_rates = sd(Germination_rate))
# #Comparing these species-level values to another way of measuring it
# germcount <- germinationdata %>% group_by(Species, Germinated01) %>% tally()
# germsp <- germcount %>%
#   group_by(Species) %>%
#   mutate(number_germinated = sum(n[Germinated01 == "1"]),
#          total_number = sum(n[Germinated01 == "1"])+sum(n[Germinated01 == "0"]),
#          germination_rate = number_germinated/total_number) %>%
#   filter(row_number() == 1)
# #Selecting only relevant columns
# germsp <- germsp %>% select(Species, germination_rate)
# ################ THIS GIVES VERYYYY DIFFERENT VALUES, TRUST THIS MORE.
# #And split between sun/shade
# germshadecount <- germinationdata %>% group_by(Species, Cover, Germinated01) %>% tally()
# germshadesp <- germshadecount %>%
#   group_by(Species, Cover) %>%
#   mutate(number_germinated = sum(n[Germinated01 == "1"]),
#          total_number = sum(n[Germinated01 == "1"])+sum(n[Germinated01 == "0"]),
#          germination_rate = number_germinated/total_number) %>%
#   filter(row_number() == 1)
# #Selecting only relevant columns
# germshadesp <- germshadesp %>% select(Species, Cover, germination_rate)
# 
# 
# ggplot(germinationdata, aes(x=Species, y=Germination_rate)) + 
#   geom_jitter(colour="forestgreen", size=1.6, alpha = 0.2,width = 0.2)+
#   geom_point(stat="summary", size=2) +
#   geom_errorbar(stat="summary", width=0,size=0.8)+
#   theme_classic()
# 
# #Defaults to mean and se
# #ggplot(germinationdata, aes(x=Cover, y=Germination_rate)) + 
# #  geom_jitter(colour="forestgreen", size=1.6, alpha = 0.2,width = 0.2)+
# #  geom_point(stat="summary", size=1.5) +
# #  geom_errorbar(stat="summary", width=0,size=0.8)+
# #  theme_classic()+
# #  facet_wrap(~Species)
# 
# ggplot(germinationdata, aes(x = Cover, y = Germination_rate))+
#   geom_boxplot()+
#   geom_jitter(colour = "steelblue3", alpha = 0.2, width = 0.2)+
#   theme_classic()+
#   facet_wrap(~Species)
# 
# #Series of t-tests
# 
# for (i in 1:length(specieslist)){
#   print(specieslist[i])
#   print(
#     t.test(Germination_rate ~ Cover, data = filter(germinationdata, Species == specieslist[i])))
# }
# #All significantly different except ARCA. All germinated more in the sun (except ARCA).
# 
# ggplot(germinationdata, aes(x = log(cc_percentage+1), y = Germination_rate))+
#   geom_point(alpha = 0.4)+
#   theme_classic()+
#   geom_smooth(method="lm")+
#   facet_wrap(~Species)
# 
# r.squaredGLMM(germmodelcARCA)
# r.squaredGLMM(germmodelcHYGL)
# r.squaredGLMM(germmodelcLARO)
# r.squaredGLMM(germmodelcPEAI)
# r.squaredGLMM(germmodelcPLDE)
# r.squaredGLMM(germmodelcPOLE)
# r.squaredGLMM(germmodelcTRCY)
# r.squaredGLMM(germmodelcTROR)
# r.squaredGLMM(germmodelcVERO)
####

## Trying Isaac's plots ###
source("Isaac_functions.R")

###create a list of the dataframes which only have C data
species.c.list<-list(germmodARCA, germmodHYGL, germmodLARO, germmodPEAI, germmodPLDE, germmodPOLE, germmodTRCY, germmodTROR, germmodVERO)

###create a list of names for figure headings
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

dev.off()
pdf("Output/Figures/testing2610.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(2,6,2,1))
#pty
#A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region
par(pty="s")

for(i in 1:length(species.c.list)){
  plotted.data<-as.data.frame(species.c.list[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC1, pch=19, col="grey60", ylab="Germination rate", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "PC1", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer(percent_germ~std_PC1+std_PC2+std_PC3+(1|Site/Plot), plotted.data)
  if(summary(model)$coefficients[2,3]<0.05){
    x.for.plot<-plotted.data$std_PC1
    y.for.plot<-plotted.data$percent_germ
    pred.data<-data.frame(intercept=rep(1, 100), x.for.plot=seq.func(x.for.plot))
    pred<-lmer.predict(model, newdat=pred.data, se.mult=1.96, binom=T, poisson=F)
    plot.CI.func(x.for.plot=seq.func(x.for.plot), pred=pred$y, upper=pred$phi, lower=pred$plo, env.colour="black", env.trans=40, line.colour="black", line.type=1, line.weight=3)
  }
}
dev.off()
#WHY NO LINES COMING UP???? ###### Not working

#Decide what should and shouldn't be logged
par(mfrow=c(3,3))
#ARCA
x_to_plot<-seq.func(germARCA$std_PC1)
with(germARCA, plot(log_germ ~ std_PC1))
arcapreddata <- with(germmodARCA, data.frame(1, x_to_plot, 0, 0))
arcapred <- glmm.predict(mod = germmodARCA, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)
#HYGL
x_to_plot<-seq.func(germHYGL$std_PC1)
with(germHYGL, plot(log_germ ~ std_PC1))
arcapreddata <- with(germmodHYGL, data.frame(1, x_to_plot, 0, 0))
arcapred <- glmm.predict(mod = germmodHYGL, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)
#LARO
x_to_plot<-seq.func(germLARO$std_PC1)
with(germLARO, plot(log_germ ~ std_PC1))
arcapreddata <- with(germmodLARO, data.frame(1, x_to_plot, 0, 0))
arcapred <- glmm.predict(mod = germmodLARO, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)
#PEAI
x_to_plot<-seq.func(germPEAI$std_PC1)
with(germPEAI, plot(log_germ ~ std_PC1))
arcapreddata <- with(germmodPEAI, data.frame(1, x_to_plot, 0, 0))
arcapred <- glmm.predict(mod = germmodPEAI, newdat = arcapreddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)


with(trordata, plot(log_germ ~ std_PC1))
trorpreddata <- with(trordata, data.frame(1, x_to_plot, 0, 0, 0, 0, 0))
trorpred <- glmm.predict(mod = trorsurvmod2, newdat = trorpreddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = arcapred$y, upper = arcapred$upper, lower = arcapred$lower, env.colour = "blue", env.trans = 50, line.colour = "blue", line.weight = 2, line.type = 1)
