##### Created 21/04/22
# This script will contain the code needed to generate figures and tables in the paper
# WA Perenjori Watering Experiment
# Alexandra Catling

#### Loading packages and data ####
#Data imported from data preparation sheet
source("data_preparation.R")
#vitaldata has all datasets combined 
# germination, survival, seed production, neighbour info, traits, quad factors
# one row per subplot with seeds sown - NAs are very important since, e.g. survival info is only on germinated subplots

#Functions file
source("R_functions/functions.R")

#Packages
library(ggplot2)
library(DHARMa)
library(glmmTMB)
library(kableExtra)
library(car)
library(sjPlot)
library(gridExtra)

#ggplot theme for ease of plotting. Use theme_classic()+ my_theme
my_theme <- theme(axis.title.x = element_text(size = 14, face = 'bold'),
                  axis.title.y = element_text(size = 14, face = 'bold'),
                  axis.text = element_text(size = 14),
                  strip.text.x = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.title = element_blank())

#### Determining whether each species needs quadratic term ####
### Continuous variables PC1, PC2, neighbour abundance

# Could do it by species, but would rather use loops:
# #ARCA - linear PC1 survival
# arcasurvquad1 <- glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, arcadata)
# summary(arcasurvquad1)

### Germination ###
## PEAI - quadratic PC1 germination
#LARO - quadratic PC1 germination
## All other species - linear PC1 germination
for (i in 1:length(specieslist)){
  print(specieslist[i])
  germquadPC1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(germquadPC1))
}
#Plotting if quadratic to double check
#This produces a plot that plots quadratic responses if the quadratic term is < 0.05, otherwise linear plot
dev.off()
pdf("Output/Figures/germ_PC1_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC1, pch=19, col="grey60", ylab="Probability of germination", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC1) 
  model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, plotted.data)
  if(summary(model)$coefficients[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

## PLDE - quadratic PC2 germination
## HYGL - quadratic PC2 germination
## ARCA - quadratic PC2 germination
## All other species - linear PC2 germination
for (i in 1:length(specieslist)){
  print(specieslist[i])
  germquadPC2 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(germquadPC2))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/germ_PC2_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$percent_germ~plotted.data$std_PC2, pch=19, col="grey60", ylab="Probability of germination", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  model<-glmer(cbind(total_germ, total_no_germ)~std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(cbind(total_germ, total_no_germ)~std_PC2 + (1|Site/Plot), family = binomial, plotted.data)
  if(summary(model)$coefficients[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()
#POLE not converging for linear model
# poletestmod <- glmer(cbind(total_germ, total_no_germ)~std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = binomial, poledata)
# summary(poletestmod)
# poletestmod2 <- glmer(cbind(total_germ, total_no_germ)~std_PC2 + (1|Site/Plot), family = binomial, poledata)
# summary(poletestmod2)
# poletestmod2 <- glmer(cbind(total_germ, total_no_germ)~std_PC2 + (1|Site/Plot), family = binomial, poledata)
# 
# ggplot(poledata, aes(x = std_PC2, y = percent_germ))+
#   geom_jitter(alpha = 0.3)+
#   geom_smooth()+
#   theme_classic()

### Survival ###
## PLDE - quadratic survival ~ PC1
## All other species - linear PC1 survival
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survquadPC1 <- glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                           family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(survquadPC1))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/surv_PC1_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC1, pch=19, col="grey60", ylab="Probability of survival", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC1) 
  model<-glmer(surv_to_produce_seeds~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, plotted.data)
  if(summary(model)$coefficients[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

## All species - linear PC2 survival
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survquadPC2 <- glmer(surv_to_produce_seeds ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(survquadPC2))
}
#Plotting if quadratic to double check
#TRCY linear model not converging
dev.off()
pdf("Output/Figures/surv_PC2_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_PC2, pch=19, col="grey60", ylab="Probability of survival", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  model<-glmer(surv_to_produce_seeds~std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(surv_to_produce_seeds~std_PC2 + (1|Site/Plot), family = binomial, plotted.data)
  if(summary(model)$coefficients[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

### TRCY - quadratic neighbour abundance survival
## All other species linear
# #ARCA linear model not converging
# ggplot(arcadata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
#   geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
#   geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ poly(x, 2))+
#   theme_classic()
# 
# hyglsurvnamodel2 <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + (1|Site/Plot), 
#                          family = binomial, arcadata)

for (i in 1:length(specieslist)){
  print(specieslist[i])
  survquadNA <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(survquadNA))
}
#Plotting if quadratic to double check
#Added a line of code to plot * if significant 
dev.off()
pdf("Output/Figures/surv_NA_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_logp1_totalabund, pch=19, col="grey60", ylab="Probability of survival", xlab="Neighbour abundance (log plus 1, standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_logp1_totalabund) 
  model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, plotted.data)
  if(summary(model)$coefficients[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
    text(x = 2,y = 0.9,"*", cex = 10)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

### Fecundity ###
#Note that I have to use different datasets - seedmodeldata/seedarca/species.list.f
# arcaspquad1 <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, seedarca)
# summary(arcaspquad1)
# coef(summary(arcaspquad1))$cond[3,4]
#Figured out two ways to extract the p value from glmmTMB models
# 1: two steps:
#test <- coef(summary(arcaspquad1))[["cond"]]
#test[3,4]
#2: one step:
#test2<-coef(summary(arcaspquad1))$cond[3,4]

## VERO - quadratic PC1 fecundity
## All other species - linear PC1 fecundity
for (i in 1:length(specieslist)){
  print(specieslist[i])
  seedquadPC1 <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                       family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(seedquadPC1))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/SP_PC1_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_PC1, pch=19, col="grey60", ylab="Number of seeds produced", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC1) 
  model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  model2<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, plotted.data)
  if(coef(summary(model))$cond[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

## All species - linear PC2 fecundity
for (i in 1:length(specieslist)){
  print(specieslist[i])
  seedquadPC2 <- glmmTMB(No_viable_seeds_grouped ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                         family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(seedquadPC2))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/SP_PC2_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_PC2, pch=19, col="grey60", ylab="Number of seeds produced", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  model<-glmmTMB(No_viable_seeds_grouped~std_PC2 + I(std_PC2^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  model2<-glmmTMB(No_viable_seeds_grouped~std_PC2 + (1|Site/Plot), family = nbinom2, plotted.data)
  if(coef(summary(model))$cond[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

## All species - linear neighbour abundance fecundity
for (i in 1:length(specieslist)){
  print(specieslist[i])
  seedquadPC2 <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), 
                         family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(seedquadPC2))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/SP_NA_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_logp1_totalabund, pch=19, col="grey60", ylab="Number of seeds produced", xlab="Neighbour abundance (log plus 1, standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_logp1_totalabund) 
  model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  model2<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, plotted.data)
  if(coef(summary(model))$cond[3,4]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

### Population growth rate, lambda ###

## PEAI - quadratic PC1 lambda
## All other species - linear PC1 lambda
for (i in 1:length(specieslist)){
  print(specieslist[i])
  lambdaquadPC1 <- lmer(log_lambda_p1 ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                       data = filter(popdata, Species == specieslist[i]))
  print(summary(lambdaquadPC1))
}
#Plotting if quadratic to double check
#This produces a plot that plots quadratic responses if the quadratic term is < 0.05, otherwise linear plot
#summary(mod1)$coefficients[3,5] - these are the coefficient coordinates we want

dev.off()
pdf("Output/Figures/lambda_PC1_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.l)){
  plotted.data<-as.data.frame(species.list.l[i])
  plot(plotted.data$log_lambda_p1~plotted.data$std_PC1, pch=19, col="grey60", ylab="Population growth rate (logged plus 1)", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC1) 
  model<-lmer(log_lambda_p1~std_PC1 + I(std_PC1^2) + (1|Site/Plot), plotted.data)
  model2<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), plotted.data)
  if(summary(model)$coefficients[3,5]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

#This is what peai should look like, significant quadratic term
# plot(log_lambda_p1 ~ std_PC1, lambdapeai)
# model <- lmer(log_lambda_p1~std_PC1 + I(std_PC1^2) + (1|Site/Plot), lambdapeai)
# x_to_plot<-seq.func(lambdapeai$std_PC1) 
# preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
# plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
# plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

## All species - linear PC2 lambda
for (i in 1:length(specieslist)){
  print(specieslist[i])
  lambdaquadPC2 <- lmer(log_lambda_p1 ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                        data = filter(popdata, Species == specieslist[i]))
  print(summary(lambdaquadPC2))
}
#Plotting if quadratic to double check

dev.off()
pdf("Output/Figures/lambda_PC2_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.l)){
  plotted.data<-as.data.frame(species.list.l[i])
  plot(plotted.data$log_lambda_p1~plotted.data$std_PC2, pch=19, col="grey60", ylab="Population growth rate (logged plus 1)", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  model<-lmer(log_lambda_p1~std_PC2 + I(std_PC2^2) + (1|Site/Plot), plotted.data)
  model2<-lmer(log_lambda_p1~std_PC2 + (1|Site/Plot), plotted.data)
  if(summary(model)$coefficients[3,5]<0.05){
    preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
    plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }else{
    preddata <- with(model2, data.frame(1, x_to_plot))
    plotted.pred <- glmm.predict(mod = model2, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
    plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
  }
}
dev.off()

#Can't look at lambda ~ continuous neighbour abundance

### Making a table with quadratic significance values ####
## Germination
#Create table to put these values in
germquadtable <- matrix(ncol=2, nrow = 9)
#define column names and row names of matrix
colnames(germquadtable) <- c('PC1', 'PC2')
rownames(germquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
#convert matrix to table 
germquadtable <- as.data.frame(germquadtable)

## Extracting values for PC1
for (i in 1:length(specieslist)){
  model <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                 family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  nam <- paste0("germPC1pvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,4])
}
##Replace values with appropriate ones
germquadtable[1,1] <- germPC1pvalueARCA
germquadtable[2,1] <- germPC1pvalueHYGL
germquadtable[3,1] <- germPC1pvalueLARO
germquadtable[4,1] <- germPC1pvaluePEAI
germquadtable[5,1] <- germPC1pvaluePLDE
germquadtable[6,1] <- germPC1pvaluePOLE
germquadtable[7,1] <- germPC1pvalueTRCY
germquadtable[8,1] <- germPC1pvalueTROR
germquadtable[9,1] <- germPC1pvalueVERO
## Extracting values for PC2
for (i in 1:length(specieslist)){
  model <- glmer(cbind(total_germ, total_no_germ) ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                 family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  nam <- paste0("germPC2pvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,4])
}
##Replace values with appropriate ones
germquadtable[1,2] <- germPC2pvalueARCA
germquadtable[2,2] <- germPC2pvalueHYGL
germquadtable[3,2] <- germPC2pvalueLARO
germquadtable[4,2] <- germPC2pvaluePEAI
germquadtable[5,2] <- germPC2pvaluePLDE
germquadtable[6,2] <- germPC2pvaluePOLE
germquadtable[7,2] <- germPC2pvalueTRCY
germquadtable[8,2] <- germPC2pvalueTROR
germquadtable[9,2] <- germPC2pvalueVERO

#No neighbour abundance analysis here

## Survival
#Create table to put these values in
survquadtable <- matrix(ncol=3, nrow = 9)
#define column names and row names of matrix
colnames(survquadtable) <- c('PC1', 'PC2', 'N.A.')
rownames(survquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
#convert matrix to table 
survquadtable <- as.data.frame(survquadtable)

## Extracting values for PC1
for (i in 1:length(specieslist)){
  model <- glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                 family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  nam <- paste0("survPC1pvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,4])
}
##Replace values with appropriate ones
survquadtable[1,1] <- survPC1pvalueARCA
survquadtable[2,1] <- survPC1pvalueHYGL
survquadtable[3,1] <- survPC1pvalueLARO
survquadtable[4,1] <- survPC1pvaluePEAI
survquadtable[5,1] <- survPC1pvaluePLDE
survquadtable[6,1] <- survPC1pvaluePOLE
survquadtable[7,1] <- survPC1pvalueTRCY
survquadtable[8,1] <- survPC1pvalueTROR
survquadtable[9,1] <- survPC1pvalueVERO
## Extracting values for PC2
for (i in 1:length(specieslist)){
  model <- glmer(surv_to_produce_seeds ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                 family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  nam <- paste0("survPC2pvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,4])
}
##Replace values with appropriate ones
survquadtable[1,2] <- survPC2pvalueARCA
survquadtable[2,2] <- survPC2pvalueHYGL
survquadtable[3,2] <- survPC2pvalueLARO
survquadtable[4,2] <- survPC2pvaluePEAI
survquadtable[5,2] <- survPC2pvaluePLDE
survquadtable[6,2] <- survPC2pvaluePOLE
survquadtable[7,2] <- survPC2pvalueTRCY
survquadtable[8,2] <- survPC2pvalueTROR
survquadtable[9,2] <- survPC2pvalueVERO

## Extracting values for neighbour abundance
#Note that HYGL doesn't converge here, not quadratic
for (i in 1:length(specieslist)){
  model <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), 
                 family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  nam <- paste0("survNApvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,4])
}
##Replace values with appropriate ones
survquadtable[1,3] <- survNApvalueARCA
survquadtable[2,3] <- survNApvalueHYGL
survquadtable[3,3] <- survNApvalueLARO
survquadtable[4,3] <- survNApvaluePEAI
survquadtable[5,3] <- survNApvaluePLDE
survquadtable[6,3] <- survNApvaluePOLE
survquadtable[7,3] <- survNApvalueTRCY
survquadtable[8,3] <- survNApvalueTROR
survquadtable[9,3] <- survNApvalueVERO

## Fecundity
#Create table to put these values in
fecundityquadtable <- matrix(ncol=3, nrow = 9)
#define column names and row names of matrix
colnames(fecundityquadtable) <- c('PC1', 'PC2', 'N.A.')
rownames(fecundityquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
#convert matrix to table 
fecundityquadtable <- as.data.frame(fecundityquadtable)

## Extracting values for PC1
for (i in 1:length(specieslist)){
  model <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                   family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  nam <- paste0("seedPC1pvalue", specieslist[i])
  assign(nam, coef(summary(model))$cond[3,4])
}
##Replace values with appropriate ones
fecundityquadtable[1,1] <- seedPC1pvalueARCA
fecundityquadtable[2,1] <- seedPC1pvalueHYGL
fecundityquadtable[3,1] <- seedPC1pvalueLARO
fecundityquadtable[4,1] <- seedPC1pvaluePEAI
fecundityquadtable[5,1] <- seedPC1pvaluePLDE
fecundityquadtable[6,1] <- seedPC1pvaluePOLE
fecundityquadtable[7,1] <- seedPC1pvalueTRCY
fecundityquadtable[8,1] <- seedPC1pvalueTROR
fecundityquadtable[9,1] <- seedPC1pvalueVERO
## Extracting values for PC2
for (i in 1:length(specieslist)){
  model <- glmmTMB(No_viable_seeds_grouped ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                   family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  nam <- paste0("seedPC2pvalue", specieslist[i])
  assign(nam, coef(summary(model))$cond[3,4])
}
##Replace values with appropriate ones
fecundityquadtable[1,2] <- seedPC2pvalueARCA
fecundityquadtable[2,2] <- seedPC2pvalueHYGL
fecundityquadtable[3,2] <- seedPC2pvalueLARO
fecundityquadtable[4,2] <- seedPC2pvaluePEAI
fecundityquadtable[5,2] <- seedPC2pvaluePLDE
fecundityquadtable[6,2] <- seedPC2pvaluePOLE
fecundityquadtable[7,2] <- seedPC2pvalueTRCY
fecundityquadtable[8,2] <- seedPC2pvalueTROR
fecundityquadtable[9,2] <- seedPC2pvalueVERO

## Extracting values for neighbour abundance
for (i in 1:length(specieslist)){
  model <- glmmTMB(No_viable_seeds_grouped ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), 
                         family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  nam <- paste0("seedNApvalue", specieslist[i])
  assign(nam, coef(summary(model))$cond[3,4])
}
##Replace values with appropriate ones
fecundityquadtable[1,3] <- seedNApvalueARCA
fecundityquadtable[2,3] <- seedNApvalueHYGL
fecundityquadtable[3,3] <- seedNApvalueLARO
fecundityquadtable[4,3] <- seedNApvaluePEAI
fecundityquadtable[5,3] <- seedNApvaluePLDE
fecundityquadtable[6,3] <- seedNApvaluePOLE
fecundityquadtable[7,3] <- seedNApvalueTRCY
fecundityquadtable[8,3] <- seedNApvalueTROR
fecundityquadtable[9,3] <- seedNApvalueVERO

## Lambda
#Create table to put these values in
lambdaquadtable <- matrix(ncol=2, nrow = 9)
#define column names and row names of matrix
colnames(lambdaquadtable) <- c('PC1', 'PC2')
rownames(lambdaquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
#convert matrix to table 
lambdaquadtable <- as.data.frame(lambdaquadtable)

## Extracting values for PC1
for (i in 1:length(specieslist)){
  model <- lmer(log_lambda_p1 ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
                   data = filter(popdata, Species == specieslist[i]))
  nam <- paste0("lambdaPC1pvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,5])
}
##Replace values with appropriate ones
lambdaquadtable[1,1] <- lambdaPC1pvalueARCA
lambdaquadtable[2,1] <- lambdaPC1pvalueHYGL
lambdaquadtable[3,1] <- lambdaPC1pvalueLARO
lambdaquadtable[4,1] <- lambdaPC1pvaluePEAI
lambdaquadtable[5,1] <- lambdaPC1pvaluePLDE
lambdaquadtable[6,1] <- lambdaPC1pvaluePOLE
lambdaquadtable[7,1] <- lambdaPC1pvalueTRCY
lambdaquadtable[8,1] <- lambdaPC1pvalueTROR
lambdaquadtable[9,1] <- lambdaPC1pvalueVERO
## Extracting values for PC2
for (i in 1:length(specieslist)){
  model <- lmer(log_lambda_p1 ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                data = filter(popdata, Species == specieslist[i]))
  nam <- paste0("lambdaPC2pvalue", specieslist[i])
  assign(nam, summary(model)$coefficients[3,5])
}
##Replace values with appropriate ones
lambdaquadtable[1,2] <- lambdaPC2pvalueARCA
lambdaquadtable[2,2] <- lambdaPC2pvalueHYGL
lambdaquadtable[3,2] <- lambdaPC2pvalueLARO
lambdaquadtable[4,2] <- lambdaPC2pvaluePEAI
lambdaquadtable[5,2] <- lambdaPC2pvaluePLDE
lambdaquadtable[6,2] <- lambdaPC2pvaluePOLE
lambdaquadtable[7,2] <- lambdaPC2pvalueTRCY
lambdaquadtable[8,2] <- lambdaPC2pvalueTROR
lambdaquadtable[9,2] <- lambdaPC2pvalueVERO


### Binding germ, surv and fecund tables together ###
vitalquadpvalues <- cbind(germquadtable, survquadtable, fecundityquadtable, lambdaquadtable)

#Plotting with kableR
# To use italics, <i> and </i>. To use bold, <b> and </b>.
#Or digits = 2 to change rounding of everything to 2, want to change this so PEAI is 4.0e-5 and TRCY is 3.0e-2 and all others are 2 dp
## Can use scales::pvalue or p_value scales to change to < 0.01, < 0.001 etc. do this*
#Not sure how because my dataframe has duplicate names #  mutate(scales::pvalue[,1]) %>%
#Could manually rename couple of values that are too long? # test[4,1] <- "<0.001". Keeps bolding but converts them to character so need to figure that out
#test[2,5] <- "na"
#test$N.A. <- as.numeric(test$N.A.)

vitalquadpvalues %>% 
  kbl(caption = "<b>Supplementary 1</b>. Model output <i>p</i>-values for quadratic terms of fixed effects. N.A. denotes neighbour abundance. <b>Bolded</b> values indicate a <i>p</i>-value of <0.05.", digits = c(5, 2, 2, 2, 3, 2, 2, 2, 2, 2)) %>%
  kable_classic(full_width = F, html_font = "Times") %>%
  column_spec(1, italic = T) %>%
  #row_spec(0, bold = T) %>%
  column_spec(2, bold = ifelse(vitalquadpvalues[[1]] <0.05, TRUE, FALSE)) %>%
  column_spec(3, bold = ifelse(vitalquadpvalues[[2]] <0.05, TRUE, FALSE)) %>%
  column_spec(4, bold = ifelse(vitalquadpvalues[[3]] <0.05, TRUE, FALSE)) %>%
  column_spec(5, bold = ifelse(vitalquadpvalues[[4]] <0.05, TRUE, FALSE)) %>%
  column_spec(6, bold = ifelse(vitalquadpvalues[[5]] <0.05, TRUE, FALSE)) %>%
  column_spec(7, bold = ifelse(vitalquadpvalues[[6]] <0.05, TRUE, FALSE)) %>%
  column_spec(8, bold = ifelse(vitalquadpvalues[[7]] <0.05, TRUE, FALSE)) %>%
  column_spec(9, bold = ifelse(vitalquadpvalues[[8]] <0.05, TRUE, FALSE)) %>%
  column_spec(10, bold = ifelse(vitalquadpvalues[[9]] <0.05, TRUE, FALSE)) %>%
  column_spec(11, bold = ifelse(vitalquadpvalues[[10]] <0.05, TRUE, FALSE)) %>%
  add_header_above(c("", "Germination" = 2, "Survival" = 3, "Fecundity" = 3, "Lambda" = 2))
  #add_header_above(c("", "p value quadratic term" = 8))
#Struggling to save this using save_kable and as_image() atm.
#Can copy it from the Viewer using copy to clipboard, maintain aspect ratio, first value 500

#### Modelling vital rates in response to all factors #### 
### Germination ###
## ARCA - has quadratic germ ~ PC2
arcagermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                          family = binomial, arcadata)
arcagermmod1dharma <- simulateResiduals(arcagermfinalmod)
plot(arcagermmod1dharma)
#Not great
vif(arcagermfinalmod)
summary(arcagermfinalmod)
testDispersion(arcagermfinalmod)

# ggplot(arcadata, aes(x = std_PC1, y = percent_germ))+
#   geom_jitter(alpha=0.4)+
#   geom_smooth()+
#   theme_classic()

## HYGL- has quadratic germ ~ PC2
hyglgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
                          family = binomial, hygldata)
hyglgermmod1dharma <- simulateResiduals(hyglgermfinalmod)
plot(hyglgermmod1dharma)
#Not great residuals
vif(hyglgermfinalmod)
summary(hyglgermfinalmod)

## LARO - has quadratic germ ~ PC1
larogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 +
                            I(std_PC1^2) + (1|Site/Plot), 
                          family = binomial, larodata)
larogermmod1dharma <- simulateResiduals(larogermfinalmod)
plot(larogermmod1dharma)
#not great
vif(larogermfinalmod)
summary(larogermfinalmod)

## PEAI - has quadratic germ ~ PC1
peaigermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + 
                            I(std_PC1^2) + (1|Site/Plot), 
                          family = binomial, peaidata)
peaigermmod1dharma <- simulateResiduals(peaigermfinalmod)
plot(peaigermmod1dharma)
#Not great residuals, dispersion test and KS test significant
vif(peaigermfinalmod)
summary(peaigermfinalmod)

## PLDE - has quadratic germ ~ PC2
pldegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + 
                            I(std_PC2^2)  + (1|Site/Plot), 
                          family = binomial, pldedata)
pldegermmod1dharma <- simulateResiduals(pldegermfinalmod)
plot(pldegermmod1dharma)
#Not too bad residuals, KS, dispersion and outlier tests significant though
vif(pldegermfinalmod)
summary(pldegermfinalmod)

## POLE
polegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, poledata)
polegermmod1dharma <- simulateResiduals(polegermfinalmod)
plot(polegermmod1dharma)
#good
vif(polegermfinalmod)
summary(polegermfinalmod)

## TRCY
trcygermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, trcydata)
trcygermmod1dharma <- simulateResiduals(trcygermfinalmod)
plot(trcygermmod1dharma)
#terrible
vif(trcygermfinalmod)
summary(trcygermfinalmod)

## TROR
trorgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, trordata)
trorgermmod1dharma <- simulateResiduals(trorgermfinalmod)
plot(trorgermmod1dharma)
#Not great residuals, KS test significant
vif(trorgermfinalmod)
summary(trorgermfinalmod)

## VERO
verogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, verodata)
verogermmod1dharma <- simulateResiduals(verogermfinalmod)
plot(verogermmod1dharma)
#Not great residuals, KS and dispersion tests significant
vif(verogermfinalmod)
summary(verogermfinalmod)

### Survival ###
## ARCA - quadratic surv ~ NA
arcasurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + 
                        I(std_logp1_totalabund^2) + Treatment:I(std_logp1_totalabund^2) + std_PC1:I(std_logp1_totalabund^2) + (1|Site/Plot), 
                      family = binomial, arcadata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
arcasurvmod1dharma <- simulateResiduals(arcasurvmod1)
plot(arcasurvmod1dharma)
#good
vif(arcasurvmod1)
summary(arcasurvmod1)
# Model simplification step - Significant Dry:PC1, removing all others
arcasurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            I(std_logp1_totalabund^2) + Treatment:std_PC1 + (1|Site/Plot), 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), arcadata)
arcasurvfinalmoddharma <- simulateResiduals(arcasurvfinalmod)
plot(arcasurvfinalmoddharma)
#good 
summary(arcasurvfinalmod)
#NA^2 is not significant but NA is - trying to plot it to see if NA truly is sig
# with(arcadata, plot(jitter(surv_to_produce_seeds, amount = 0.05) ~ std_logp1_totalabund))
# curve(plogis(cbind(1, 0, 0, 0, 0, x, 0, 0, 0*0, 0*0)%*%fixef(arcasurvfinalmod)), add = TRUE)
# ggplot(arcadata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
#   geom_jitter()+
#   geom_smooth()+
#   theme_classic()

## HYGL - needs optimiser to  converge
hyglsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
hyglsurvmod1dharma <- simulateResiduals(hyglsurvmod1)
plot(hyglsurvmod1dharma)
#good
vif(hyglsurvmod1)
summary(hyglsurvmod1)
# Model simplification step - Significant Dry:PC1, Wet:PC1, removing others
hyglsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            Treatment:std_PC1 + (1|Site/Plot), 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
summary(hyglsurvfinalmod)

#### Meeting John 22/06/22
#### Plotting interaction between two continuous variables
# PC1 is x 
#with(hygldata, plot(jitter(surv_to_produce_seeds, amount = 0.05) ~ std_PC1, col = ifelse(std_logp1_totalabund>0, "red", "blue")))
#curve(plogis(cbind(1, 0, 0, x, 0, 1, 0, 0*x, 0*x, 0*1, 0*1, x*1)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "red")
#curve(plogis(cbind(1, 0, 0, x, 0, -1, 0, 0*x, 0*x, 0*-1, 0*-1, x*-1)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "blue")

## LARO
larosurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, larodata)
larosurvmod1dharma <- simulateResiduals(larosurvmod1)
plot(larosurvmod1dharma)
#good
vif(larosurvmod1)
summary(larosurvmod1)
# Model simplification step - nothing significant, removing all
larosurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = binomial, larodata)
larosurvfinalmoddharma <- simulateResiduals(larosurvfinalmod)
plot(larosurvfinalmoddharma)
#good
summary(larosurvfinalmod)

# Testing why surv ~ NA isn't significant for laro even though it looks it on plots
# larosimplemod <- glm(surv_to_produce_seeds ~ std_logp1_totalabund, family = "binomial", larodata)
# summary(larosimplemod)
# #just
# larosimplemod2 <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + (1|Site/Plot), family = "binomial", larodata)
# summary(larosimplemod2)
# #random make a small difference

## PEAI 
peaisurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), peaidata)
peaisurvmod1dharma <- simulateResiduals(peaisurvmod1)
plot(peaisurvmod1dharma)
#good
vif(peaisurvmod1)
summary(peaisurvmod1)
# Model simplification step - significant PC1:NA, removing others
peaisurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                          family = binomial, peaidata)
peaisurvfinalmoddharma <- simulateResiduals(peaisurvfinalmod)
plot(peaisurvfinalmoddharma)
#good
summary(peaisurvfinalmod)

## PLDE - quadratic PC1
pldesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + 
                        I(std_PC1^2) + Treatment:I(std_PC1^2) + I(std_PC1^2):std_logp1_totalabund + (1|Site/Plot), 
                        family = binomial, pldedata)
pldesurvmod1dharma <- simulateResiduals(pldesurvmod1)
plot(pldesurvmod1dharma)
#okay
vif(pldesurvmod1)
summary(pldesurvmod1)
# Model simplification step - Nothing significant, removing all
pldesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            I(std_PC1^2) + (1|Site/Plot), 
                          family = binomial, pldedata)
pldesurvfinalmoddharma <- simulateResiduals(pldesurvfinalmod)
plot(pldesurvfinalmoddharma)
#good
summary(pldesurvfinalmod)

## POLE
#problem here?
polesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), poledata)
polesurvmod1dharma <- simulateResiduals(polesurvmod1)
plot(polesurvmod1dharma)
#okay
vif(polesurvmod1)
summary(polesurvmod1)
# Model simplification step - Nothing significant, removing all
polesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, poledata)
polesurvfinalmoddharma <- simulateResiduals(polesurvfinalmod)
plot(polesurvfinalmoddharma)
#okay
summary(polesurvfinalmod)

## TRCY - quadratic surv ~ NA
trcysurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + 
                        I(std_logp1_totalabund^2) + Treatment:I(std_logp1_totalabund^2) + std_PC1:I(std_logp1_totalabund^2) + (1|Site/Plot), 
                      family = binomial, trcydata)
trcysurvmod1dharma <- simulateResiduals(trcysurvmod1)
plot(trcysurvmod1dharma)
#not too bad
vif(trcysurvmod1)
summary(trcysurvmod1)
# Model simplification step - Significant Wet:PC1, Wet:NA^2, removing all others
#Where quadratic terms are significant, only keeping quadratic term (not main term as well, only adjusting curvature, not slope too)
trcysurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            I(std_logp1_totalabund^2) + Treatment:std_PC1 + Treatment:I(std_logp1_totalabund^2) + (1|Site/Plot), 
                          family = binomial, trcydata)
trcysurvfinalmoddharma <- simulateResiduals(trcysurvfinalmod)
plot(trcysurvfinalmoddharma)
vif(trcysurvfinalmod)
#not too bad
summary(trcysurvfinalmod)

## TROR - needs optimiser to converge
trorsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trordata)
trorsurvmod1dharma <- simulateResiduals(trorsurvmod1)
plot(trorsurvmod1dharma)
#good
vif(trorsurvmod1)
summary(trorsurvmod1)
# Model simplification step - Significant Wet:PC1, removing all others
trorsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            Treatment:std_PC1 + (1|Site/Plot), 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trordata)
trorsurvfinalmoddharma <- simulateResiduals(trorsurvfinalmod)
plot(trorsurvfinalmoddharma)
#not too bad
summary(trorsurvfinalmod)

## VERO
verosurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, verodata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
verosurvmod1dharma <- simulateResiduals(verosurvmod1)
plot(verosurvmod1dharma)
#good
vif(verosurvmod1)
summary(verosurvmod1)
# Model simplification step - Significant Wet:PC1, removing all others
#Still needs optimiser to converge
verosurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            Treatment:std_PC1 + (1|Site/Plot), 
                          family = binomial, verodata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
verosurvfinalmoddharma <- simulateResiduals(verosurvfinalmod)
plot(verosurvfinalmoddharma)
#good
summary(verosurvfinalmod)

### Fecundity ###
## ARCA
arcaseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedarca)
arcaseedmod1dharma <- simulateResiduals(arcaseedmod1)
plot(arcaseedmod1dharma)
#good
summary(arcaseedmod1)
# Model simplification step - Significant Wet:NA, removing all others
arcaseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                              Treatment:std_logp1_totalabund + (1|Site/Plot), 
                            family = nbinom2, seedarca)
arcaseedfinalmoddharma <- simulateResiduals(arcaseedfinalmod)
plot(arcaseedfinalmoddharma)
#good
summary(arcaseedfinalmod)
r.squaredGLMM(arcaseedfinalmod)
# ### PC2 and neighbour abund are correlated, so seeing how model differs when I remove each
# arcaseedfinalmod2 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_logp1_totalabund + Dodder01 + 
#                               Treatment:std_logp1_totalabund + (1|Site/Plot), 
#                             family = nbinom2, seedarca)
# arcaseedfinalmoddharma2 <- simulateResiduals(arcaseedfinalmod2)
# plot(arcaseedfinalmoddharma2)
# #good
# summary(arcaseedfinalmod2)
# r.squaredGLMM(arcaseedfinalmod2)
# #
# arcaseedfinalmod3 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + Dodder01 + 
#                                (1|Site/Plot), 
#                              family = nbinom2, seedarca)
# arcaseedfinalmoddharma3 <- simulateResiduals(arcaseedfinalmod3)
# plot(arcaseedfinalmoddharma3)
# #good
# summary(arcaseedfinalmod3)
# r.squaredGLMM(arcaseedfinalmod3)
# AIC(arcaseedfinalmod, arcaseedfinalmod2, arcaseedfinalmod3)

## HYGL
hyglseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedhygl)
hyglseedmod1dharma <- simulateResiduals(hyglseedmod1)
plot(hyglseedmod1dharma)
#good
summary(hyglseedmod1)
# Model simplification step - Nothing significant, removing all
hyglseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedhygl)
hyglseedfinalmoddharma <- simulateResiduals(hyglseedfinalmod)
plot(hyglseedfinalmoddharma)
#good
summary(hyglseedfinalmod)

## LARO
laroseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                          family = nbinom2, seedlaro)
laroseedmod1dharma <- simulateResiduals(laroseedmod1)
plot(laroseedmod1dharma)
#good
summary(laroseedmod1)
# Model simplification step - Significant PC1:NA, removing all others
laroseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                              std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                            family = nbinom2, seedlaro)
laroseedfinalmoddharma <- simulateResiduals(laroseedfinalmod)
plot(laroseedfinalmoddharma)
#good
summary(laroseedfinalmod)

## PEAI
peaiseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedpeai)
peaiseedmod1dharma <- simulateResiduals(peaiseedmod1)
plot(peaiseedmod1dharma)
#good
summary(peaiseedmod1)
# Model simplification step - Significant Wet:NA, removing others
peaiseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                              Treatment:std_logp1_totalabund + (1|Site/Plot), seedpeai)
peaiseedfinalmoddharma <- simulateResiduals(peaiseedfinalmod)
plot(peaiseedfinalmoddharma)
#good
summary(peaiseedfinalmod)

## PLDE
pldeseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedplde)
pldeseedmod1dharma <- simulateResiduals(pldeseedmod1)
plot(pldeseedmod1dharma)
#good
summary(pldeseedmod1)
# Model simplification step - Nothing significant, removing all
pldeseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedplde)
pldeseedfinalmoddharma <- simulateResiduals(pldeseedfinalmod)
plot(pldeseedfinalmoddharma)
#good
summary(pldeseedfinalmod)

## POLE - can't use this data, only 19 observations 
#Not converging with glmer.nb or glmm TMB
# poleseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
#                           Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
#                         family = nbinom2, seedpole)
# # Model simplification step - Removing all interactions, seeing if simple model will fit
# poleseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
#                             family = nbinom2, seedpole)
# poleseedfinalmoddharma <- simulateResiduals(poleseedfinalmod)
# plot(poleseedfinalmoddharma)
#terrible residuals

## TRCY
trcyseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedtrcy)
trcyseedmod1dharma <- simulateResiduals(trcyseedmod1)
plot(trcyseedmod1dharma)
#good
summary(trcyseedmod1)
# Model simplification step - Nothing significant, removing all
trcyseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                             + (1|Site/Plot), 
                            family = nbinom2, seedtrcy)
trcyseedfinalmoddharma <- simulateResiduals(trcyseedfinalmod)
plot(trcyseedfinalmoddharma)
#good
summary(trcyseedfinalmod)

## TROR
trorseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedtror)
trorseedmod1dharma <- simulateResiduals(trorseedmod1)
plot(trorseedmod1dharma)
#good
summary(trorseedmod1)
# Model simplification step - Significant Dry:PC1, removing all others
trorseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                              Treatment:std_PC1 + (1|Site/Plot), 
                            family = nbinom2, seedtror)
trorseedfinalmoddharma <- simulateResiduals(trorseedfinalmod)
plot(trorseedfinalmoddharma)
#not too bad
summary(trorseedfinalmod)

## VERO - significant fecundity ~ PC1
veroseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + 
                          I(std_PC1^2) + Treatment:I(std_PC1^2) + std_logp1_totalabund:I(std_PC1^2) + (1|Site/Plot), 
                        family = nbinom2, seedvero)
veroseedmod1dharma <- simulateResiduals(veroseedmod1)
plot(veroseedmod1dharma)
#good
summary(veroseedmod1)
# Model simplification step - Nothing significant, removing all
veroseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                              I(std_PC1^2) + (1|Site/Plot), 
                            family = nbinom2, seedvero)
veroseedfinalmoddharma <- simulateResiduals(veroseedfinalmod)
plot(veroseedfinalmoddharma)
#good
summary(veroseedfinalmod)

### Lambda ###
## ARCA
arcalambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdamod1dharma <- simulateResiduals(arcalambdamod1)
plot(arcalambdamod1dharma)
#good
vif(arcalambdamod1)
summary(arcalambdamod1)
#Model simplification step - No significant interactions, removing all interactions
arcalambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdafinalmoddharma <- simulateResiduals(arcalambdafinalmod)
plot(arcalambdafinalmoddharma)
#good
summary(arcalambdafinalmod)

## HYGL
hygllambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdamod1dharma <- simulateResiduals(hygllambdamod1)
plot(hygllambdamod1dharma)
#okay
vif(hygllambdamod1)
summary(hygllambdamod1)
#Model simplification step - No interactions significant, removing all
hygllambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdafinalmoddharma <- simulateResiduals(hygllambdafinalmod)
plot(hygllambdafinalmoddharma)
#okay
summary(hygllambdafinalmod)

##LARO
larolambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdalaro)
larolambdamod1dharma <- simulateResiduals(larolambdamod1)
plot(larolambdamod1dharma)
#okay
vif(larolambdamod1)
summary(larolambdamod1)
#Model simplification step - Significant PC1:neighbours, removing others
larolambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdalaro)
larolambdafinalmoddharma <- simulateResiduals(larolambdafinalmod)
plot(larolambdafinalmoddharma)
#not great
summary(larolambdafinalmod)

##PEAI - has quadratic lambda ~ PC1
peailambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + 
                         I(std_PC1^2) + Treatment:I(std_PC1^2) + I(std_PC1^2):neighbours01 + (1|Site/Plot), lambdapeai)
peailambdamod1dharma <- simulateResiduals(peailambdamod1)
plot(peailambdamod1dharma)
#okay, dispersion test significant
vif(peailambdamod1)
summary(peailambdamod1)
#Model simplification step - Significant Wet:PC1, removing all other interactions
peailambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + 
                             I(std_PC1^2) + Treatment:std_PC1 + (1|Site/Plot), lambdapeai)
peailambdafinalmoddharma <- simulateResiduals(peailambdafinalmod)
plot(peailambdafinalmoddharma)
#good
summary(peailambdafinalmod)

##PLDE
pldelambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdamod1dharma <- simulateResiduals(pldelambdamod1)
plot(pldelambdamod1dharma)
#very bad residuals
vif(pldelambdamod1)
summary(pldelambdamod1)
#Model simplification step - No significant interactions, removing all interactions
pldelambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdafinalmoddharma <- simulateResiduals(pldelambdafinalmod)
plot(pldelambdafinalmoddharma)
#bad residuals
summary(pldelambdafinalmod)
#with(lambdaplde, plot(log_lambda_p1 ~ std_PC1))

##POLE
polelambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdapole)
polelambdamod1dharma <- simulateResiduals(polelambdamod1)
plot(polelambdamod1dharma)
#not good
vif(polelambdamod1)
summary(polelambdamod1)
#Model simplification step - No significant interactions, removing all interactions
polelambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + (1|Site/Plot), lambdapole)
polelambdafinalmoddharma <- simulateResiduals(polelambdafinalmod)
plot(polelambdafinalmoddharma)
#not good - not enough data? 18 points. Residuals bad
summary(polelambdafinalmod)

##TRCY
trcylambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdamod1dharma <- simulateResiduals(trcylambdamod1)
plot(trcylambdamod1dharma)
#good
vif(trcylambdamod1)
summary(trcylambdamod1)
#Model simplification step - Significant PC1:neighbours, removing all other interactions
trcylambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + 
                             std_PC1:neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdafinalmoddharma <- simulateResiduals(trcylambdafinalmod)
plot(trcylambdafinalmoddharma)
#not bad
summary(trcylambdafinalmod)

##TROR
trorlambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdatror)
trorlambdamod1dharma <- simulateResiduals(trorlambdamod1)
plot(trorlambdamod1dharma)
#good
vif(trorlambdamod1)
summary(trorlambdamod1)
#Model simplification step - No significant interactions, removing all interactions
trorlambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + (1|Site/Plot), lambdatror)
trorlambdafinalmoddharma <- simulateResiduals(trorlambdafinalmod)
plot(trorlambdafinalmoddharma)
#good
summary(trorlambdafinalmod)

##VERO
verolambdamod1 <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 +
                         Treatment:std_PC1 + Treatment:neighbours01 + std_PC1:neighbours01 + (1|Site/Plot), lambdavero)
verolambdamod1dharma <- simulateResiduals(verolambdamod1)
plot(verolambdamod1dharma)
#good
vif(verolambdamod1)
summary(verolambdamod1)
#Model simplification step - Significant PC1:neighbours, removing all other interactions
verolambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + 
                             std_PC1:neighbours01 + (1|Site/Plot), lambdavero)
verolambdafinalmoddharma <- simulateResiduals(verolambdafinalmod)
plot(verolambdafinalmoddharma)
#good
summary(verolambdafinalmod)

#### Plotting vital rate responses to continuous and categorical neighbours ####
ggplot(popdata, aes(x = neighbours01, y = log_lambda_p1))+
  geom_boxplot()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)
#Fecundity response to categorical neighbours - PLOT level
ggplot(popdata, aes(x = neighbours01, y = log(plot_fecundity+1)))+
  geom_boxplot()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)
#Fecundity response to categorical neighbours - individual level
#Need to fix neighbour NAs - 20 or so.
ggplot(data = seedmodeldata %>% filter(complete.cases(NeighboursYN)), aes(x = NeighboursYN, y = No_viable_seeds_grouped+1))+
  geom_boxplot()+
  scale_y_log10()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)

#Survival response to categorical neighbours - PLOT level
# Informed by these models, used for calculating lambda: 
#e.g. peaipopsurvmodel <- glmer(surv_to_produce_seeds ~ NeighboursYN + (NeighboursYN|Site/Plot),family = binomial, peaidata)
ggplot(popdata, aes(x = neighbours01, y = plot_survival))+
  geom_boxplot()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)
  
#Survival response to categorical neighbours - individual level
#Can't plot this with a boxplot or curve! Need catplot
# vitaldata %>% filter(!is.na(NeighboursYN)) %>% ggplot(aes(x = NeighboursYN, y = surv_to_produce_seeds, group = NeighboursYN))+
#   #geom_smooth(method = "glm", method.args=list(family="binomial"))+
#   #geom_boxplot(aes(group = NeighboursYN))+
#   geom_jitter(alpha = 0.2, width = 0.05, height = 0.05)+
#   theme_classic()+
#   my_theme+
#   facet_wrap(~Species)
### cat plots of survival from package interactions
library(interactions)
arcacatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, arcadata)
hyglcatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, hygldata)
larocatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, larodata)
peaicatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, peaidata)
pldecatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, pldedata)
polecatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, poledata)
trcycatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, trcydata)
trorcatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, trordata)
verocatmod <- glmer(surv_to_produce_seeds ~ NeighboursYN + (1|Site/Plot), family = binomial, verodata)

dev.off()
pdf("Output/Figures/catplot_survival_nbh.pdf", width=21, height=21)
a <- cat_plot(arcacatmod, pred = NeighboursYN, data = arcadata, main = "ARCA", x.label = "")
b <- cat_plot(hyglcatmod, pred = NeighboursYN, data = hygldata, main = "HYGL", x.label = "")
c <- cat_plot(larocatmod, pred = NeighboursYN, data = larodata, main = "LARO", x.label = "")
d <- cat_plot(peaicatmod, pred = NeighboursYN, data = peaidata, main = "PEAI", x.label = "")
e <- cat_plot(pldecatmod, pred = NeighboursYN, data = pldedata, main = "PLDE", x.label = "")
f <- cat_plot(polecatmod, pred = NeighboursYN, data = poledata, main = "POLE", x.label = "")
g <- cat_plot(trcycatmod, pred = NeighboursYN, data = trcydata, main = "TRCY")
h <- cat_plot(trorcatmod, pred = NeighboursYN, data = trordata, main = "TROR")
i <- cat_plot(verocatmod, pred = NeighboursYN, data = verodata, main = "VERO")
grid.arrange(a, b, c, d, e, f, g, h, i, nrow = 3, ncol =3)
dev.off()

### Trying survival another way, gives the same plots as cat_plot
library(ggeffects)
arca <- ggpredict(arcacatmod, "NeighboursYN")
hygl <- ggpredict(hyglcatmod, "NeighboursYN")
laro <- ggpredict(larocatmod, "NeighboursYN")
peai <- ggpredict(peaicatmod, "NeighboursYN")
plde <- ggpredict(pldecatmod, "NeighboursYN")
pole <- ggpredict(polecatmod, "NeighboursYN")
trcy <- ggpredict(trcycatmod, "NeighboursYN")
tror <- ggpredict(trorcatmod, "NeighboursYN")
vero <- ggpredict(verocatmod, "NeighboursYN")
grid.arrange(arca, hygl, laro, peai, plde, pole, trcy, tror, vero, nrow = 3, ncol =3)

arcaplot <- plot(arca)
hyglplot <- plot(hygl)
laroplot <- plot(laro)
peaiplot <- plot(peai)
pldeplot <- plot(plde)
poleplot <- plot(pole)
trcyplot <- plot(trcy)
trorplot <- plot(tror)
veroplot <- plot(vero)

grid.arrange(arcaplot, hyglplot, laroplot, peaiplot, pldeplot, poleplot, trcyplot, trorplot, veroplot)

### Making panels of each species' responses to neighbours and PC1 ####
#### NA ####
#8x3, 8 species, survival, fecundity, lambda responses
#Need to say in caption that neighbour abundance is logged plus 1, standardised

#frame = FALSE is a nice way to remove box around plot
#can use log="y" to make y axis logged but need to plus one and for model curve to match

#population growth rate is log plus 1
#neighbour abundance is log plus 1

dev.off()
pdf("Output/Figures/panel_NA.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))
#square doesn't work, too small
#par(pty="s")
#ARCA
#Survival - linear model won't converge with Site included
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|plotid), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
stripchart(log_lambda_p1 ~ neighbours01, lambdaarca, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex = 2, add = TRUE)
text(x = 1.5, y = 3.3, "*", cex = 4, col = "red")

#Adding HYGL
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
stripchart(log_lambda_p1 ~ neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 4, col = "red")

#Adding LARO
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedlaro)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedlaro)
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
stripchart(log_lambda_p1 ~ neighbours01, lambdalaro, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 1.85, "*", cex = 4, col = "red")

#Adding PEAI
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedpeai)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedpeai)
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
stripchart(log_lambda_p1 ~ neighbours01, lambdapeai, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)

#Adding PLDE
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, pldedata)
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedplde)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedplde)
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 85,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
stripchart(log_lambda_p1 ~ neighbours01, lambdaplde, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.45, "*", cex = 4, col = "red")

#Adding TRCY
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
stripchart(log_lambda_p1 ~ neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 4, col = "red")

#Adding TROR
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, trordata)
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), log = "y", xlim=c(-1,2.5), ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtror)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtror)
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
stripchart(log_lambda_p1 ~ neighbours01, lambdatror, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.5, "*", cex = 4, col = "red")

#Adding VERO
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), log = "y", xlim=c(-1,2.5), ylim=c(1,220), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedvero)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedvero)
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 200,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab="", names = c("No neighbours", "Neighbours"), col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
stripchart(log_lambda_p1 ~ neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.3, "*", cex = 4, col = "red")

###Overall text
##x labels
mtext("Neighbour abundance", adj = 0.14, side = 1, cex = 2.2,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, cex = 2.2, outer = TRUE)
mtext("Neighbour presence", adj = 0.92, side = 1, cex = 2.2, outer = TRUE)
##y labels
#mtext("Probability of survival", cex = 3, outer=TRUE, side = 2)
#Figure out how to put these in the middle - actually not sure I want to
#mtext("Number of seeds produced", cex = 3, outer=TRUE, side = 2)
#mtext("Population growth rate (log plus 1)", cex = 3, outer=TRUE, side = 2)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 3)
mtext("Fecundity", outer=TRUE, adj = 0.54, side = 3, cex = 3)
mtext("Lambda", outer=TRUE, adj=0.90, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", adj = -0.13, side = 3, cex = 3, outer = TRUE)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 12, side = 3, cex = 2, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 25, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 35, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 45, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 44, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 65, side = 3, cex = 2, outer = TRUE)
mtext(~italic("V. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

dev.off()

#### PC1 ####
#Adding in germination
#Lambda no longer boxplot, continuous

dev.off()
pdf("Output/Figures/panel_PC1.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))
#ARCA
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 4)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdaarca)
x_to_plot<-seq.func(lambdaarca$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding HYGL
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity - this one is better not logged but logging for the sake of consistency!
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdahygl)
x_to_plot<-seq.func(lambdahygl$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding LARO
#Germination - quadratic
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedlaro)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedlaro)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 57,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdalaro)
x_to_plot<-seq.func(lambdalaro$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding PEAI
#Germination - quadratic I(std_PC1^2)
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1,800), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedpeai)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedpeai)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 490,"*", cex = 4, col = "red")
#Lambda - quadratic
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
model<-lmer(log_lambda_p1~std_PC1 + I(std_PC1^2) + (1|Site/Plot), lambdapeai)
x_to_plot<-seq.func(lambdapeai$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding PLDE
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, pldedata)
x_to_plot<-seq.func(pldedata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival - quadratic
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
model<-glmer(surv_to_produce_seeds~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, pldedata)
x_to_plot<-seq.func(pldedata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedplde)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedplde)
x_to_plot<-seq.func(pldedata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 85,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdaplde)
x_to_plot<-seq.func(lambdaplde$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding TRCY
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdatrcy)
x_to_plot<-seq.func(lambdatrcy$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding TROR
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, trordata)
x_to_plot<-seq.func(trordata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, trordata)
x_to_plot<-seq.func(trordata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtror)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtror)
x_to_plot<-seq.func(trordata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 55,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdatror)
x_to_plot<-seq.func(lambdatror$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding VERO
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity - quadratic
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedvero)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, seedvero)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 4, col = "red")

###Overall text
##x labels
mtext("std_PC1", adj = 0.11, side = 1, cex = 2.2,outer = TRUE)
mtext("std_PC1", adj = 0.4, side = 1, cex = 2.2,outer = TRUE)
mtext("std_PC1", adj = 0.65, side = 1, cex = 2.2, outer = TRUE)
mtext("std_PC1", adj = 0.92, side = 1, cex = 2.2, outer = TRUE)
##y labels
#mtext("Probability of survival", cex = 3, outer=TRUE, side = 2)
#Figure out how to put these in the middle - actually not sure I want to
#mtext("Number of seeds produced", cex = 3, outer=TRUE, side = 2)
#mtext("Population growth rate (log plus 1)", cex = 3, outer=TRUE, side = 2)
##main labels
mtext("Germination", outer=TRUE, adj=0.09,side = 3, cex = 3)
mtext("Survival", outer=TRUE, adj=0.4,side = 3, cex = 3)
mtext("Fecundity", outer=TRUE, adj = 0.65, side = 3, cex = 3)
mtext("Lambda", outer=TRUE, adj=0.93, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.13, side = 3, cex = 3, )
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 12, side = 3, cex = 2, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 25, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 35, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 45, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 44, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 65, side = 3, cex = 2, outer = TRUE)
mtext(~italic("V. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

dev.off()

########### Organising by similarity of responses and key examples PC1 ######

dev.off()
pdf("Output/Figures/panel_PC1_reorganised.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))

#Adding LARO
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedlaro)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedlaro)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 57,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdalaro)
x_to_plot<-seq.func(lambdalaro$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding TROR
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, trordata)
x_to_plot<-seq.func(trordata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, trordata)
x_to_plot<-seq.func(trordata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtror)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtror)
x_to_plot<-seq.func(trordata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 55,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdatror)
x_to_plot<-seq.func(lambdatror$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding PLDE
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, pldedata)
x_to_plot<-seq.func(pldedata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival - quadratic
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
model<-glmer(surv_to_produce_seeds~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, pldedata)
x_to_plot<-seq.func(pldedata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedplde)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedplde)
x_to_plot<-seq.func(pldedata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 85,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdaplde)
x_to_plot<-seq.func(lambdaplde$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding PEAI
#Germination - quadratic I(std_PC1^2)
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1,800), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedpeai)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedpeai)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 490,"*", cex = 4, col = "red")
#Lambda - quadratic
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
model<-lmer(log_lambda_p1~std_PC1 + I(std_PC1^2) + (1|Site/Plot), lambdapeai)
x_to_plot<-seq.func(lambdapeai$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#ARCA
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 4)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdaarca)
x_to_plot<-seq.func(lambdaarca$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding HYGL
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity - this one is better not logged but logging for the sake of consistency!
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdahygl)
x_to_plot<-seq.func(lambdahygl$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding TRCY
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdatrcy)
x_to_plot<-seq.func(lambdatrcy$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding VERO
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity - quadratic
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedvero)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, seedvero)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 4, col = "red")

###Overall text
##x labels
mtext("std_PC1", adj = 0.11, side = 1, cex = 2.2,outer = TRUE)
mtext("std_PC1", adj = 0.4, side = 1, cex = 2.2,outer = TRUE)
mtext("std_PC1", adj = 0.65, side = 1, cex = 2.2, outer = TRUE)
mtext("std_PC1", adj = 0.92, side = 1, cex = 2.2, outer = TRUE)
##y labels
#mtext("Probability of survival", cex = 3, outer=TRUE, side = 2)
#Figure out how to put these in the middle - actually not sure I want to
#mtext("Number of seeds produced", cex = 3, outer=TRUE, side = 2)
#mtext("Population growth rate (log plus 1)", cex = 3, outer=TRUE, side = 2)
##main labels
mtext("Germination", outer=TRUE, adj=0.09,side = 3, cex = 3)
mtext("Survival", outer=TRUE, adj=0.4,side = 3, cex = 3)
mtext("Fecundity", outer=TRUE, adj = 0.65, side = 3, cex = 3)
mtext("Lambda", outer=TRUE, adj=0.93, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.13, side = 3, cex = 3)
mtext(~italic("L. rosea"), adj = -0.15, padj= 5, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 15, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 25, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 35, side = 3, cex = 2, outer = TRUE)
mtext(~italic("A. calendula"), adj = -0.15, padj= 45, side = 3, cex = 2, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 44, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 51, side = 3, cex = 2, outer = TRUE)
mtext(~italic("V. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

dev.off()

### Pulling out four examples
#Arca - neutral
#Trcy - one positive but neutral lambda
#Vero - one positive and positive lambda
#Laro - counteracting

dev.off()
pdf("Output/Figures/panel_PC1_examples.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))
#ARCA
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 4)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdaarca)
x_to_plot<-seq.func(lambdaarca$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding TRCY
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdatrcy)
x_to_plot<-seq.func(lambdatrcy$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding VERO
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity - quadratic
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedvero)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, seedvero)
x_to_plot<-seq.func(verodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 4, col = "red")

#Adding LARO
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 4, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedlaro)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedlaro)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 57,"*", cex = 4, col = "red")
#Lambda
plot(log_lambda_p1 ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
model<-lmer(log_lambda_p1~std_PC1 + (1|Site/Plot), lambdalaro)
x_to_plot<-seq.func(lambdalaro$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

###Overall text
##x labels
mtext("std_PC1", adj = 0.11, side = 1, padj = -37.5, cex = 2.2,outer = TRUE)
mtext("std_PC1", adj = 0.4, side = 1, padj = -37.5, cex = 2.2,outer = TRUE)
mtext("std_PC1", adj = 0.65, side = 1, padj = -37.5, cex = 2.2, outer = TRUE)
mtext("std_PC1", adj = 0.92, side = 1, padj = -37.5, cex = 2.2, outer = TRUE)
##main labels
mtext("Germination", outer=TRUE, adj=0.09,side = 3, cex = 3)
mtext("Survival", outer=TRUE, adj=0.4,side = 3, cex = 3)
mtext("Fecundity", outer=TRUE, adj = 0.65, side = 3, cex = 3)
mtext("Lambda", outer=TRUE, adj=0.93, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.13, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 13, side = 3, cex = 2, outer = TRUE)
mtext(~italic("V. rosea"), adj = -0.15, padj= 25, side = 3, cex = 2, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 35, side = 3, cex = 2, outer = TRUE)

dev.off()


########### Organising by similarity of responses and key examples NA ######

dev.off()
pdf("Output/Figures/panel_NA_reorganised.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))

#Adding HYGL
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
stripchart(log_lambda_p1 ~ neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 4, col = "red")

#Adding LARO
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedlaro)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedlaro)
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
stripchart(log_lambda_p1 ~ neighbours01, lambdalaro, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 1.85, "*", cex = 4, col = "red")

#Adding PEAI
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedpeai)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedpeai)
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
stripchart(log_lambda_p1 ~ neighbours01, lambdapeai, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.1, "*", cex = 4, col = "red")

#Adding TROR
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, trordata)
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), log = "y", xlim=c(-1,2.5), ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtror)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtror)
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
stripchart(log_lambda_p1 ~ neighbours01, lambdatror, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.5, "*", cex = 4, col = "red")

#Adding PLDE
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, pldedata)
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedplde)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedplde)
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 85,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
stripchart(log_lambda_p1 ~ neighbours01, lambdaplde, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.45, "*", cex = 4, col = "red")


#Adding VERO
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), log = "y", xlim=c(-1,2.5), ylim=c(1,220), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedvero)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedvero)
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 200,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab="", names = c("No neighbours", "Neighbours"), col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
stripchart(log_lambda_p1 ~ neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.3, "*", cex = 4, col = "red")

#ARCA
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
stripchart(log_lambda_p1 ~ neighbours01, lambdaarca, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex = 2, add = TRUE)
text(x = 1.5, y = 3.3, "*", cex = 4, col = "red")

#Adding TRCY
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
stripchart(log_lambda_p1 ~ neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 4, col = "red")

###Overall text
##x labels
mtext("Neighbour abundance", adj = 0.14, side = 1, cex = 2.2,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, cex = 2.2, outer = TRUE)
mtext("Neighbour presence", adj = 0.92, side = 1, cex = 2.2, outer = TRUE)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 3)
mtext("Fecundity", outer=TRUE, adj = 0.54, side = 3, cex = 3)
mtext("Lambda", outer=TRUE, adj=0.90, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", adj = -0.13, side = 3, cex = 3, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 5, side = 3, cex = 2, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 14, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 25, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 35, side = 3, cex = 2, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 45, side = 3, cex = 2, outer = TRUE)
mtext(~italic("V. rosea"), adj = -0.15, padj= 55, side = 3, cex = 2, outer = TRUE)
mtext(~italic("A. calendula"), adj = -0.15, padj= 65, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 60, side = 3, cex = 2, outer = TRUE)

dev.off()

####### Pulling out four examples ### 
#How to choose which?? Most data?
#Hygl - non-sig, sig lambda
#Vero - Sig fecundity,  sig lambda
#Trcy - quad survival, sig lambda

dev.off()
pdf("Output/Figures/panel_NA_examples.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))

#Adding HYGL
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
stripchart(log_lambda_p1 ~ neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 4, col = "red")

#Adding VERO
#Survival
plot(surv_to_produce_seeds ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|Site/Plot), family = binomial, verodata)
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), log = "y", xlim=c(-1,2.5), ylim=c(1,220), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedvero)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedvero)
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 200,"*", cex = 4, col = "red")
#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names = NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
stripchart(log_lambda_p1 ~ neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.3, "*", cex = 4, col = "red")

#Adding TRCY
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 4, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
##Testing with full model instead here** 
#Why is it making me put 0s for random effects too?? but not above? check this* *stuck*
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                              I(std_PC1^2) + (1|Site/Plot),family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, 0, 0, 0, x_to_plot, 0, 0^2, 0))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Lambda
boxplot(log_lambda_p1 ~ neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names = c("No neighbours", "Neighbours"), col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
stripchart(log_lambda_p1 ~ neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 4, col = "red")

###Overall text
##x labels
mtext("Neighbour abundance", adj = 0.12, side = 1, padj = -47, cex = 2.2,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, padj = -47, cex = 2.2, outer = TRUE)
mtext("Neighbour presence", adj = 0.92, side = 1, padj = -47, cex = 2.2, outer = TRUE)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 3)
mtext("Fecundity", outer=TRUE, adj = 0.54, side = 3, cex = 3)
mtext("Lambda", outer=TRUE, adj=0.90, side = 3, cex = 3)
mtext("Species", adj = -0.13, padj = 0.1, side = 3, cex = 3, outer = TRUE)
#Use mxtext for species names
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 4, side = 3, cex = 2, outer = TRUE)
mtext(~italic("V. rosea"), adj = -0.15, padj= 15, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 20, side = 3, cex = 2, outer = TRUE)

dev.off()

#### Creating a table of final model output ####
#Need to run models first, run all of 'modelling vital rates in response to all factors'

## GERMINATION ##
## Extracting values for all in a loop
germ_model_list <- list(arcagermfinalmod, hyglgermfinalmod, larogermfinalmod, peaigermfinalmod, pldegermfinalmod, polegermfinalmod, trcygermfinalmod, trorgermfinalmod, verogermfinalmod)
effects = lapply(1:length(germ_model_list), function(x) {
  as.data.frame(coef(summary(germ_model_list[[x]]))) %>% mutate(Species=paste0(x))})
germ_effects_table <- do.call("rbind", effects)

#Make rownames a column 
germ_effects_table <- cbind(Effect = rownames(germ_effects_table), germ_effects_table)
#Remove rownames
rownames(germ_effects_table) <- NULL

#Renaming effects since loop adding values to ends
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, '(Intercept)')] <- 'Intercept'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'std_PC1')] <- 'std_PC1'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'std_PC2')] <- 'std_PC2'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'I(std_PC1^2)')] <- 'std_PC1^2'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'I(std_PC2^2)')] <- 'std_PC2^2'

germ_effects_table <- within(germ_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_effects_table <- within(germ_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_effects_table <- within(germ_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_effects_table <- within(germ_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_effects_table <- within(germ_effects_table, Species[Species == '5'] <- 'Plantago debilis')
germ_effects_table <- within(germ_effects_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_effects_table <- within(germ_effects_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_effects_table <- within(germ_effects_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_effects_table <- within(germ_effects_table, Species[Species == '9'] <- 'Velleia rosea')

#Renaming columns
germ_effects_table <- germ_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')

#Making column with Estimate (+/- SE) and p value asterisks all combined
#germ_effects_table$collated <- sprintf("%1.1f ± %1.1f", germ_effects_table$Estimate, germ_effects_table$SE)

#Add column for asterisks based on below function
germ_effects_table <- germ_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                              p_value <0.001~"***",
                                                              p_value <0.01~"**",
                                                              p_value <0.05~"*"))
germ_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", germ_effects_table$Estimate, germ_effects_table$SE, germ_effects_table$p_asterisks)


germ_effects_kbl <- germ_effects_table %>% select(Species, Effect, collated)
#test <- germ_effects_table %>% pivot_wider(names_from = Species, values_from = collated)
#test2 <- germ_effects_table %>% pivot_wider(names_from = Species, values_from = c(Estimate, SE, p_value, collated))
# test3 <- germ_effects_table %>% group_by(Species) %>%  
#   mutate(row = row_number()) %>%  
#   pivot_wider(names_from = Species, values_from = c(Estimate, SE, p_value, collated)) %>%
#   select(-row)
germ_effects_kbl <- germ_effects_kbl %>% group_by(Species) %>% mutate(row = row_number()) %>%
  pivot_wider(names_from = Species, values_from = collated) %>% select(-row)

#view(germ_effects_kbl)

#germ_effects_kbl[3,5] <- cell_spec(germ_effects_kbl[3,5], "latex", bold = T)

#Plotting with kableR
germ_effects_kbl %>% mutate(Effect = c("Intercept", "PC1", "PC2", "PC2^2", "PC1^2")) %>%
kbl(align = 'lccccccccc', caption = "<b>Supplementary X</b>. Germination effects table with Estimate ± SE. Asterisks denote significance: * p<0.05, ** p<0.01, *** p<0.001") %>%
    kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
    add_header_above(c("Germination" = 1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
    row_spec(0, italic = T) %>%
    column_spec(1, italic = F) %>%
  column_spec(1:10, width = 4)
  #column_spec(2, bold = ifelse(ends_with("*") <0.05, TRUE, FALSE))
  
#Figure out how to make ^2 superscript and
# how to bold rows containing asterisks

## SURVIVAL ##  
surv_model_list <- list(arcasurvfinalmod, hyglsurvfinalmod, larosurvfinalmod, peaisurvfinalmod, pldesurvfinalmod, trcysurvfinalmod, trorsurvfinalmod, verosurvfinalmod)
effects = lapply(1:length(surv_model_list), function(x) {
  as.data.frame(coef(summary(surv_model_list[[x]]))) %>% mutate(Species=paste0(x))})
surv_effects_table <- do.call("rbind", effects)

surv_effects_table <- cbind(Effect = rownames(surv_effects_table), surv_effects_table)
rownames(surv_effects_table) <- NULL

surv_effects_table$Effect[startsWith(surv_effects_table$Effect, '(Intercept)')] <- 'Intercept'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_PC1:std_logp1_totalabund')] <- 'PC1:Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry:std_PC1')] <- 'Dry:PC1'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet:std_PC1')] <- 'Wet:PC1'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry:I(std_logp1_totalabund^2)')] <- 'Dry:Neighbour abundance^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet:I(std_logp1_totalabund^2)')] <- 'Wet:Neighbour abundance^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry:std_logp1_totalabund')] <- 'Dry:Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet:std_logp1_totalabund')] <- 'Wet:Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'Dodder01')] <- 'Dodder'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'I(std_PC1^2)')] <- 'PC1^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'I(std_PC2^2)')] <- 'PC2^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'I(std_logp1_totalabund^2)')] <- 'Neighbour abundance^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_logp1_totalabund')] <- 'Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_PC1')] <- 'PC1'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_PC2')] <- 'PC2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet')] <- 'Wet'

surv_effects_table <- within(surv_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
surv_effects_table <- within(surv_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
surv_effects_table <- within(surv_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
surv_effects_table <- within(surv_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
surv_effects_table <- within(surv_effects_table, Species[Species == '5'] <- 'Plantago debilis')
surv_effects_table <- within(surv_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
surv_effects_table <- within(surv_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
surv_effects_table <- within(surv_effects_table, Species[Species == '8'] <- 'Velleia rosea')

surv_effects_table <- surv_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')
surv_effects_table <- surv_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
surv_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", surv_effects_table$Estimate, surv_effects_table$SE, surv_effects_table$p_asterisks)

surv_effects_kbl <- surv_effects_table %>% select(Species, Effect, collated)

#Not sure why Dry and Wet are duplicating, extensively troubleshooted and not sure, but grouping by effects works well
surv_effects_kbl <- surv_effects_kbl %>% group_by(Species, Effect) %>% mutate(row = row_number()) %>%
  pivot_wider(names_from = Species, values_from = collated) %>% select(-row)

#Plotting with kableR
surv_effects_kbl %>% kbl(align = 'lcccccccc', caption = "<b>Supplementary X</b>. Survival effects table with Estimate ± SE. Asterisks denote significance: * p<0.05, ** p<0.01, *** p<0.001") %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("Survival" = 1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

## FECUNDITY ##  
seed_model_list <- list(arcaseedfinalmod, hyglseedfinalmod, laroseedfinalmod, peaiseedfinalmod, pldeseedfinalmod, trcyseedfinalmod, trorseedfinalmod, veroseedfinalmod)
#glmmTMB model coefs need to be extracted slightly differently
effects = lapply(1:length(seed_model_list), function(x) {
  as.data.frame(coef(summary(seed_model_list[[x]]))[["cond"]]) %>% mutate(Species=paste0(x))})
seed_effects_table <- do.call("rbind", effects)

seed_effects_table <- cbind(Effect = rownames(seed_effects_table), seed_effects_table)
rownames(seed_effects_table) <- NULL

seed_effects_table$Effect[startsWith(seed_effects_table$Effect, '(Intercept)')] <- 'Intercept'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_PC1:std_logp1_totalabund')] <- 'PC1:Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry:std_PC1')] <- 'Dry:PC1'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet:std_PC1')] <- 'Wet:PC1'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry:I(std_logp1_totalabund^2)')] <- 'Dry:Neighbour abundance^2'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet:I(std_logp1_totalabund^2)')] <- 'Wet:Neighbour abundance^2'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry:std_logp1_totalabund')] <- 'Dry:Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet:std_logp1_totalabund')] <- 'Wet:Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'Dodder01')] <- 'Dodder'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'I(std_PC1^2)')] <- 'PC1^2'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'I(std_PC2^2)')] <- 'PC2^2'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'I(std_logp1_totalabund^2)')] <- 'Neighbour abundance^2'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_logp1_totalabund')] <- 'Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_PC1')] <- 'PC1'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_PC2')] <- 'PC2'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet')] <- 'Wet'

seed_effects_table <- within(seed_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
seed_effects_table <- within(seed_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
seed_effects_table <- within(seed_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
seed_effects_table <- within(seed_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
seed_effects_table <- within(seed_effects_table, Species[Species == '5'] <- 'Plantago debilis')
seed_effects_table <- within(seed_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
seed_effects_table <- within(seed_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
seed_effects_table <- within(seed_effects_table, Species[Species == '8'] <- 'Velleia rosea')

seed_effects_table <- seed_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')
seed_effects_table <- seed_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
seed_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", seed_effects_table$Estimate, seed_effects_table$SE, seed_effects_table$p_asterisks)

seed_effects_kbl <- seed_effects_table %>% select(Species, Effect, collated)

seed_effects_kbl <- seed_effects_kbl %>% group_by(Species) %>% mutate(row = row_number()) %>%
  pivot_wider(names_from = Species, values_from = collated) %>% select(-row)

#Plotting with kableR
seed_effects_kbl %>% kbl(align = 'lcccccccc', caption = "<b>Supplementary X</b>. Fecundity effects table with Estimate ± SE. Asterisks denote significance: * p<0.05, ** p<0.01, *** p<0.001") %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("Fecundity" = 1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

## LAMBDA ##  
lambda_model_list <- list(arcalambdafinalmod, hygllambdafinalmod, larolambdafinalmod, peailambdafinalmod, pldelambdafinalmod, trcylambdafinalmod, trorlambdafinalmod, verolambdafinalmod)
#glmmTMB model coefs need to be extracted slightly differently
effects = lapply(1:length(lambda_model_list), function(x) {
  as.data.frame(coef(summary(lambda_model_list[[x]]))) %>% mutate(Species=paste0(x))})
lambda_effects_table <- do.call("rbind", effects)

lambda_effects_table <- cbind(Effect = rownames(lambda_effects_table), lambda_effects_table)
rownames(lambda_effects_table) <- NULL

lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, '(Intercept)')] <- 'Intercept'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_PC1:neighbours01neighbours')] <- 'PC1:Neighbour presence'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry:I(std_PC1^2)')] <- 'Dry:PC1^2'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet:I(std_PC1^2)')] <- 'Wet:PC1^2'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'I(std_PC1^2)')] <- 'PC1^2'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'neighbours01neighbours')] <- 'Neighbour presence'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_PC1')] <- 'PC1'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_PC2')] <- 'PC2'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet')] <- 'Wet'

lambda_effects_table <- within(lambda_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '5'] <- 'Plantago debilis')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '8'] <- 'Velleia rosea')

lambda_effects_table <- lambda_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|t|)')
lambda_effects_table <- lambda_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                                p_value <0.001~"***",
                                                                                p_value <0.01~"**",
                                                                                p_value <0.05~"*"))
lambda_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", lambda_effects_table$Estimate, lambda_effects_table$SE, lambda_effects_table$p_asterisks)

lambda_effects_kbl <- lambda_effects_table %>% select(Species, Effect, collated)

lambda_effects_kbl <- lambda_effects_kbl %>% group_by(Species) %>% mutate(row = row_number()) %>%
  pivot_wider(names_from = Species, values_from = collated) %>% select(-row)

#Plotting with kableR

lambda_effects_kbl %>% kbl(align = 'lcccccccc', caption = "<b>Supplementary X</b>. Population growth rate effects table with Estimate ± SE. Asterisks denote significance: * p<0.05, ** p<0.01, *** p<0.001") %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("Lambda" = 1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1, "n = "=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)
  
######## Creating stacked bar chart with counts of effect responses ##########
# Produce one graph per demographic rate - 4 total
## Germination ##
germtally <- germ_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                   pos = sum(p_value<0.05 & Estimate>0),
                                                                   neg = sum(p_value<0.05 & Estimate<0))
#Manually adjusting for quadratic terms
#checking significance of linear terms where quadratic terms are ns
# ggplot(arcadata, aes(y = percent_germ, x = std_PC1))+
#   geom_jitter(alpha=0.2)+
#   geom_smooth(formula = y ~ poly(x, 2))+
#   theme_classic()

#convex - positive quadratic, high left and right, low middle
#concave - negative
#PC1: 1 pos_quad (-1 pos peai). PC2: 2 neg_quad (-1 pos plde, -1 pos laro)
germtally <- within(germtally, pos[Effect == 'std_PC1'] <- '5')
germtally$convex <- 0
germtally <- within(germtally, convex[Effect == 'std_PC1'] <- '2')
germtally$concave <- 0
germtally <- within(germtally, concave[Effect == 'std_PC2'] <- '2')
germtally <- within(germtally, ns[Effect == 'std_PC2'] <- '6')
germtally <- within(germtally, convex[Effect == 'std_PC2'] <- '1')
germtally <- within(germtally, pos[Effect == 'std_PC2'] <- '0')

germtally <- subset(germtally, Effect == "std_PC1" | Effect == "std_PC2")

germtally$pos <- as.numeric(germtally$pos)
germtally$convex <- as.numeric(germtally$convex)
germtally$concave <- as.numeric(germtally$concave)
germtally$ns <- as.numeric(germtally$ns)
#Pivotting longer for plotting
germtally_long <- germtally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom? 
#NS, pos, neg, pos_quad, neg_quad
germtally_long$response_type <- factor(germtally_long$response_type, level = c("concave", "convex", "neg", "pos", "ns"))
## Plotting as proportional bar chart
a <- ggplot(germtally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 1))+
  ylab("Count of effect responses")+
  scale_x_discrete(labels = c("PC1", "PC2"))+
  theme_classic()+
  scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "convex" = "grey20", "concave" = "grey0"))+
  ggtitle("(a)")

## Survival ##
# Need to fix my problem where I only included quadratic terms in a couple of interactions
#Dry:NA and Wet:NA need to be fixed ** do this**
survtally <- surv_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                   pos = sum(p_value<0.05 & Estimate>0),
                                                                   neg = sum(p_value<0.05 & Estimate<0))
#PC1^2 (one plde, but ns and PC1 ns too) and PC2^2 (none) good as is 
#ARCA has ns quadratic but significant linear response to neighbour abundance
#NA: 1 pos_quad (-1 neg trcy). 

#plotting to make sure this linear significant relationship is real, yes content that it is
# ggplot(arcadata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
#   geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
#   geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ x)+
#   theme_classic()
# testmod <- glm(surv_to_produce_seeds ~ std_logp1_totalabund, family = "binomial", arcadata)
# summary(testmod)

survtally <- within(survtally, neg[Effect == 'Neighbour abundance'] <- '1')
survtally$pos_quad <- 0
survtally <- within(survtally, pos_quad[Effect == 'Neighbour abundance'] <- '1')
survtally <- filter(survtally, !Effect == "Neighbour abundance^2", !Effect == "Dry:Neighbour abundance^2", !Effect == "Wet:Neighbour abundance^2", !Effect == "Intercept", !Effect == "PC1^2")

survtally$neg <- as.numeric(survtally$neg)
survtally$pos_quad <- as.numeric(survtally$pos_quad)
#Pivotting longer for plotting
survtally_long <- survtally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom
survtally_long$response_type <- factor(survtally_long$response_type, level = c("neg_quad", "pos_quad", "neg", "pos", "ns"))

#Reorder the effects?
#c("PC1", "PC2", "Neighbour abundance", "Dry", "Wet", "Dodder", "PC1:Neighbour abundance", .....)

## Plotting as proportional bar chart
b <- ggplot(survtally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 1))+
  ylab("Count of effect responses")+
 scale_x_discrete(labels = c("Dodder", "Dry", "Dry:N.A.", "Dry:PC1", "N.A.", "PC1", "PC1:N.A.", "PC2", "Wet", "Wet:N.A.", "Wet:PC1"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "neg" = "grey40", "pos_quad" = "grey20"))+
  ggtitle("(b)")

## Fecundity ##
seedtally <- seed_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                   pos = sum(p_value<0.05 & Estimate>0),
                                                                   neg = sum(p_value<0.05 & Estimate<0))
#PC1: vero and trcy ns PC1 and PC^2
seedtally <- filter(seedtally, !Effect == "Intercept", !Effect == "PC1^2")
#Pivotting longer for plotting
seedtally_long <- seedtally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom
seedtally_long$response_type <- factor(seedtally_long$response_type, level = c("neg_quad", "pos_quad", "neg", "pos", "ns"))

#Reorder the effects?
#c("PC1", "PC2", "Neighbour abundance", "Dry", "Wet", "Dodder", "PC1:Neighbour abundance", .....)

## Plotting as proportional bar chart
c <- ggplot(seedtally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 1))+
  ylab("Count of effect responses")+
 scale_x_discrete(labels = c("Dodder", "Dry", "Dry:N.A.", "Dry:PC1", "N.A.", "PC1", "PC1:N.A.", "PC2", "Wet", "Wet:N.A.", "Wet:PC1"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "neg" = "grey40"))+
  ggtitle("(c)")

## Lambda ##
# Need to fix my problem where I only included quadratic terms in a couple of interactions
#Dry:PC1^2 peai and Dry:PC1^2 ** do this**

lambdatally <- lambda_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                       pos = sum(p_value<0.05 & Estimate>0),
                                                                       neg = sum(p_value<0.05 & Estimate<0))
#PEAI has ns PC1^2 but significant linear response to PC1
#plotting to make sure this linear significant relationship is real, yes content that it is
# ggplot(lambdapeai, aes(x = std_PC1, y = lambda))+
#   geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
#   geom_smooth(method = "lm")+
#   theme_classic()
# testmod <- lm(lambda ~ std_PC1, lambdapeai)
# summary(testmod)
lambdatally <- filter(lambdatally, !Effect == "Dry:PC1^2", !Effect == "Wet:PC1^2", !Effect == "Intercept", !Effect == "PC1^2")
#Pivotting longer for plotting
lambdatally_long <- lambdatally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom
lambdatally_long$response_type <- factor(lambdatally_long$response_type, level = c("neg", "pos", "ns"))

#Reorder the effects?
#c("PC1", "PC2", "Neighbour abundance", "Dry", "Wet", "Dodder", "PC1:Neighbour abundance", .....)

## Plotting as proportional bar chart
d <- ggplot(lambdatally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 1))+
  ylab("Count of effect responses")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "neg" = "grey40"))+
  ggtitle("(d)")

grid.arrange(a, b, c, d)

####### Making a map of Perenjori! ####
library(ggplot2)
library(dplyr)
library(sf)
#Read in the SA2 shapefile downloaded from the ABS
#Data from ABS localities
# https://data.gov.au/dataset/ds-wa-d2dc22c6-0840-448c-819f-b6fb21411517/details?q=
#Following this guide: https://medium.com/analytics-vidhya/mapping-australia-in-r-6ce092c48b49

ausplotdata <- read_sf("Data_and_code_from_others/SA2_2016_AUST.shp")
head(ausplotdata)

#filter the Australian SA2 shapefile for only WA
waplotdata <- ausplotdata %>% filter(STE_NAME16 == "Western Australia")

ggplot()+
  geom_sf(data = ausplotdata)+
  xlab("Longitude")+
  ylab("Latitude") +
  xlim(110,155)+
  theme_classic()

ggplot()+
  geom_sf(data = ausplotdata)+
  geom_sf(data = waplotdata, fill = "blue") +
  xlab("Longitude")+
  ylab("Latitude") +
  xlim(110,155)+
  theme_classic()

ggplot()+
  geom_sf(data = waplotdata) +
  xlab("Longitude")+
  ylab("Latitude") +
  theme_classic()

#import a shapefile of state boundaries
aus_state_data <- read_sf("Data_and_code_from_others/STE_2016_AUST.shp")

ggplot()+
  geom_sf(data = aus_state_data)+
  theme_classic()

#make a new dataset of cities in Australia (google the locations)
#
#"West Perenjori Nature Reserve", -29.46703, 116.20600
wa_cities <- tribble(
  ~city, ~lat, ~long, 
  "Perth",-31.953512, 115.857048,
  "Perenjori", -29.443172, 116.288301)


#convert those two columns to geometry column with the st_as_sf() function. Google Maps uses the coordinate reference system 4326 (the GPS system).
wa_cities_geometry <- wa_cities %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

ggplot() +
  geom_sf(data=aus_state_data)+
  geom_sf(data=wa_cities_geometry, size=1)+
  geom_sf_label(data =wa_cities_geometry, aes(label = city))+
  xlim(110,155)+
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#Using ggrepel package to try and offset city labels and give them points
library(ggrepel)
dev.off()
pdf("Output/Figures/Perenjori_map.pdf")
ggplot() + geom_sf(data = aus_state_data, fill = "white") + 
  geom_text_repel(data= wa_cities,aes(x=long, y=lat, label=city), fontface = "bold",nudge_x = c(-3,3), nudge_y = c(-3,3)) +
  geom_point(data = wa_cities, aes(x = long, y = lat), size = 3) +  
  xlim(110,155)+
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()
dev.off()

state_labels <- tribble(
  ~state, ~lat, ~long, 
  "Western Australia",-26.519050, 122.165609,
  "Northern Territory", -19.567005, 133.587342,
  "Victoria", -36.999911, 144.079037,
  "New South Wales", -32.979835, 146.336442,
  "Australian Capital Territory", -35.430785, 148.959320,
  "Queensland", -22.880978, 144.214005,
  "Tasmania", -41.900339, 146.462846,
  "South Australia", -30.070128, 134.832533)
state_labels_geometry <- state_labels %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

#with state labels
ggplot() + geom_sf(data = aus_state_data, fill = "white") + 
  geom_text_repel(data= wa_cities,aes(x=long, y=lat, label=city), fontface = "bold",nudge_x = c(-3,3), nudge_y = c(-3,3)) +
  geom_point(data = wa_cities, aes(x = long, y = lat), size = 3)+
  geom_sf_label(data = state_labels_geometry, aes(label = state), size = 2)+
  xlim(110,155)+
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#making dataframe for my 4 open and 4 shaded blocks (colour coded)

block_coords <- tribble(
  ~block, ~lat, ~long, 
  "Open 1", -29.4791851, 116.1986948,
  "Open 2", -29.4794714,	116.1987961,
  "Open 3", -29.4791268,	116.1982687,
  "Open 4", -29.479611,	116.197514,
  "Shade 1", -29.4792905,	116.1984353,
  "Shade 2", -29.4790544,	116.1982767,
  "Shade 3", -29.4795126,	116.19695,
  "Shade 4", -29.4798626,	116.1968421)



