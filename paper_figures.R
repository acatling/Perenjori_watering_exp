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
library(lmerTest)
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

## LARO - quadratic PC2 germination
## PLDE - quadratic PC2 germination
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

### ARCA - quadratic neighbour abundance survival
### TRCY - quadratic neighbour abundance survival
## All other species linear
#HYGL isn't converging now?! ugh says it has 100 observations but filter(Species == "HYGL" & !is.na(Total_abundance)) gives 102 observations...
ggplot(hygldata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
  geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
  geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ poly(x, 2))+
  theme_classic()
#It's not even remotely quadratic, could that be the problem? Do not include quadratic term
hyglsurvnamodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), 
                         family = binomial, hygldata)
summary(hyglsurvnamodel)

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
arcaspquad1 <- glmmTMB(No_viable_seeds_grouped ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = nbinom2, seedarca)
summary(arcaspquad1)
coef(summary(arcaspquad1))$cond[3,4]
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

### Making a table with these significance values ####
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
vitalquadpvalues %>% kbl(caption = "<b>Supplementary 1</b>. Model output <i>p</i>-values for quadratic terms of fixed effects. N.A. denotes neighbour abundance. <b>Bolded</b> values indicate a <i>p</i>-value of <0.05.", digits = c(5, 2, 2, 2, 3, 2, 2, 2, 2, 2)) %>%
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
## ARCA
arcagermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                        family = binomial, arcadata)
arcagermmod1dharma <- simulateResiduals(arcagermmod1)
plot(arcagermmod1dharma)
#Not too bad
vif(arcagermmod1)
summary(arcagermmod1)
arcagermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = binomial, arcadata)
summary(arcagermfinalmod)

## HYGL
hyglgermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = binomial, hygldata)
hyglgermmod1dharma <- simulateResiduals(hyglgermmod1)
plot(hyglgermmod1dharma)
#Not great residuals
vif(hyglgermmod1)
summary(hyglgermmod1)
hyglgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, hygldata)

## LARO - has quadratic germ ~ PC2
larogermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 +
                        I(std_PC2^2) + (1|Site/Plot), 
                      family = binomial, larodata)
larogermmod1dharma <- simulateResiduals(larogermmod1)
plot(larogermmod1dharma)
#not too bad
vif(larogermmod1)
summary(larogermmod1)
larogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 +
                            I(std_PC2^2) + (1|Site/Plot), 
                          family = binomial, larodata)

## PEAI - has quadratic germ ~ PC1
peaigermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + 
                        I(std_PC1^2) + (1|Site/Plot), 
                      family = binomial, peaidata)
peaigermmod1dharma <- simulateResiduals(peaigermmod1)
plot(peaigermmod1dharma)
#Not too bad residuals, dispersion test and KS test significant though
vif(peaigermmod1)
summary(peaigermmod1)
peaigermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + 
                            I(std_PC1^2) + (1|Site/Plot), 
                          family = binomial, peaidata)

## PLDE - has quadratic germ ~ PC2
pldegermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + 
                        I(std_PC2^2) + (1|Site/Plot), 
                      family = binomial, pldedata)
pldegermmod1dharma <- simulateResiduals(pldegermmod1)
plot(pldegermmod1dharma)
#Not too bad residuals, KS, dispersion and outlier tests significant though
vif(pldegermmod1)
summary(pldegermmod1)
pldegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + 
                            I(std_PC2^2)  + (1|Site/Plot), 
                          family = binomial, pldedata)

## POLE
polegermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = binomial, poledata)
polegermmod1dharma <- simulateResiduals(polegermmod1)
plot(polegermmod1dharma)
#good
vif(polegermmod1)
summary(polegermmod1)
polegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, poledata)

## TRCY
trcygermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = binomial, trcydata)
trcygermmod1dharma <- simulateResiduals(trcygermmod1)
plot(trcygermmod1dharma)
#Not too bad residuals, KS, dispersion and outlier tests significant though
vif(trcygermmod1)
summary(trcygermmod1)
trcygermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, trcydata)

## TROR
trorgermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = binomial, trordata)
trorgermmod1dharma <- simulateResiduals(trorgermmod1)
plot(trorgermmod1dharma)
#Not great residuals, KS test significant
vif(trorgermmod1)
summary(trorgermmod1)
trorgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, trordata)

## VERO
verogermmod1 <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                      family = binomial, verodata)
verogermmod1dharma <- simulateResiduals(verogermmod1)
plot(verogermmod1dharma)
#Not great residuals, KS and dispersion tests significant
vif(verogermmod1)
summary(verogermmod1)
verogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, verodata)

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

## HYGL - needs optimiser to  converge
hyglsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
hyglsurvmod1dharma <- simulateResiduals(hyglsurvmod1)
plot(hyglsurvmod1dharma)
#good
vif(hyglsurvmod1)
summary(hyglsurvmod1)
# Model simplification step - Significant Dry:PC1, Wet:PC1, Wet:NA, PC1:NA, keeping all
hyglsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
summary(hyglsurvfinalmod)

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
larosimplemod <- glm(surv_to_produce_seeds ~ std_logp1_totalabund, family = "binomial", larodata)
summary(larosimplemod)
#almost but not quite
larosimplemod2 <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + (1|Site/Plot), family = "binomial", larodata)
summary(larosimplemod2)
#random effects barely make a difference

## PEAI 
peaisurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), peaidata)
peaisurvmod1dharma <- simulateResiduals(peaisurvmod1)
plot(peaisurvmod1dharma)
#good
vif(peaisurvmod1)
summary(peaisurvmod1)
# Model simplification step - Nothing significant, removing all
peaisurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
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
#good
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
#Problem here*, with estimates
polesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), poledata)
polesurvmod1dharma <- simulateResiduals(polesurvmod1)
plot(polesurvmod1dharma)
#not too bad
vif(polesurvmod1)
summary(polesurvmod1)
# Model simplification step - Nothing significant, removing all
# polesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
#                           family = binomial, poledata)
# polesurvfinalmoddharma <- simulateResiduals(polesurvfinalmod)
# plot(polesurvfinalmoddharma)
# #good
# summary(polesurvfinalmod)

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
trcysurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            I(std_logp1_totalabund^2) + Treatment:std_PC1 + Treatment:I(std_logp1_totalabund^2) + (1|Site/Plot), 
                          family = binomial, trcydata)
trcysurvfinalmoddharma <- simulateResiduals(trcysurvfinalmod)
plot(trcysurvfinalmoddharma)
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
# Model simplification step - Nothing significant, removing all
peaiseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedpeai)
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
poleseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                          Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedpole)
# Model simplification step - Removing all interactions, seeing if simple model will fit
poleseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedpole)
poleseedfinalmoddharma <- simulateResiduals(poleseedfinalmod)
plot(poleseedfinalmoddharma)
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
                              I(std_PC1^2) + (1|Site/Plot), 
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
#Model simplification step - No significant interactions, removing all interactions
larolambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + (1|Site/Plot), lambdalaro)
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
#Model simplification step - Significant Wet:PC1^2, removing all other interactions
peailambdafinalmod <- lmer(log_lambda_p1 ~ Treatment + std_PC1 + std_PC2 + neighbours01 + 
                             I(std_PC1^2) + Treatment:I(std_PC1^2) + (1|Site/Plot), lambdapeai)
peailambdafinalmoddharma <- simulateResiduals(peailambdafinalmod)
plot(peailambdafinalmoddharma)
#good
summary(peailambdafinalmod)

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

###########################################################
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
