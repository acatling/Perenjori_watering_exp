##### Created 21/04/22
# This script will contain the code needed to generate figures and tables in the paper
# WA Perenjori Watering Experiment
# Alexandra Catling

#### Loading packages and data ####
#Data imported from data preparation sheet
source("data_preparation.R")
#vitaldata has all datasets combined 
# germination, survival, seed production, neighbour info, quad factors
# one row per subplot with seeds sown - NAs are very important since, e.g. survival info is only on germinated subplots

#Functions file
source("functions.R")

#Packages
library(ggplot2)
library(ggrepel)
library(MuMIn)
library(DHARMa)
library(glmmTMB)
library(kableExtra)
library(grid)
library(car)
library(sjPlot)
library(gridExtra)
library(corrplot)

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
  plot(plotted.data$percent_germ~plotted.data$std_PC2, pch=19, col="grey60", ylab="Probability of germination", xlab="PC2 (standardised)", cex.lab=2, .axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13)
  title(main=bquote(italic(.(species.name.list[i]))), .main=2.5)
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

## All species - linear PC2 survival, (LARO = 0.053)
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
## All other species linear, (ARCA = 0.052)
# #ARCA linear model not converging, will converge with 1|plotid
# ggplot(arcadata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
#   geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
#   geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ poly(x, 2))+
#   theme_classic()
# arcanamodel <- glmer(surv_to_produce_seeds ~ std_logp1_totalabund + (1|Site/Plot),
#                         family = binomial, arcadata)
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
  lambdaquadPC1 <- lmer(log_lambda ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
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
  plot(plotted.data$log_lambda~plotted.data$std_PC1, pch=19, col="grey60", ylab="Population growth rate (logged plus 1)", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC1) 
  model<-lmer(log_lambda~std_PC1 + I(std_PC1^2) + (1|Site/Plot), plotted.data)
  model2<-lmer(log_lambda~std_PC1 + (1|Site/Plot), plotted.data)
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
plot(log_lambda ~ std_PC1, lambdapeai)
model <- lmer(log_lambda~std_PC1 + I(std_PC1^2) + (1|Site/Plot), lambdapeai)
x_to_plot<-seq.func(lambdapeai$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

## All species - linear PC2 lambda
for (i in 1:length(specieslist)){
  print(specieslist[i])
  lambdaquadPC2 <- lmer(log_lambda ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
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
  plot(plotted.data$log_lambda~plotted.data$std_PC2, pch=19, col="grey60", ylab="Population growth rate (logged plus 1)", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_PC2) 
  model<-lmer(log_lambda~std_PC2 + I(std_PC2^2) + (1|Site/Plot), plotted.data)
  model2<-lmer(log_lambda~std_PC2 + (1|Site/Plot), plotted.data)
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
rownames(germquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Goodenia rosea")
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
colnames(survquadtable) <- c('PC1', 'PC2', 'N. Ab.')
rownames(survquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Goodenia rosea")
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
survquadtable[6,1] <- NA
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
survquadtable[6,2] <- NA
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
survquadtable[6,3] <- NA
survquadtable[7,3] <- survNApvalueTRCY
survquadtable[8,3] <- survNApvalueTROR
survquadtable[9,3] <- survNApvalueVERO

## Fecundity
#Create table to put these values in
fecundityquadtable <- matrix(ncol=3, nrow = 9)
#define column names and row names of matrix
colnames(fecundityquadtable) <- c('PC1', 'PC2', 'N. Ab.')
rownames(fecundityquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Goodenia rosea")
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
fecundityquadtable[6,1] <- NA
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
fecundityquadtable[6,2] <- NA
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
fecundityquadtable[6,3] <- NA
fecundityquadtable[7,3] <- seedNApvalueTRCY
fecundityquadtable[8,3] <- seedNApvalueTROR
fecundityquadtable[9,3] <- seedNApvalueVERO

## Lambda
#Create table to put these values in
lambdaquadtable <- matrix(ncol=2, nrow = 9)
#define column names and row names of matrix
colnames(lambdaquadtable) <- c('PC1', 'PC2')
rownames(lambdaquadtable) <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Goodenia rosea")
#convert matrix to table 
lambdaquadtable <- as.data.frame(lambdaquadtable)

## Extracting values for PC1
for (i in 1:length(specieslist)){
  model <- lmer(log_lambda ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), 
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
lambdaquadtable[6,1] <- NA
lambdaquadtable[7,1] <- lambdaPC1pvalueTRCY
lambdaquadtable[8,1] <- lambdaPC1pvalueTROR
lambdaquadtable[9,1] <- lambdaPC1pvalueVERO
## Extracting values for PC2
for (i in 1:length(specieslist)){
  model <- lmer(log_lambda ~ std_PC2 + I(std_PC2^2) + (1|Site/Plot), 
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
lambdaquadtable[6,2] <- NA
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
  kbl(caption = "<b>Supplementary Information 3</b>. Model output <i>p</i>-values for quadratic terms of fixed effects. 
      N. Ab. denotes neighbour abundance and NA signifies insufficient data to run models. <b>Bolded</b> values indicate a <i>p</i>-value of <0.05.", digits = c(4, 3, 2, 2, 3, 2, 2, 2, 2, 2)) %>%
  kable_classic(full_width = F, html_font = "Times") %>%
  column_spec(1, italic = T) %>%
  #row_spec(3, bold = ifelse(vitalquadpvalues[1,] <0.05, TRUE, FALSE)) %>%
  column_spec(2, bold = ifelse(vitalquadpvalues[[1]] >0.05 | is.na(vitalquadpvalues[[1]]), FALSE, TRUE)) %>%
  column_spec(3, bold = ifelse(vitalquadpvalues[[2]] >0.05 | is.na(vitalquadpvalues[[2]]), FALSE, TRUE)) %>%
  column_spec(4, bold = ifelse(vitalquadpvalues[[3]] >0.05 | is.na(vitalquadpvalues[[3]]), FALSE, TRUE)) %>%
  column_spec(5, bold = ifelse(vitalquadpvalues[[4]] >0.05 | is.na(vitalquadpvalues[[4]]), FALSE, TRUE)) %>%
  column_spec(6, bold = ifelse(vitalquadpvalues[[5]] >0.05 | is.na(vitalquadpvalues[[5]]), FALSE, TRUE)) %>%
  column_spec(7, bold = ifelse(vitalquadpvalues[[6]] >0.05 | is.na(vitalquadpvalues[[6]]), FALSE, TRUE)) %>%
  column_spec(8, bold = ifelse(vitalquadpvalues[[7]] >0.05 | is.na(vitalquadpvalues[[7]]), FALSE, TRUE)) %>%
  column_spec(9, bold = ifelse(vitalquadpvalues[[8]] >0.05 | is.na(vitalquadpvalues[[8]]), FALSE, TRUE)) %>%
  column_spec(10, bold = ifelse(vitalquadpvalues[[9]] >0.05 | is.na(vitalquadpvalues[[9]]), FALSE, TRUE)) %>%
  add_header_above(c("", "Emergence" = 2, "Survival" = 3, "Fecundity" = 3, "Population growth" = 2))
  #add_header_above(c("", "p value quadratic term" = 8))
#Struggling to save this using save_kable and as_image() atm.
#Can copy it from the Viewer using copy to clipboard, maintain aspect ratio, first value 500

#### Modelling vital rates in response to all factors #### 
#### Germination ####

### Individual level random effects can help with overdispersion
# Try adding row ID for each subplot

## ARCA - has quadratic germ ~ PC2
arcagermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC2^2) + (1|Site/Plot/rowID), 
                          family = binomial, arcadata)
arcagermmod1dharma <- simulateResiduals(arcagermfinalmod)
plot(arcagermmod1dharma)
vif(arcagermfinalmod)
summary(arcagermfinalmod)
testDispersion(arcagermfinalmod)

## HYGL- has quadratic germ ~ PC2
hyglgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC2^2) + (1|Site/Plot/rowID), 
                          family = binomial, hygldata)
hyglgermmod1dharma <- simulateResiduals(hyglgermfinalmod)
plot(hyglgermmod1dharma)
vif(hyglgermfinalmod)
summary(hyglgermfinalmod)
testDispersion(hyglgermfinalmod)

## LARO - has quadratic germ ~ PC1
larogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC1^2) + (1|Site/Plot/rowID), 
                          family = binomial, larodata)
larogermmod1dharma <- simulateResiduals(larogermfinalmod)
plot(larogermmod1dharma)
vif(larogermfinalmod)
summary(larogermfinalmod)
testDispersion(larogermfinalmod)

## PEAI - has quadratic germ ~ PC1
peaigermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC1^2) + (1|Site/Plot/rowID), 
                          family = binomial, peaidata)
peaigermmod1dharma <- simulateResiduals(peaigermfinalmod)
plot(peaigermmod1dharma)
vif(peaigermfinalmod)
summary(peaigermfinalmod)
testDispersion(peaigermfinalmod)

## PLDE - has quadratic germ ~ PC2
pldegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC2^2)  + (1|Site/Plot/rowID), 
                          family = binomial, pldedata)
pldegermmod1dharma <- simulateResiduals(pldegermfinalmod)
plot(pldegermmod1dharma)
vif(pldegermfinalmod)
summary(pldegermfinalmod)
testDispersion(pldegermfinalmod)

## POLE
polegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot/rowID), 
                          family = binomial, poledata)
polegermmod1dharma <- simulateResiduals(polegermfinalmod)
plot(polegermmod1dharma)
vif(polegermfinalmod)
summary(polegermfinalmod)
testDispersion(polegermfinalmod)

## TRCY
trcygermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot/rowID), 
                          family = binomial, trcydata)
trcygermmod1dharma <- simulateResiduals(trcygermfinalmod)
plot(trcygermmod1dharma)
vif(trcygermfinalmod)
summary(trcygermfinalmod)
testDispersion(trcygermfinalmod)

## TROR
trorgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot/rowID), 
                          family = binomial, trordata)
trorgermmod1dharma <- simulateResiduals(trorgermfinalmod)
plot(trorgermmod1dharma)
vif(trorgermfinalmod)
summary(trorgermfinalmod)
testDispersion(trorgermfinalmod)

## VERO
#won't converge with rowID and has terrible dispersion without it
verogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + (1|Site/Plot), 
                          family = binomial, verodata)
verogermmod1dharma <- simulateResiduals(verogermfinalmod)
plot(verogermmod1dharma)
#Not great residuals, KS test significant
vif(verogermfinalmod)
summary(verogermfinalmod)
testDispersion(verogermfinalmod)

#### Survival ####
## ARCA
arcasurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, arcadata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
arcasurvmod1dharma <- simulateResiduals(arcasurvmod1)
plot(arcasurvmod1dharma)
#good
vif(arcasurvmod1)
summary(arcasurvmod1)
# Model simplification step - Significant Dry:PC1 and PC1:NA, removing all others
arcasurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            Treatment:std_PC1 + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
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
# Model simplification step - Significant Dry:PC1 and Wet:PC1, removing others
hyglsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            Treatment:std_PC1 + (1|Site/Plot), 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
hyglsurvfinalmoddharma <- simulateResiduals(hyglsurvfinalmod)
plot(hyglsurvfinalmoddharma)
vif(hyglsurvfinalmod)
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
                            std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                          family = binomial, peaidata)
peaisurvfinalmoddharma <- simulateResiduals(peaisurvfinalmod)
plot(peaisurvfinalmoddharma)
#good
vif(peaisurvfinalmod)
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
#Huge NA vif!
summary(pldesurvmod1)
# Model simplification step - significant PC1:NA, removing others
pldesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                            I(std_PC1^2) + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                          family = binomial, pldedata)
pldesurvfinalmoddharma <- simulateResiduals(pldesurvfinalmod)
plot(pldesurvfinalmoddharma)
#good
vif(pldesurvfinalmod)
summary(pldesurvfinalmod)

## POLE
#problem here - not using this model anyway
polesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                      family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), poledata)
# polesurvmod1dharma <- simulateResiduals(polesurvmod1)
# plot(polesurvmod1dharma)
# #okay
# vif(polesurvmod1)
# summary(polesurvmod1)
# # Model simplification step - Nothing significant, removing all
# polesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
#                           family = binomial, poledata)
# polesurvfinalmoddharma <- simulateResiduals(polesurvfinalmod)
# plot(polesurvfinalmoddharma)
# #okay
# summary(polesurvfinalmod)

## TRCY - quadratic surv ~ NA
trcysurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + 
                        I(std_logp1_totalabund^2) + Treatment:I(std_logp1_totalabund^2) + std_PC1:I(std_logp1_totalabund^2) + (1|Site/Plot), 
                      family = binomial, trcydata)
trcysurvmod1dharma <- simulateResiduals(trcysurvmod1)
plot(trcysurvmod1dharma)
#not too bad - some high NA^2, makes sense
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
#good
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

#### Fecundity ####
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
#I know glmmTMB says it has a convergence issue but when I run the glmer.nb (no issues), they give
#almost exactly the same coefficients and p-values, so keeping glmmTMB
peaiseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
                              Treatment:std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedpeai)
# peaiseedfinalmod2 <- glmer.nb(No_viable_seeds_grouped ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + 
#                                Treatment:std_logp1_totalabund + (1|Site/Plot), seedpeai)
#coef(summary(peaiseedfinalmod))
peaiseedfinalmoddharma <- simulateResiduals(peaiseedfinalmod)
plot(peaiseedfinalmoddharma)
#good
summary(peaiseedfinalmod)

#Not that much data with neighbours and a lot of zeros - problem?
#PC1 and PC2 seem fine
# ggplot(seedpeai, aes(x = std_logp1_totalabund, y = No_viable_seeds_grouped))+
#   geom_point(alpha=0.3)+
#   theme_classic()

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
                             + (1|Site/Plot), family = nbinom2, seedtrcy)
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

#### Lambda ####
## ARCA
hist(lambdaplde$log_lambda)
arcalambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdamod1dharma <- simulateResiduals(arcalambdamod1)
plot(arcalambdamod1dharma)
#good
vif(arcalambdamod1)
summary(arcalambdamod1)
#Model simplification step - No significant interactions, removing all interactions
arcalambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdafinalmoddharma <- simulateResiduals(arcalambdafinalmod)
plot(arcalambdafinalmoddharma)
#good
summary(arcalambdafinalmod)

## HYGL
hygllambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdamod1dharma <- simulateResiduals(hygllambdamod1)
plot(hygllambdamod1dharma)
#okay
vif(hygllambdamod1)
summary(hygllambdamod1)
#Model simplification step - No interactions significant, removing all
hygllambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdafinalmoddharma <- simulateResiduals(hygllambdafinalmod)
plot(hygllambdafinalmoddharma)
#okay
summary(hygllambdafinalmod)

##LARO
larolambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdalaro)
larolambdamod1dharma <- simulateResiduals(larolambdamod1)
plot(larolambdamod1dharma)
#okay
vif(larolambdamod1)
summary(larolambdamod1)
#Model simplification step - Significant PC1:neighbours, removing others
larolambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdalaro)
larolambdafinalmoddharma <- simulateResiduals(larolambdafinalmod)
plot(larolambdafinalmoddharma)
#not great
summary(larolambdafinalmod)

##PEAI - has quadratic lambda ~ PC1
peailambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + 
                         I(std_PC1^2) + Treatment:I(std_PC1^2) + I(std_PC1^2):Neighbours01 + (1|Site/Plot), lambdapeai)
peailambdamod1dharma <- simulateResiduals(peailambdamod1)
plot(peailambdamod1dharma)
#okay, dispersion test significant
vif(peailambdamod1)
summary(peailambdamod1)
#Model simplification step - Significant Dry:PC1 and PC1:neighbours, removing all other interactions
peailambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + 
                             I(std_PC1^2) + Treatment:std_PC1 + std_PC1:Neighbours01 + (1|Site/Plot), lambdapeai)
peailambdafinalmoddharma <- simulateResiduals(peailambdafinalmod)
plot(peailambdafinalmoddharma)
#good
summary(peailambdafinalmod)

##PLDE
pldelambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdamod1dharma <- simulateResiduals(pldelambdamod1)
plot(pldelambdamod1dharma)
#good
vif(pldelambdamod1)
summary(pldelambdamod1)
#Model simplification step - No significant interactions, removing all interactions
pldelambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdafinalmoddharma <- simulateResiduals(pldelambdafinalmod)
plot(pldelambdafinalmoddharma)
#okay
summary(pldelambdafinalmod)

##POLE - not enough data
# polelambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
#                          Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdapole)
# polelambdamod1dharma <- simulateResiduals(polelambdamod1)
# plot(polelambdamod1dharma)
# #not good
# vif(polelambdamod1)
# summary(polelambdamod1)
# #Model simplification step - No significant interactions, removing all interactions
# polelambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + (1|Site/Plot), lambdapole)
# polelambdafinalmoddharma <- simulateResiduals(polelambdafinalmod)
# plot(polelambdafinalmoddharma)
# #not good - not enough data? 18 points. Residuals bad
# summary(polelambdafinalmod)

##TRCY
trcylambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdamod1dharma <- simulateResiduals(trcylambdamod1)
plot(trcylambdamod1dharma)
#good
vif(trcylambdamod1)
summary(trcylambdamod1)
#Model simplification step - Significant Treatment:neighbours and PC1:neighbours, removing all other interactions
trcylambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + 
                             Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdafinalmoddharma <- simulateResiduals(trcylambdafinalmod)
plot(trcylambdafinalmoddharma)
#good
summary(trcylambdafinalmod)

##TROR
trorlambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdatror)
trorlambdamod1dharma <- simulateResiduals(trorlambdamod1)
plot(trorlambdamod1dharma)
#good
vif(trorlambdamod1)
summary(trorlambdamod1)
#Model simplification step - No significant interactions, removing all interactions
trorlambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + (1|Site/Plot), lambdatror)
trorlambdafinalmoddharma <- simulateResiduals(trorlambdafinalmod)
plot(trorlambdafinalmoddharma)
#good
summary(trorlambdafinalmod)

##VERO
verolambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
                         Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdavero)
verolambdamod1dharma <- simulateResiduals(verolambdamod1)
plot(verolambdamod1dharma)
#good
vif(verolambdamod1)
summary(verolambdamod1)
#Model simplification step - Significant Treatment:PC1 and PC1:neighbours, removing all other interactions
verolambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + 
                             Treatment:std_PC1 + std_PC1:Neighbours01 + (1|Site/Plot), lambdavero)
verolambdafinalmoddharma <- simulateResiduals(verolambdafinalmod)
plot(verolambdafinalmoddharma)
#good
summary(verolambdafinalmod)

#### Plotting vital rate responses to continuous and categorical neighbours ####
ggplot(popdata, aes(x = Neighbours01, y = log_lambda))+
  geom_boxplot()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species, scales = "free")
#Fecundity response to categorical neighbours - individual level
#Need to fix neighbour NAs - 20 or so.
ggplot(data = seedmodeldata %>% filter(complete.cases(Neighbours01)), aes(x = Neighbours01, y = No_viable_seeds_grouped+1))+
  geom_boxplot()+
  scale_y_log10()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)

#Survival response to categorical neighbours - PLOT level
ggplot(popdata, aes(x = Neighbours01, y = plot_survival))+
  geom_boxplot()+
  geom_jitter(alpha = 0.2, width = 0.1)+
  theme_classic()+
  my_theme+
  facet_wrap(~Species)

### Making panels of each species' responses to neighbours and PC1 NOT FROM MODELS ####
#### Neighbours ####
#8x3, 8 species, survival, fecundity, lambda (categorical) responses
#frame = FALSE is a nice way to remove box around plot
#can use log="y" to make y axis logged but need to plus one and for model curve to match

#population growth rate is log
#neighbour abundance is log plus 1

dev.off()
pdf("Output/Figures/panel_NA_old.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(5, 20, 5, 1), mar =c(3.5,6,1,1))
#ARCA
#Survival - linear model won't converge with Site included, using plotid instead
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + (1|plotid), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 100,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate (log)", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
stripchart(log_lambda ~ Neighbours01, lambdaarca, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex = 2, add = TRUE)
text(x = 1.5, y = 2.6, "*", cex = 10, col = "red")

#Adding HYGL
#Survival
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
stripchart(log_lambda ~ Neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)

#Adding LARO
#Survival
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
stripchart(log_lambda ~ Neighbours01, lambdalaro, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)

## Testing significance to neighbours alone - ns
# test <- lmer(log_lambda ~ Neighbours01 + (1|Site/Plot), lambdalaro)
# summary(test)

#Adding PEAI
#Survival
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
stripchart(log_lambda ~ Neighbours01, lambdapeai, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)

#Adding PLDE
#Survival
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, pldedata)
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
text(x = 2.4,y = 70,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
stripchart(log_lambda ~ Neighbours01, lambdaplde, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.1, "*", cex = 10, col = "red")

#Adding TRCY
#Survival - quadratic
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
stripchart(log_lambda ~ Neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)

## Testing significance to neighbours alone - ns
# test <- lmer(log_lambda ~ Neighbours01 + (1|Site/Plot), lambdatrcy)
# summary(test)

#Adding TROR
#Survival
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trordata)
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
stripchart(log_lambda ~ Neighbours01, lambdatror, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.3, "*", cex = 6, col = "red")

#Adding VERO
#Survival
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, verodata)
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
text(x = 2.4,y = 175,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab="", names = c("Absent", "Present"), col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
stripchart(log_lambda ~ Neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)

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
mtext(~italic("G. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

dev.off()

#### PC1 ####
#Adding in germination
#Lambda no longer boxplot, continuous

dev.off()
pdf("Output/Figures/panel_PC1.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(10, 20, 5, 1), mar =c(3.5,10,1,1))
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
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdaarca)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")

# # summary(model)
# ggplot(hygldata, aes(x = std_PC1, y = surv_to_produce_seeds, colour = Treatment))+
#   geom_jitter(alpha = 0.8, width = 0.05, height = 0.05)+
#   geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ x)+
#   theme_classic()
# #This doesn't look significant, not sure how/why - interactions?
# with(hygldata, plot(jitter(surv_to_produce_seeds, amount = 0.05) ~ std_PC1))
# #dry
# x_to_plot<-seq.func(hygldata$std_PC1)
# curve(plogis(cbind(1, 1, 0, x, 0, 0, 0, 1*x, 0*x)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "red")
# #wet
# curve(plogis(cbind(1, 0, 1, x, 0, 0, 0, 0*x, 1*x)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "blue")
# #ambient
# curve(plogis(cbind(1, 0, 0, x, 0, 0, 0, 0*x, 0*x)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "green")
# 
# #Is hygl significant ~ PC1?
# testhyglmod <- glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
# summary(testhyglmod)
# #I think not, just signif int with watering

#Fecundity - this one is better not logged but logging for the sake of consistency!
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdahygl)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 57,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdalaro)
x_to_plot<-seq.func(lambdalaro$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding PEAI
#Germination - quadratic
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 490,"*", cex = 6, col = "red")
#Lambda - quadratic
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
model<-lmer(log_lambda~std_PC1 + I(std_PC1^2) + (1|Site/Plot), lambdapeai)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdaplde)
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
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdatrcy)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 55,"*", cex = 6, col = "red")

## Is this signif? Int with Treatment
with(seedtror, plot(jitter(surv_to_produce_seeds, amount = 0.05) ~ std_PC1))
#dry
x_to_plot<-seq.func(seedtror$std_PC1)
curve(plogis(cbind(1, 1, 0, x, 0, 0, 0, 1*x, 0*x)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "red")
#wet
curve(plogis(cbind(1, 0, 1, x, 0, 0, 0, 0*x, 1*x)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "blue")
#ambient
curve(plogis(cbind(1, 0, 0, x, 0, 0, 0, 0*x, 0*x)%*%fixef(hyglsurvfinalmod)), add = TRUE, col = "green")

#Is hygl significant ~ PC1?
testtrormod <- glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtror)
summary(testtrormod)

#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdatror)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 6, col = "red")

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
mtext(~italic("G. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

dev.off()

#### Testing lambda plot more model terms #####
#This works
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda~std_PC1 + Neighbours01 + (1|Site/Plot), lambdatrcy)
x_to_plot<-seq.func(lambdatrcy$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, 0))
#no neighbours
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "red", env.trans = 50, line.colour = "red", line.weight = 2, line.type = 1)
#yes neighbours
preddata <- with(model, data.frame(1, x_to_plot, 1))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "green", env.trans = 50, line.colour = "green", line.weight = 2, line.type = 1)

### Interactions PC1 and neighbours?
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda~std_PC1 + Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdatrcy)
x_to_plot<-seq.func(lambdatrcy$std_PC1)
#interaction and with neighbours
preddata <- with(model, data.frame(1, x_to_plot, 1, x_to_plot*1))
#working
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
# interaction and without neighbours
preddata <- with(model, data.frame(1, x_to_plot, 0, x_to_plot*0))
#working
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
### All terms of full model
# trcylambdafinalmod <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + 
#                              Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdatrcy)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 + 
              Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdatrcy)
x_to_plot<-seq.func(lambdatrcy$std_PC1)
#ambient watering, mean PC2, no nbhs
preddata <- with(model, data.frame(1, 0, 0, x_to_plot, 0, 0, 0*0, 0*0, x_to_plot*0))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "red", env.trans = 50, line.colour = "red", line.weight = 2, line.type = 1)
#ambient watering, with nbhs
preddata <- with(model, data.frame(1, 0, 0, x_to_plot, 0, 1, 0*1, 0*1, x_to_plot*1))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "green", env.trans = 50, line.colour = "green", line.weight = 2, line.type = 1)


### PC1 one plot lambda ####
#### Making the panel with lambda calculated once per plot (not separated by presence/absence neighbours)
### Not updated with significance***
dev.off()
pdf("Output/Figures/panel_PC1_one_plot_26092022.pdf", width=21, height=21)
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
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 6)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdaarca_plot)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity - this one is better not logged but logging for the sake of consistency!
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdahygl_plot)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 57,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdalaro_plot)
x_to_plot<-seq.func(lambdalaro$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

#Adding PEAI
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, peaidata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, peaidata)
x_to_plot<-seq.func(peaidata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 490,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdapeai_plot)
x_to_plot<-seq.func(lambdapeai$std_PC1)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 85,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdaplde_plot)
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
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdatrcy_plot)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 55,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdatror_plot)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero_plot)
model<-lmer(log_lambda~std_PC1 + (1|Site), lambdavero_plot)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)

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
mtext(~italic("G. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 57,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdalaro)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 55,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdatror)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 85,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdaplde)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 490,"*", cex = 6, col = "red")
#Lambda - quadratic
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
model<-lmer(log_lambda~std_PC1 + I(std_PC1^2) + (1|Site/Plot), lambdapeai)
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
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 6)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdaarca)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
#Survival
plot(surv_to_produce_seeds ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, hygldata)
model<-glmer(surv_to_produce_seeds~std_PC1 + (1|Site/Plot), family = binomial, hygldata)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity - this one is better not logged but logging for the sake of consistency!
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), ylim=c(1, 90), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedhygl)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedhygl)
x_to_plot<-seq.func(hygldata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdahygl)
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
text(x = 2.5,y = 0.9,"*", cex = 6)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdatrcy)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 6, col = "red")

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
mtext(~italic("G. rosea"), adj = -0.15, padj= 76, side = 3, cex = 2, outer = TRUE)

dev.off()

### Pulling out four examples ####
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
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 6)
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdaarca)
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
text(x = 2.5,y = 0.9,"*", cex = 6)
#Fecundity
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", ylim=c(1, 100), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_PC1 + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 75,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdatrcy)
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
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdavero)
x_to_plot<-seq.func(lambdavero$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 4.2,"*", cex = 6, col = "red")

#Adding LARO
#Germination
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Germination fraction", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, larodata)
model<-glmer(cbind(total_germ, total_no_germ)~std_PC1 + (1|Site/Plot), family = binomial, larodata)
x_to_plot<-seq.func(larodata$std_PC1)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 1.5,y = 0.9,"*", cex = 6, col = "red")
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
text(x = 1.5,y = 57,"*", cex = 6, col = "red")
#Lambda
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab="Population growth rate", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
model<-lmer(log_lambda~std_PC1 + (1|Site/Plot), lambdalaro)
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
mtext(~italic("G. rosea"), adj = -0.15, padj= 25, side = 3, cex = 2, outer = TRUE)
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
stripchart(log_lambda ~ Neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 6, col = "red")

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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdalaro)
stripchart(log_lambda ~ Neighbours01, lambdalaro, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 1.85, "*", cex = 6, col = "red")

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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdapeai)
stripchart(log_lambda ~ Neighbours01, lambdapeai, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.1, "*", cex = 6, col = "red")

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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatror)
stripchart(log_lambda ~ Neighbours01, lambdatror, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.5, "*", cex = 6, col = "red")

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
text(x = 2.4,y = 85,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaplde)
stripchart(log_lambda ~ Neighbours01, lambdaplde, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.45, "*", cex = 6, col = "red")


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
text(x = 2.4,y = 200,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab="", names = c("No neighbours", "Neighbours"), col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
stripchart(log_lambda ~ Neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.3, "*", cex = 6, col = "red")

#ARCA
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab = NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, arcadata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, arcadata)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedarca)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedarca)
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.4,y = 260,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdaarca)
stripchart(log_lambda ~ Neighbours01, lambdaarca, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex = 2, add = TRUE)
text(x = 1.5, y = 3.3, "*", cex = 6, col = "red")

#Adding TRCY
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
#Fecundity
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), log = "y", ylim=c(1, 90), pch=19, col=alpha("grey60", 0.3), ylab="Number of seeds produced", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, seedtrcy)
model<-glmmTMB(No_viable_seeds_grouped~std_logp1_totalabund + (1|Site/Plot), family = nbinom2, seedtrcy)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
stripchart(log_lambda ~ Neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 6, col = "red")

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
mtext(~italic("G. rosea"), adj = -0.15, padj= 55, side = 3, cex = 2, outer = TRUE)
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names= NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdahygl)
stripchart(log_lambda ~ Neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 6, col = "red")

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
text(x = 2.4,y = 200,"*", cex = 6, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names = NA, col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdavero)
stripchart(log_lambda ~ Neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 4.3, "*", cex = 6, col = "red")

#Adding TRCY
#Survival - quadratic
plot(surv_to_produce_seeds ~ jitter(std_logp1_totalabund, 15), xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="Probability of survival", xlab=NA, tck=-0.01, cex= 2, cex.lab = 1.5, cex.axis = 1.5, trcydata)
model<-glmer(surv_to_produce_seeds~std_logp1_totalabund + I(std_logp1_totalabund^2) + (1|Site/Plot), family = binomial, trcydata)
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
preddata <- with(model, data.frame(1, x_to_plot, x_to_plot^2))
plotted.pred <- glmm.predict(mod = model, newdat = preddata, se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey1", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
text(x = 2.5,y = 0.9,"*", cex = 6, col = "red")
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
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="Population growth rate", xlab=NA, names = c("No neighbours", "Neighbours"), col = "white", cex= 2, cex.lab = 1.5, cex.axis = 1.5, lambdatrcy)
stripchart(log_lambda ~ Neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2, add = TRUE)
text(x = 1.5, y = 2.95, "*", cex = 6, col = "red")

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
mtext(~italic("G. rosea"), adj = -0.15, padj= 15, side = 3, cex = 2, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.15, padj= 20, side = 3, cex = 2, outer = TRUE)

dev.off()

#### Creating a table of final model outputs ####
#Need to run models first, run all of 'modelling vital rates in response to all factors'
#### Extracting marginal and conditional r squared values for table ####
#Extracting values for theoretical marginal R squared
germ_model_list <- list(arcagermfinalmod, hyglgermfinalmod, larogermfinalmod, peaigermfinalmod, pldegermfinalmod, polegermfinalmod, trcygermfinalmod, trorgermfinalmod, verogermfinalmod)
rsquared = lapply(1:length(germ_model_list), function(x) {
  as.data.frame(r.squaredGLMM(germ_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
germ_rsquared_table <- do.call("rbind", rsquared)

#Rename species
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '9'] <- 'Goodenia rosea')

## GERMINATION ####
## Extracting values for all in a loop
germ_model_list <- list(arcagermfinalmod, hyglgermfinalmod, larogermfinalmod, peaigermfinalmod, pldegermfinalmod, polegermfinalmod, trcygermfinalmod, trorgermfinalmod, verogermfinalmod)
effects = lapply(1:length(germ_model_list), function(x) {
  as.data.frame(coef(summary(germ_model_list[[x]]))) %>% mutate(Species=paste0(x))})
germ_effects_table <- do.call("rbind", effects)

#make rownames a column
germ_effects_table <- cbind(Effect = rownames(germ_effects_table), germ_effects_table)
rownames(germ_effects_table) <- NULL

#Rename species
germ_effects_table <- within(germ_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_effects_table <- within(germ_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_effects_table <- within(germ_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_effects_table <- within(germ_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_effects_table <- within(germ_effects_table, Species[Species == '5'] <- 'Plantago debilis')
germ_effects_table <- within(germ_effects_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_effects_table <- within(germ_effects_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_effects_table <- within(germ_effects_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_effects_table <- within(germ_effects_table, Species[Species == '9'] <- 'Goodenia rosea')

#Renaming effects since loop adding values to ends
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, '(Intercept)')] <- 'Intercept'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'std_PC1')] <- 'PC1'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'std_PC2')] <- 'PC2'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'I(std_PC1^2)')] <- 'PC1^2'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'I(std_PC2^2)')] <- 'PC2^2'

germ_effects_table <- within(germ_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_effects_table <- within(germ_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_effects_table <- within(germ_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_effects_table <- within(germ_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_effects_table <- within(germ_effects_table, Species[Species == '5'] <- 'Plantago debilis')
germ_effects_table <- within(germ_effects_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_effects_table <- within(germ_effects_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_effects_table <- within(germ_effects_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_effects_table <- within(germ_effects_table, Species[Species == '9'] <- 'Goodenia rosea')

#Renaming columns
germ_effects_table <- germ_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')

#Making column with Estimate (+/- SE) and p value asterisks all combined
#Add column for asterisks based on below function
germ_effects_table <- germ_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                              p_value <0.001~"***",
                                                              p_value <0.01~"**",
                                                              p_value <0.05~"*"))
germ_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", germ_effects_table$Estimate, germ_effects_table$SE, germ_effects_table$p_asterisks)

#Join with r squared values
#germ_effects_table <- left_join(germ_effects_table, germ_rsquared_table)

germ_effects_kbl <- germ_effects_table %>% select(Species, Effect, collated)
germ_effects_kbl <- germ_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)

#Plotting with kableR
germ_effects_kbl %>% kbl(align = 'lccccccccc') %>%
    kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.16/0.82"=1, "0.34/0.92"=1, "0.29/0.84"=1, "0.29/0.89"=1, "0.25/0.91"=1, "0.21/0.96"=1, "0.09/0.92"=1, "0.30/0.90"=1, "0.51/0.85"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Emergence" = 1, "n=190"=1, "n=189"=1, "n=192"=1, "n=192"=1, "n=92"=1, "n=185"=1, "n=192"=1, "n=191"=1, "n=191"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
    column_spec(1, italic = F) %>%
  column_spec(1:10, width = 4)

## SURVIVAL #### 
### r squared values
surv_model_list <- list(arcasurvfinalmod, hyglsurvfinalmod, larosurvfinalmod, peaisurvfinalmod, pldesurvfinalmod, trcysurvfinalmod, trorsurvfinalmod, verosurvfinalmod)
rsquared = lapply(1:length(surv_model_list), function(x) {
  as.data.frame(r.squaredGLMM(surv_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
surv_rsquared_table <- do.call("rbind", rsquared)

#Rename species
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '7'] <- 'Trachymene ornata')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '8'] <- 'Goodenia rosea')


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
surv_effects_table <- within(surv_effects_table, Species[Species == '8'] <- 'Goodenia rosea')

surv_effects_table <- surv_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')
surv_effects_table <- surv_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
surv_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", surv_effects_table$Estimate, surv_effects_table$SE, surv_effects_table$p_asterisks)

surv_effects_kbl <- surv_effects_table %>% select(Species, Effect, collated)

#Not sure why Dry and Wet are duplicating, extensively troubleshooted and not sure, but grouping by effects works well
surv_effects_kbl <- surv_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)

#Plotting with kableR
surv_effects_kbl %>% kbl(align = 'lcccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.43/0.51"=1, "0.72/0.76"=1, "0.05/0.06"=1, "0.10/0.27"=1, "0.26/0.26"=1, "0.30/0.31"=1, "0.11/0.31"=1, "0.17/0.36"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Survival" = 1, "n=163"=1, "n=104"=1, "n=163"=1, "n=162"=1, "n=76"=1, "n=159"=1, "n=160"=1, "n=143"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

## FECUNDITY #### 
seed_model_list <- list(arcaseedfinalmod, hyglseedfinalmod, laroseedfinalmod, peaiseedfinalmod, pldeseedfinalmod, trcyseedfinalmod, trorseedfinalmod, veroseedfinalmod)
rsquared = lapply(1:length(seed_model_list), function(x) {
  as.data.frame(r.squaredGLMM(seed_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
seed_rsquared_table <- do.call("rbind", rsquared)

#Rename species
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '7'] <- 'Trachymene ornata')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '8'] <- 'Goodenia rosea')


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
seed_effects_table <- within(seed_effects_table, Species[Species == '8'] <- 'Goodenia rosea')

seed_effects_table <- seed_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')
seed_effects_table <- seed_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
seed_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", seed_effects_table$Estimate, seed_effects_table$SE, seed_effects_table$p_asterisks)

seed_effects_kbl <- seed_effects_table %>% select(Species, Effect, collated)

seed_effects_kbl <- seed_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)

#Plotting with kableR
seed_effects_kbl %>% kbl(align = 'lcccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.26/0.69"=1, "0.24/0.24"=1, "0.25/0.39"=1, "0.43/0.43"=1, "0.38/0.48"=1, "0.18/0.19"=1, "0.26/0.26"=1, "0.23/0.47"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Seed production" = 1, "n=55"=1, "n=41"=1, "n=84"=1, "n=79"=1, "n=38"=1, "n=115"=1, "n=82"=1, "n=96"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

## LAMBDA ####  
lambda_model_list <- list(arcalambdafinalmod, hygllambdafinalmod, larolambdafinalmod, peailambdafinalmod, pldelambdafinalmod, trcylambdafinalmod, trorlambdafinalmod, verolambdafinalmod)
rsquared = lapply(1:length(lambda_model_list), function(x) {
  as.data.frame(r.squaredGLMM(lambda_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
lambda_rsquared_table <- do.call("rbind", rsquared)

#Rename species
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '7'] <- 'Trachymene ornata')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '8'] <- 'Goodenia rosea')

lambda_model_list <- list(arcalambdafinalmod, hygllambdafinalmod, larolambdafinalmod, peailambdafinalmod, pldelambdafinalmod, trcylambdafinalmod, trorlambdafinalmod, verolambdafinalmod)
effects = lapply(1:length(lambda_model_list), function(x) {
  as.data.frame(coef(summary(lambda_model_list[[x]]))) %>% mutate(Species=paste0(x))})
lambda_effects_table <- do.call("rbind", effects)

lambda_effects_table <- cbind(Effect = rownames(lambda_effects_table), lambda_effects_table)
rownames(lambda_effects_table) <- NULL

lambda_effects_table <- within(lambda_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '5'] <- 'Plantago debilis')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '8'] <- 'Goodenia rosea')

lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, '(Intercept)')] <- 'Intercept'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_PC1:Neighbours01Neighbours1')] <- 'PC1:Neighbour presence'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'I(std_PC1^2)')] <- 'PC1^2'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'Neighbours01Neighbours1')] <- 'Neighbour presence'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_PC1')] <- 'PC1'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_PC2')] <- 'PC2'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry:std_PC1')] <- 'Dry:PC1'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet:std_PC1')] <- 'Wet:PC1'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry:Neighbours01Neighbours1')] <- 'Dry:Neighbour presence'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet:Neighbours01Neighbours1')] <- 'Wet:Neighbour presence'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet')] <- 'Wet'

lambda_effects_table <- lambda_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|t|)')
lambda_effects_table <- lambda_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                   p_value <0.001~"***",
                                                   p_value <0.01~"**",
                                                   p_value <0.05~"*"))
lambda_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", lambda_effects_table$Estimate, lambda_effects_table$SE, lambda_effects_table$p_asterisks)
lambda_effects_kbl <- lambda_effects_table %>% select(Species, Effect, collated)
lambda_effects_kbl <- lambda_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)
#Plotting with kableR
lambda_effects_kbl %>% kbl(align = 'lcccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.18/0.23"=1, "0.25/0.50"=1, "0.21/0.61"=1, "0.37/0.66"=1, "0.17/0.80"=1, "0.21/0.64"=1, "0.09/0.36"=1, "0.25/0.54"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Population growth" = 1, "n=48"=1, "n=48"=1, "n=48"=1, "n=48"=1, "n=24"=1, "n=48"=1, "n=48"=1, "n=48"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
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
germtally <- within(germtally, pos[Effect == 'PC1'] <- '5')
germtally$convex <- 0
germtally <- within(germtally, convex[Effect == 'PC1'] <- '2')
germtally$concave <- 0
germtally <- within(germtally, concave[Effect == 'PC2'] <- '2')
germtally <- within(germtally, ns[Effect == 'PC2'] <- '6')
germtally <- within(germtally, convex[Effect == 'PC2'] <- '1')
germtally <- within(germtally, pos[Effect == 'PC2'] <- '0')

germtally <- subset(germtally, Effect == "PC1" | Effect == "PC2")

germtally$pos <- as.numeric(germtally$pos)
germtally$convex <- as.numeric(germtally$convex)
germtally$concave <- as.numeric(germtally$concave)
germtally$ns <- as.numeric(germtally$ns)
#Pivotting longer for plotting
germtally_long <- germtally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at the top 
#NS, pos, neg, pos_quad, neg_quad
germtally_long$response_type <- factor(germtally_long$response_type, level = c("ns", "concave", "convex", "pos", "neg"))
## Plotting as proportional bar chart
#lots of spaces in xlab to try and match up sizes of different plots
a <- ggplot(germtally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 2))+
  ylab("Count of response type")+
  xlab("")+
  scale_x_discrete(labels = c("                       PC1", "                     PC2"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=12, angle = 90, hjust=0.95, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=12, face = "bold"))+
  theme(legend.position="none")+
  scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "convex" = "grey20", "concave" = "grey0"))+
  ggtitle("A) Emergence")

## Survival ##
survtally <- surv_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                   pos = sum(p_value<0.05 & Estimate>0),
                                                                   neg = sum(p_value<0.05 & Estimate<0))
# ggplot(arcadata, aes(x = std_logp1_totalabund, y = surv_to_produce_seeds))+
#   geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
#   geom_smooth(method = "glm", method.args=list(family="binomial"), formula = y ~ x)+
#   theme_classic()
# testmod <- glm(surv_to_produce_seeds ~ std_logp1_totalabund, family = "binomial", arcadata)
# summary(testmod)
# testmod <- glm(surv_to_produce_seeds ~ Treatment, family = "binomial", trordata)
# summary(testmod)

survtally <- within(survtally, neg[Effect == 'Neighbour abundance'] <- '1')
survtally$convex <- 0
survtally <- within(survtally, convex[Effect == 'Neighbour abundance'] <- '1')
survtally <- filter(survtally, !Effect == "Dodder", !Effect == "Neighbour abundance^2", !Effect == "Dry:Neighbour abundance^2", !Effect == "Wet:Neighbour abundance^2", !Effect == "Intercept", !Effect == "PC1^2")

survtally$neg <- as.numeric(survtally$neg)
survtally$convex <- as.numeric(survtally$convex)
#Pivotting longer for plotting
survtally_long <- survtally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom
survtally_long$response_type <- factor(survtally_long$response_type, level = c("ns", "convex", "pos", "neg"))

#Reorder the effects
survtally_long$Effect <- factor(survtally_long$Effect, level = c("PC1", "PC2", "Neighbour abundance", "Dry", "Wet", "PC1:Neighbour abundance", "Dry:PC1", "Wet:PC1"))

## Plotting as proportional bar chart
b <- ggplot(survtally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 2))+
  ylab("")+
  xlab("")+
  scale_x_discrete(labels = c("PC1", "PC2", "N. abundance", "Dry", "Wet", "PC1:N. abundance", "Dry:PC1", "Wet:PC1"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=12, angle = 90, hjust=0.95, vjust=0.5), axis.text.y=element_text(size=12))+
  theme(legend.position="none")+
  scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "neg" = "grey40", "convex" = "grey20"))+
  ggtitle("B) Survival")

## Fecundity ##
seedtally <- seed_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                   pos = sum(p_value<0.05 & Estimate>0),
                                                                   neg = sum(p_value<0.05 & Estimate<0))
seedtally <- filter(seedtally, !Effect == "Dry:Neighbour abundance", !Effect == "Wet:Neighbour abundance", !Effect == "Wet:PC1", !Effect == "Intercept", !Effect == "PC1^2")
#Pivotting longer for plotting
seedtally_long <- seedtally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom
seedtally_long$response_type <- factor(seedtally_long$response_type, level = c("ns", "pos", "neg"))
#Reorder the effects
seedtally_long$Effect <- factor(seedtally_long$Effect, level = c("PC1", "PC2", "Neighbour abundance", "Dry", "Wet", "Dodder", "PC1:Neighbour abundance", "Dry:PC1", "Wet:PC1"))

## Plotting as proportional bar chart
c <- ggplot(seedtally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 2))+
  ylab("Count of response type")+
  xlab("Effect")+
  scale_x_discrete(labels = c("PC1", "PC2", "N. abundance", "Dry", "Wet", "Dodder", " PC1:N. abundance", "Dry:PC1", "Wet:PC1"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=12, angle = 90, hjust=0.95, vjust=0.5), axis.text.y=element_text(size=12), axis.title=element_text(size=12, face = "bold"))+
  theme(legend.position="none")+
  scale_fill_manual(values = c("ns" = "grey80", "pos" = "grey60", "neg" = "grey40"))+
  ggtitle("C) Seed production")

## Lambda ##
lambdatally <- lambda_effects_table %>% group_by(Effect) %>% summarise(ns = sum(p_value>0.05),
                                                                       pos = sum(p_value<0.05 & Estimate>0),
                                                                       neg = sum(p_value<0.05 & Estimate<0))
#plotting to make sure this linear significant relationship is real, yes content that it is
# ggplot(lambdapeai, aes(x = std_PC1, y = lambda))+
#   geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
#   geom_smooth(method = "lm")+
#   theme_classic()
# testmod <- lm(lambda ~ Treatment, lambdavero)
# summary(testmod)
lambdatally <- filter(lambdatally, !Effect == "Intercept", !Effect == "PC1^2", !Effect == "Wet:Neighbour presence")
#Pivotting longer for plotting
lambdatally_long <- lambdatally %>% pivot_longer(!Effect, names_to = "response_type", values_to = "count")
#Reordering response types so that non-significant is at bottom
lambdatally_long$response_type <- factor(lambdatally_long$response_type, level = c("ns", "pos", "neg"))
#Reorder the effects
lambdatally_long$Effect <- factor(lambdatally_long$Effect, level = c("PC1", "PC2", "Neighbour presence", "Dry", "Wet", "PC1:Neighbour presence", "Dry:Neighbour presence", "Dry:PC1", "Wet:PC1"))

## Plotting as proportional bar chart
#Create graph
d <- ggplot(lambdatally_long, aes(x = Effect, y = count))+
  geom_bar(aes(fill = response_type), position = "stack", stat = "identity")+
  scale_y_continuous(breaks = seq(0,10, by = 2))+
  ylab("")+
  xlab("Effect")+
  scale_x_discrete(labels = c("PC1", "PC2", "N. presence", "Dry", "Wet", "PC1:N. presence", "Dry:N. presence", "Dry:PC1", "Wet:PC1"))+
  theme_classic()+
  labs(fill="response type")+
  theme(axis.text.x = element_text(size=14, angle = 90, hjust=0.95, vjust=0.5), 
        legend.position = "top", axis.text.y=element_text(size=12), legend.text = element_text(size=12), legend.title = element_text(size=12), axis.title=element_text(size=12, face = "bold"))+
scale_fill_manual(values = c("ns" = "grey80", "positive" = "grey60", "negative" = "grey40", "convex" = "grey20", "concave" = "grey0"))+
  ggtitle("D) Population growth")
#Function needed to extract legend
#Save legend
legend <- get_legend(d)
#Removing legend from plot
d <- d + theme(legend.position="none")

#Heights needs to correspond to number of rows, widths to number of columns
# grid.arrange(legend, a, b, c, d, ncol=2, nrow = 3, 
#              layout_matrix = rbind(c(1,1), c(2,3), c(4,5)),
#              widths = c(2.7, 2.7), heights = c(0.2, 2.5, 2.5))

#Now with universal x and y labels
#a <- a + ylab("")
#d <- d + xlab("")
#c <- c + ylab("") + xlab("")

#y_title <- expression(paste(bold("Count of response type               Count of response type")))
#x_title <- expression(paste(bold("Effect                                                    Effect")))

grid.arrange(legend, a, b, c, d, ncol=2, nrow = 3, 
             layout_matrix = rbind(c(1,1), c(2,3), c(4,5)),
             widths = c(1, 1), heights = c(0.15, 1, 1))
             #left=textGrob(y_title, rot=90, gp=gpar(fontsize=14)),
            # bottom=textGrob(x_title, gp=gpar(fontsize=14)))
####### Making a map of Perenjori! ####
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
  geom_sf(data = aus_state_data, fill = "white")+
  theme_classic()

#make a new dataset of cities in Australia (google the locations)
#
#"West Perenjori Nature Reserve", -29.46703, 116.20600
wa_cities <- tribble(
  ~city, ~lat, ~long, 
  "Perth",-31.953512, 115.857048,
  "Perenjori", -29.443172, 116.288301)

wa_cities2 <- tribble(
  ~city, ~lat, ~long, 
  "Perenjori", -29.443172, 116.288301)

#convert those two columns to geometry column with the st_as_sf() function. Google Maps uses the coordinate reference system 4326 (the GPS system).
wa_cities_geometry <- wa_cities %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

wa_cities_geometry2 <- wa_cities2 %>% 
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
  #geom_text_repel(data= wa_cities,aes(x=long, y=lat, label=city), fontface = "bold",nudge_x = c(-3,3), nudge_y = c(-3,3)) +
  #geom_point(data = wa_cities, aes(x = long, y = lat), size = 3)+
  geom_sf_label(data = state_labels_geometry, aes(label = state), size = 2)+
  xlim(110,155)+
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#making dataframe for my 4 open and 4 shaded blocks (colour coded)

block_coords <- tribble(
  ~block, ~lat, ~long, ~openvshade,
  "Open 1", -29.4791851, 116.1986948, "open",
  "Open 2", -29.4794714,	116.1987961,"open",
  "Open 3", -29.4791268,	116.1982687,"open",
  "Open 4", -29.479611,	116.197514,"open",
  "Shade 1", -29.4792905,	116.1984353, "shaded",
  "Shade 2", -29.4790544,	116.1982767,"shaded",
  "Shade 3", -29.4795126,	116.19695,"shaded",
  "Shade 4", -29.4798626,	116.1968421, "shaded")

## Cutting this to just WA and Perenjori label
ggplot() + geom_sf(data = aus_state_data, fill = "white") + 
  geom_text_repel(data= wa_cities2,aes(x=long, y=lat, label=city)) +
  geom_point(data = wa_cities2, aes(x = long, y = lat), size = 3) +  
  xlim(110,130)+
  ylim(-36,-27)+
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()+
  theme(axis.title.x = element_text(size = 14),
axis.title.y = element_text(size = 14),
axis.text = element_text(size = 10))

##Plotting blocks
ggplot(block_coords, aes(group = openvshade))+
geom_point(aes(x = long, y = lat, shape=openvshade), size = 5, stroke = 1.5) +
  scale_shape_manual(values=c(1, 16))+
  xlim(116.1965,116.199)+
  ylim(-29.479,-29.480)+
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_classic()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position="top")


#### Making predictive plot of effect of watering treatment ####
### Survival 
#mean number of neighbours (min value for std_totalabund is -0.79)
#no dodder
#mean PC1 and PC2
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            Treatment:std_PC1 + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), arcadata)
## predicts for dry, ambient and wet 
arca_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, c(1, 0, 0)*0, c(0, 0, 1)*0, 0*0), 
                           se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
## add the category
arca_water_pred$treatment<-c("dry", "ambient", "wet")
#add species name
arca_water_pred$Species <- 'arca'
#assign x values
arca_water_pred$x <- c(0.5, 1, 1.5)

plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, arca_water_pred)
arrows(x0=arca_water_pred$x, y0=arca_water_pred$lower, x1=arca_water_pred$x, y1=arca_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)


#How much data do we have
#arca_no_nbh <- arcadata %>% filter(Neighbours01 == 0)
# ggplot(arcadata, aes(x = Treatment, y= surv_to_produce_seeds))+
#   geom_boxplot()+
#   geom_jitter(alpha=0.4)+
#   theme_classic()

#hygl
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                 Treatment:std_PC1 + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
hygl_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, c(1, 0, 0)*0, c(0, 0, 1)*0), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_water_pred$treatment<-c("dry", "ambient", "wet")
hygl_water_pred$Species <- 'hygl'
hygl_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, hygl_water_pred)
arrows(x0=hygl_water_pred$x, y0=hygl_water_pred$lower, x1=hygl_water_pred$x, y1=hygl_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)
#laro
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), larodata)
laro_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_water_pred$treatment<-c("dry", "ambient", "wet")
laro_water_pred$Species <- 'laro'
laro_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, laro_water_pred)
arrows(x0=laro_water_pred$x, y0=laro_water_pred$lower, x1=laro_water_pred$x, y1=laro_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)
#peai
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                 std_PC1:std_logp1_totalabund + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), peaidata)
peai_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, 0*0), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_water_pred$treatment<-c("dry", "ambient", "wet")
peai_water_pred$Species <- 'peai'
peai_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, peai_water_pred)
arrows(x0=peai_water_pred$x, y0=peai_water_pred$lower, x1=peai_water_pred$x, y1=peai_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)
#plde
pldesurvfinalmod
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                 I(std_PC1^2) + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), pldedata)
plde_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, 0^2, 0*0), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_water_pred$treatment<-c("dry", "ambient", "wet")
plde_water_pred$Species <- 'plde'
plde_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, plde_water_pred)
arrows(x0=plde_water_pred$x, y0=plde_water_pred$lower, x1=plde_water_pred$x, y1=plde_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)
#trcy
trcysurvfinalmod
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                 I(std_logp1_totalabund^2) + Treatment:std_PC1 + Treatment:I(std_logp1_totalabund^2) + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trcydata)
trcy_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, 0^2, c(1, 0, 0)*0, c(0, 0, 1)*0, c(1, 0, 0)*0^2, c(0, 0, 1)*0^2), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_water_pred$treatment<-c("dry", "ambient", "wet")
trcy_water_pred$Species <- 'trcy'
trcy_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, trcy_water_pred)
arrows(x0=trcy_water_pred$x, y0=trcy_water_pred$lower, x1=trcy_water_pred$x, y1=trcy_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)
#tror
trorsurvfinalmod
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                 Treatment:std_PC1 + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trordata)
tror_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, c(1, 0, 0)*0, c(0, 0, 1)*0), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_water_pred$treatment<-c("dry", "ambient", "wet")
tror_water_pred$Species <- 'tror'
tror_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, tror_water_pred)
arrows(x0=tror_water_pred$x, y0=tror_water_pred$lower, x1=tror_water_pred$x, y1=tror_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)
#vero
verosurvfinalmod
model <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                 Treatment:std_PC1 + (1|Site/Plot), 
               family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), verodata)
vero_water_pred<-glmm.predict(mod=model, newdat=data.frame(1, c(1, 0, 0), c(0, 0, 1), 0, 0, 0, 0, c(1, 0, 0)*0, c(0, 0, 1)*0), 
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_water_pred$treatment<-c("dry", "ambient", "wet")
vero_water_pred$Species <- 'vero'
vero_water_pred$x <- c(0.5, 1, 1.5)
plot(y ~ x, ylab="", xlab="", ylim=c(0, 1), tck=-0.01, vero_water_pred)
arrows(x0=vero_water_pred$x, y0=vero_water_pred$lower, x1=vero_water_pred$x, y1=vero_water_pred$upper, code=3, angle=90, length=0.1, lwd=2)

##Combine them
pred_water_all <- rbind(arca_water_pred, hygl_water_pred, laro_water_pred, peai_water_pred,
                     plde_water_pred, trcy_water_pred, tror_water_pred, vero_water_pred)
#Reorder to dry, ambient, wet
pred_water_all$treatment <- factor(pred_water_all$treatment, level = c("dry", "ambient", "wet"))

ggplot(pred_water_all, aes(x = Species, y = y, colour=treatment))+
  geom_point(position = position_dodge(0.8), cex=2.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.3), position = position_dodge(0.8), cex=1)+
  ylab("Probability of survival")+
  xlab("Species")+
  theme_classic()+
  my_theme+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = "italic"))

### Plotting interaction between watering and PC1 for survival ####
#Just for species that have PC1:Treatment term

dev.off()
pdf("Output/Figures/PC1_treatment_int.pdf", width=20, height=7)
par(mfrow=c(1,5), oma = c(12, 6, 4, 1), mar =c(2,2,1,1))

#arca
x_to_plot_dry <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcasurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb, x_to_plot_amb*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet, x_to_plot_wet*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry, x_to_plot_dry*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 6, col = "red")
#hygl
x_to_plot_dry <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglsurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.4,y = 0.95,"*", cex = 6, col = "red")
text(x = 1,y = 0.95,"*", cex = 6, col = "blue")
#trcy
x_to_plot_dry <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcysurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0^2, 0*x_to_plot_amb, 0*x_to_plot_amb, 0*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0^2, 0*x_to_plot_wet, 1*x_to_plot_wet, 0*0^2, 1*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 0^2, 1*x_to_plot_dry, 0*x_to_plot_dry, 1*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 6, col = "blue")
#tror
x_to_plot_dry <- seq.func(trordata$std_PC1[trordata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(trordata$std_PC1[trordata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(trordata$std_PC1[trordata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorsurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 6, col = "blue")
#vero
x_to_plot_dry <- seq.func(verodata$std_PC1[verodata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(verodata$std_PC1[verodata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(verodata$std_PC1[verodata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verosurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#labels
mtext("PC1 (std)", adj = 0.07, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.29, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.51, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.72, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.95, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line = 2)
mtext(~italic("A. calendula"), adj = 0.04, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = 0.27, padj=0.2, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = 0.5, padj=0.2, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = 0.71, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = 0.95, side = 3, cex = 2.5, outer = TRUE)
reset()
legend("bottom", title=NULL, horiz=T, legend=c("Dry", "Ambient", "Wet"),
       col=c("#CC79A7", "black", "#0072B2"), pch=19, cex=3, bty="n")
dev.off()

####### Same as above but no G. rosea
dev.off()
pdf("Output/Figures/PC1_treatment_int_3.pdf", width=20, height=7)
par(mfrow=c(1,5), oma = c(12, 6, 4, 1), mar =c(2,2,1,1))

#arca
x_to_plot_dry <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcasurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb, x_to_plot_amb*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet, x_to_plot_wet*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry, x_to_plot_dry*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 6, col = "red")
#hygl
x_to_plot_dry <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglsurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.4,y = 0.95,"*", cex = 6, col = "red")
text(x = 1,y = 0.95,"*", cex = 6, col = "blue")
#trcy
x_to_plot_dry <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcysurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0^2, 0*x_to_plot_amb, 0*x_to_plot_amb, 0*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0^2, 0*x_to_plot_wet, 1*x_to_plot_wet, 0*0^2, 1*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 0^2, 1*x_to_plot_dry, 0*x_to_plot_dry, 1*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 6, col = "blue")
#tror
x_to_plot_dry <- seq.func(trordata$std_PC1[trordata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(trordata$std_PC1[trordata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(trordata$std_PC1[trordata$Treatment=='Wet'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorsurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 6, col = "blue")
#labels
mtext("PC1 (std)", adj = 0.07, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.29, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.51, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("PC1 (std)", adj = 0.72, side = 1, line = 3, cex = 2.5,outer = TRUE)
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line = 2)
mtext(~italic("A. calendula"), adj = 0.04, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = 0.27, padj=0.2, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = 0.5, padj=0.2, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = 0.71, side = 3, cex = 2.5, outer = TRUE)
reset()
legend("right", title="Watering treatment", horiz=F, legend=c("Dry", "Ambient", "Wet"),
       col=c("#CC79A7", "black", "#0072B2"), pch=19, cex=2.5, bty="n")
dev.off()

## How does trcy dry:neighbour presence interaction look?
ggplot(lambdatrcy, aes(x = Neighbours01, y=lambda, colour = Treatment))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  theme_classic()

### Interaction survival, seed production and pop growth ~ PC1 and watering #####

##Data frame for wet and dry colours
colours <- data.frame(Treatment=c("Dry","Ambient","Wet"),
                         colour = c("#CC79A7","grey60","#0072B2"))
dev.off()
pdf("Output/Figures/supp_PC1_watering.pdf", width=20, height=18)
par(mfrow=c(3,4), oma = c(2, 6, 2, 1), mgp=c(5.5,1.5,0), mar =c(8,3,4,2))
##survival ##
#arca
#Merge colours with our data
arcadata <- left_join(arcadata, colours)
x_to_plot_dry <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(arcadata$std_PC1[arcadata$Treatment=='Wet'])
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(as.vector(arcadata$colour), 0.3), ylab=NA, xlab = "PC1 (std)", tck=-0.01, cex= 3, cex.axis= 3.5, cex.lab = 4, arcadata)
model <- arcasurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb, x_to_plot_amb*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet, x_to_plot_wet*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry, x_to_plot_dry*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 10, col = "red")
#hygl
hygldata <- left_join(hygldata, colours)
x_to_plot_dry <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(hygldata$std_PC1[hygldata$Treatment=='Wet'])
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(as.vector(hygldata$colour), 0.3), ylab=NA, xlab = "PC1 (std)", cex.lab=4, tck=-0.01, cex= 3, cex.axis= 3.5, hygldata)
model <- hyglsurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 10, col = "red")
text(x = 0.8,y = 0.95,"*", cex = 10, col = "blue")
#trcy
trcydata <- left_join(trcydata, colours)
x_to_plot_dry <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(trcydata$std_PC1[trcydata$Treatment=='Wet'])
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(as.vector(trcydata$colour), 0.3), ylab=NA, xlab = "PC1 (std)", cex.lab=4, tck=-0.01, cex= 3, cex.axis= 3.5, trcydata)
model <- trcysurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0^2, 0*x_to_plot_amb, 0*x_to_plot_amb, 0*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0^2, 0*x_to_plot_wet, 1*x_to_plot_wet, 0*0^2, 1*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 0^2, 1*x_to_plot_dry, 0*x_to_plot_dry, 1*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 10, col = "blue")
#tror
trordata <- left_join(trordata, colours)
x_to_plot_dry <- seq.func(trordata$std_PC1[trordata$Treatment=='Dry'])
x_to_plot_amb <- seq.func(trordata$std_PC1[trordata$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(trordata$std_PC1[trordata$Treatment=='Wet'])
plot(jitter(surv_to_produce_seeds, 0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(as.vector(trordata$colour), 0.3), ylab=NA, xlab = "PC1 (std)", cex.lab=4, tck=-0.01, cex= 3, cex.axis= 3.5, trordata)
model <- trorsurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 0.95,"*", cex = 10, col = "blue")
#labels
mtext("Probability of survival", side = 2, cex = 2.5, outer=TRUE, line = 2, adj=0.95)
mtext(~italic("A. calendula"), adj = 0.08, side = 3, padj=0.85, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = 0.36, padj=0.95, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = 0.64, padj=0.95, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = 0.91, side = 3, padj=0.85, cex = 2.5, outer = TRUE)
##seed production ##
#Placing plot centre left
par(mfg=c(2,1))
#tror
seedtror <- left_join(seedtror, colours)
x_to_plot_dry <- seq.func(seedtror$std_PC1[seedtror$Treatment=='Dry'])
x_to_plot_amb <- seq.func(seedtror$std_PC1[seedtror$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(seedtror$std_PC1[seedtror$Treatment=='Wet'])
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.6), ylim=c(1,100), log = "y", pch=19, col=alpha(as.vector(seedtror$colour), 0.3), ylab=NA, xlab="PC1 (std)", tck=-0.01, cex= 3, cex.lab = 4, cex.axis = 3.5, seedtror)
model <- trorseedfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.5,y = 70,"*", cex = 10, col = "red")
#labels
mtext("Number of viable seeds (log+1)", side = 2, cex = 2.5, outer=TRUE, line = 2)
mtext(~italic("T. ornata"), adj = 0.08, side = 3, cex = 2.5, outer = TRUE, padj=20.5)
## population growth rate ##
# peai
lambdapeai <- left_join(lambdapeai, colours)
#Lambda 
par(mfg=c(3,1))
x_to_plot_dry <- seq.func(lambdapeai$std_PC1[lambdapeai$Treatment=='Dry'])
x_to_plot_amb <- seq.func(lambdapeai$std_PC1[lambdapeai$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(lambdapeai$std_PC1[lambdapeai$Treatment=='Wet'])
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(as.vector(lambdapeai$colour), 0.3), ylab=NA, xlab="PC1 (std)", tck=-0.01, cex= 3, cex.lab = 4, cex.axis = 3.5, lambdapeai)
model <- peailambdafinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, x_to_plot_amb^2, 0*x_to_plot_amb, 0*x_to_plot_amb, x_to_plot_amb*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, x_to_plot_amb^2, 0*x_to_plot_wet, 1*x_to_plot_wet, x_to_plot_wet*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, x_to_plot_amb^2, 1*x_to_plot_dry, 0*x_to_plot_dry, x_to_plot_wet*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 4,"*", cex = 10, col = "red")
#vero
lambdavero <- left_join(lambdavero, colours)
par(mfg=c(3,2))
x_to_plot_dry <- seq.func(lambdavero$std_PC1[lambdavero$Treatment=='Dry'])
x_to_plot_amb <- seq.func(lambdavero$std_PC1[lambdavero$Treatment=='Ambient'])
x_to_plot_wet <- seq.func(lambdavero$std_PC1[lambdavero$Treatment=='Wet'])
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(as.vector(lambdavero$colour), 0.3), ylab=NA, xlab="PC1 (std)", tck=-0.01, cex= 3, cex.lab = 4, cex.axis = 3.5, lambdavero)
model <- verolambdafinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb, x_to_plot_amb*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet, x_to_plot_wet*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#dry - red
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 1, 0, x_to_plot_dry, 0, 0, 1*x_to_plot_dry, 0*x_to_plot_dry, x_to_plot_wet*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_dry, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
text(x = 1.3,y = 3.4,"*", cex = 10, col = "blue")
#labels
mtext("Population growth rate (log)", side = 2, cex = 2.5, outer=TRUE, line = 2, adj=0.03)
mtext(~italic("P. airoides"), adj = 0.08, side = 3, cex = 2.5, outer = TRUE, padj = 39.5)
mtext(~italic("G. rosea"), adj = 0.37, side = 3, cex = 2.5, outer = TRUE, padj=38.5)
reset()
legend(x = 0.65, y = 0.25, title="Watering treatment", horiz=F, legend=c("Dry", "Ambient", "Wet"),
       col=c("#CC79A7", "black", "#0072B2"), pch=19, cex=3.5, bty="n")
dev.off()

### Summarising how many plants of which species germinated early in Feb #####
febdata <- vitaldata %>% select(rowID, Site, Plot, Species, C_E_or_T, Rep, February_germination) %>%
  filter(February_germination > 0)
febdata %>% group_by(Species) %>% summarise(number_plots = n(),
                                            number_plants = sum(February_germination))


#### After JoE review - exploring single abiotic factors####
# Instead of responses ~ PC1:
#responses ~ canopy cover, soil nutrients (phosphorous, ammonium, nitrate, and potassium content), litter cover
# Total N instead? And P?
# losing soil pH

#watering treatment
#N
#P
#canopy cover

#litter cover: 0, 12.5, 25, 50
test <- vitaldata %>% group_by(Site, Plot) %>% filter(row_number()==1)
ggplot(test, aes(x=Litter_cover_percent, y=cc_percentage))+
         geom_jitter(alpha=0.3, width=0.5, height=0)+
  geom_smooth()+
  theme_classic()

lm<- lm(cc_percentage ~ Litter_cover_percent, test)
summary(lm)
r.squaredGLMM(lm)
#They are correlated. No canopy, no litter!

#Sum, total N
soildata <- soildata %>% mutate(total_n = `NH4-N`+`NO3-N`)
hist(test$log_totalN)
totalndata <- soildata %>% select(Site, Plot, total_n)
#There is one plot that is a massive outlier - but I think this 
# is the plot with orchids growing in it! Was closer to a euc
vitaldata <- left_join(vitaldata, totalndata, join_by(Site,Plot))

test <- vitaldata %>% group_by(Site, Plot) %>% filter(row_number()==1)
ggplot(test, aes(x=log(total_n.x), y=Litter_cover_percent))+
  geom_jitter(alpha=0.3, width=0.5, height=0)+
 # geom_smooth()+
  theme_classic()

#hm, total N and canopy cover are correlated too
lm<- lm(cc_percentage ~ total_n.x, test)
summary(lm)
r.squaredGLMM(lm)

ggplot(test, aes(y=cc_percentage, x=Cover))+
  geom_boxplot()+
  geom_jitter(alpha=0.3, width=0.05)+
  theme_classic()

#Sun (0-12%) and shad (58-100%) cover, ranges per plot
#site 6B: log(0.2)+log(2.1)=-0.87
#When would you log values before versus after adding them??
#negative value, so adding +1 before analysis
hist(test$log_totalN+1)

ggplot(test, aes(y=log_totalN+1, x=Cover))+
  geom_boxplot()+
  geom_jitter(alpha=0.3, width=0.05)+
  theme_classic()

lm<- lm(log_totalN+1 ~ Cover, test)
summary(lm)
r.squaredGLMM(lm)

#soil P
hist(test$log_P)
ggplot(test, aes(y=log_P, x=Cover))+
  geom_boxplot()+
  geom_jitter(alpha=0.3, width=0.05)+
  theme_classic()

lm<- lm(log_P ~ Cover, test)
summary(lm)
r.squaredGLMM(lm)

#Phosphorous is highly correlated to shade. High P in the shade

#pH?
hist(test$pH)
ggplot(test, aes(y=pH, x=Cover))+
  geom_boxplot()+
  geom_jitter(alpha=0.3, width=0.05)+
  theme_classic()

lm<- lm(pH ~ Cover, test)
summary(lm)
r.squaredGLMM(lm)


#Should neighbour abund be continuous?
#I think not - too many zeroes, not normally distributed
#Already did this decision-making for previous version though, chose to keep it continuous
hist(seedmodeldata$std_logp1_totalabund)
hist(seedmodeldata$Total_abundance)
ggplot(test2, aes(x=Total_abundance, y=log(No_viable_seeds_grouped+1)))+
  geom_jitter(alpha=0.3, height=0.05, width=0.05)+
  theme_classic()

test2 <- vitaldata %>% filter(Total_abundance>0)
test3 <- vitaldata %>% filter(Total_abundance==0)
#534 focals with at least 1 neighbour
#650 focals with no neighbours
hist(test2$Total_abundance)

hist(vitaldata$std_N)

# ggplot(vitaldata, aes(x=Site, y= std_N))+
#   geom_jitter()+
#   theme_classic()
ggplot(seedmodeldata, aes(x=Site, y= Total_abundance))+
  geom_jitter()+
  theme_classic()

test <- vitaldata %>% filter(Total_abundance>0)
hist(log(seedmodeldata$Total_abundance+1))
hist(test$Total_abundance)


#### Exploring quadratic terms - N, pH, neighbour abund ####

# Could do it by species, but would rather use loops:
# #ARCA - linear PC1 survival
# arcasurvquad1 <- glmer(surv_to_produce_seeds ~ std_PC1 + I(std_PC1^2) + (1|Site/Plot), family = binomial, arcadata)
# summary(arcasurvquad1)

### Germination ###
## PEAI - quad N germ
## All other species - linear PC1 germination
for (i in 1:length(specieslist)){
  print(specieslist[i])
  germquadPC1 <- glmer(cbind(total_germ, total_no_germ) ~ std_N + I(std_N^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(germquadPC1))
}
#Plotting if quadratic to double check
#This produces a plot that plots quadratic responses if the quadratic term is < 0.05, otherwise linear plot
dev.off()
pdf("Output/Figures/germ_N_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$percent_germ~plotted.data$std_N, pch=19, col="grey60", ylab="Probability of germination", xlab="Total N (logged, standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_N) 
  model<-glmer(cbind(total_germ, total_no_germ)~std_N + I(std_N^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(cbind(total_germ, total_no_germ)~std_N + (1|Site/Plot), family = binomial, plotted.data)
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


## All species - linear pH germination (almost signif hygl, laro, plde)
for (i in 1:length(specieslist)){
  print(specieslist[i])
  germquadPC2 <- glmer(cbind(total_germ, total_no_germ) ~ std_pH + I(std_pH^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(germquadPC2))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/germ_pH_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$percent_germ~plotted.data$std_pH, pch=19, col="grey60", ylab="Probability of germination", xlab="PC2 (standardised)", cex.lab=2, .axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13)
  title(main=bquote(italic(.(species.name.list[i]))), .main=2.5)
  x_to_plot<-seq.func(plotted.data$std_pH) 
  model<-glmer(cbind(total_germ, total_no_germ)~std_pH + I(std_pH^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(cbind(total_germ, total_no_germ)~std_pH + (1|Site/Plot), family = binomial, plotted.data)
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
## All species - linear N survival
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survquadPC1 <- glmer(surv_to_produce_seeds ~ std_N + I(std_N^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(survquadPC1))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/surv_N_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_N, pch=19, col="grey60", ylab="Probability of survival", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_N) 
  model<-glmer(surv_to_produce_seeds~std_N + I(std_N^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(surv_to_produce_seeds~std_N + (1|Site/Plot), family = binomial, plotted.data)
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

##ARCA - quadratic pH survival 
##All other species - linear pH survival
for (i in 1:length(specieslist)){
  print(specieslist[i])
  survquadPC2 <- glmer(surv_to_produce_seeds ~ std_pH + I(std_pH^2) + (1|Site/Plot), 
                       family = binomial, data = filter(vitaldata, Species == specieslist[i]))
  print(summary(survquadPC2))
}
#Plotting if quadratic to double check
#TRCY linear model not converging
dev.off()
pdf("Output/Figures/surv_pH_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.s)){
  plotted.data<-as.data.frame(species.list.s[i])
  plot(plotted.data$surv_to_produce_seeds~plotted.data$std_pH, pch=19, col="grey60", ylab="Probability of survival", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_pH) 
  model<-glmer(surv_to_produce_seeds~std_pH + I(std_pH^2) + (1|Site/Plot), family = binomial, plotted.data)
  model2<-glmer(surv_to_produce_seeds~std_pH + (1|Site/Plot), family = binomial, plotted.data)
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

### Fecundity ###
#Note that I have to use different datasets - seedmodeldata/seedarca/species.list.f

## POLE - quadratic N fecundity (but don't analyse pole seed production data, unreliable and small sample size)
## All other species - linear N fecundity
for (i in 1:length(specieslist)){
  print(specieslist[i])
  seedquadPC1 <- glmmTMB(No_viable_seeds_grouped ~ std_N + I(std_N^2) + (1|Site/Plot), 
                         family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(seedquadPC1))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/SP_N_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_N, pch=19, col="grey60", ylab="Number of seeds produced", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_N) 
  model<-glmmTMB(No_viable_seeds_grouped~std_N + I(std_N^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  model2<-glmmTMB(No_viable_seeds_grouped~std_N + (1|Site/Plot), family = nbinom2, plotted.data)
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

## All species - linear pH fecundity
for (i in 1:length(specieslist)){
  print(specieslist[i])
  seedquadPC2 <- glmmTMB(No_viable_seeds_grouped ~ std_pH + I(std_pH^2) + (1|Site/Plot), 
                         family = nbinom2, data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(seedquadPC2))
}
#Plotting if quadratic to double check
dev.off()
pdf("Output/Figures/SP_pH_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.f)){
  plotted.data<-as.data.frame(species.list.f[i])
  plot(plotted.data$No_viable_seeds_grouped~plotted.data$std_pH, pch=19, col="grey60", ylab="Number of seeds produced", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_pH) 
  model<-glmmTMB(No_viable_seeds_grouped~std_pH + I(std_pH^2) + (1|Site/Plot), family = nbinom2, plotted.data)
  model2<-glmmTMB(No_viable_seeds_grouped~std_pH + (1|Site/Plot), family = nbinom2, plotted.data)
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

## All  species - linear N lambda
for (i in 1:length(specieslist)){
  print(specieslist[i])
  lambdaquadPC1 <- lmer(log_lambda ~ std_N + I(std_N^2) + (1|Site/Plot), 
                        data = filter(popdata, Species == specieslist[i]))
  print(summary(lambdaquadPC1))
}
#Plotting if quadratic to double check
#This produces a plot that plots quadratic responses if the quadratic term is < 0.05, otherwise linear plot
#summary(mod1)$coefficients[3,5] - these are the coefficient coordinates we want

dev.off()
pdf("Output/Figures/lambda_N_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.l)){
  plotted.data<-as.data.frame(species.list.l[i])
  plot(plotted.data$log_lambda~plotted.data$std_N, pch=19, col="grey60", ylab="Population growth rate (logged plus 1)", xlab="PC1 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_N) 
  model<-lmer(log_lambda~std_N + I(std_N^2) + (1|Site/Plot), plotted.data)
  model2<-lmer(log_lambda~std_N + (1|Site/Plot), plotted.data)
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

## All species - linear pH lambda
for (i in 1:length(specieslist)){
  print(specieslist[i])
  lambdaquadPC2 <- lmer(log_lambda ~ std_pH + I(std_pH^2) + (1|Site/Plot), 
                        data = filter(popdata, Species == specieslist[i]))
  print(summary(lambdaquadPC2))
}
#Plotting if quadratic to double check

dev.off()
pdf("Output/Figures/lambda_pH_if_quad.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(4,6,2,1))
par(pty="s")
for(i in 1:length(species.list.l)){
  plotted.data<-as.data.frame(species.list.l[i])
  plot(plotted.data$log_lambda~plotted.data$std_pH, pch=19, col="grey60", ylab="Population growth rate (logged plus 1)", xlab="PC2 (standardised)", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  x_to_plot<-seq.func(plotted.data$std_pH) 
  model<-lmer(log_lambda~std_pH + I(std_pH^2) + (1|Site/Plot), plotted.data)
  model2<-lmer(log_lambda~std_pH + (1|Site/Plot), plotted.data)
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

#### After JoE review - single abiotic factors analysis ####
#Cover (open/shaded) and total soil N and pH
#Cover instead of PC1
#pH instead of PC2
#soil N as covariate (with pH)

#### New germination models ####
#Same overall pattern as PC1 with cover

## original models
## ARCA - has quadratic germ ~ PC2
arcagermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ std_PC1 + std_PC2 + I(std_PC2^2) + (1|Site/Plot/rowID), 
                          family = binomial, arcadata)
arcagermmod1dharma <- simulateResiduals(arcagermfinalmod)
plot(arcagermmod1dharma)
vif(arcagermfinalmod)
summary(arcagermfinalmod)
testDispersion(arcagermfinalmod)

## ARCA new models
#N signif
arcagermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, arcadata)
arcagermmod1dharma <- simulateResiduals(arcagermfinalmod)
plot(arcagermmod1dharma)
vif(arcagermfinalmod)
summary(arcagermfinalmod)
testDispersion(arcagermfinalmod)

## hygl new models
#Cover signif
hyglgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, hygldata)
hyglgermmod1dharma <- simulateResiduals(hyglgermfinalmod)
plot(hyglgermmod1dharma)
vif(hyglgermfinalmod)
summary(hyglgermfinalmod)

## laro new models
#Cover signif
larogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, larodata)
larogermmod1dharma <- simulateResiduals(larogermfinalmod)
plot(larogermmod1dharma)
vif(larogermfinalmod)
summary(larogermfinalmod)

## peai new models
#Cover signif
peaigermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, peaidata)
peaigermmod1dharma <- simulateResiduals(peaigermfinalmod)
plot(peaigermmod1dharma)
vif(peaigermfinalmod)
summary(peaigermfinalmod)

## plde new models
#pH signif
pldegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, pldedata)
pldegermmod1dharma <- simulateResiduals(pldegermfinalmod)
plot(pldegermmod1dharma)
vif(pldegermfinalmod)
summary(pldegermfinalmod)

## pole new models
#N signif
polegermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, poledata)
polegermmod1dharma <- simulateResiduals(polegermfinalmod)
plot(polegermmod1dharma)
vif(polegermfinalmod)
summary(polegermfinalmod)

## trcy new models
#Cover signif
trcygermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, trcydata)
trcygermmod1dharma <- simulateResiduals(trcygermfinalmod)
plot(trcygermmod1dharma)
vif(trcygermfinalmod)
summary(trcygermfinalmod)

## tror new models
#Cover signif
trorgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, trordata)
trorgermmod1dharma <- simulateResiduals(trorgermfinalmod)
plot(trorgermmod1dharma)
vif(trorgermfinalmod)
summary(trorgermfinalmod)

## vero new models
#Cover signif
verogermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, verodata)
verogermmod1dharma <- simulateResiduals(verogermfinalmod)
plot(verogermmod1dharma)
vif(verogermfinalmod)
summary(verogermfinalmod)

#### New survival models ####

## ARCA - ORIGINAL
arcasurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                        Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, arcadata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
arcasurvmod1dharma <- simulateResiduals(arcasurvmod1)
plot(arcasurvmod1dharma)
#good
vif(arcasurvmod1)
summary(arcasurvmod1)
# Model simplification step - Significant Dry:PC1 and PC1:NA, removing all others
arcasurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
                            Treatment:std_PC1 + std_PC1:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), arcadata)
arcasurvfinalmoddharma <- simulateResiduals(arcasurvfinalmod)
plot(arcasurvfinalmoddharma)
#good
summary(arcasurvfinalmod)

## New - ARCA
arcasurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, arcadata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
arcasurvmod1dharma <- simulateResiduals(arcasurvmod1)
plot(arcasurvmod1dharma)
#good
vif(arcasurvmod1)
summary(arcasurvmod1)

#N matters for survival, higher N higher survival. pH did not

#Simplifying -- watering DOES interact with cover
# water NO interaction with neighbour abundance
# cover DOES interact with neighbour abundance
#fundamentally the same as the original model

arcasurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Cover:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), arcadata)
arcasurvfinalmoddharma <- simulateResiduals(arcasurvfinalmod)
plot(arcasurvfinalmoddharma)
#good
summary(arcasurvfinalmod)

## hygl
hyglsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, hygldata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
hyglsurvmod1dharma <- simulateResiduals(hyglsurvmod1)
plot(hyglsurvmod1dharma)
#good
vif(hyglsurvmod1)
summary(hyglsurvmod1)

ggplot(hygldata, aes(y=surv_to_produce_seeds, x=Treatment, colour = Cover))+
  geom_boxplot()+
  geom_jitter(alpha=0.3, width=0.1, height=0.2)+
  theme_classic()
#all dry plants in the shade died. HUGE estimates and SE from Dry:Shade

#N matters for survival, higher N higher survival. pH did not
#Simplyfing -- watering DOES interact with cover
# watering does interact with neighbour abundance (different from orig. model, not signif in final mod)
# cover no interaction with neighbour abundance

hyglsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Treatment:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), hygldata)
hyglsurvfinalmoddharma <- simulateResiduals(hyglsurvfinalmod)
plot(hyglsurvfinalmoddharma)
vif(hyglsurvfinalmod)
#good
summary(hyglsurvfinalmod)
r.squaredGLMM(hyglsurvfinalmod)

## laro
larosurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, larodata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
larosurvmod1dharma <- simulateResiduals(larosurvmod1)
plot(larosurvmod1dharma)
#good
vif(larosurvmod1)
summary(larosurvmod1)

#Nothing signif
#fundamentally the same as the original model

larosurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), larodata)
larosurvfinalmoddharma <- simulateResiduals(larosurvfinalmod)
plot(larosurvfinalmoddharma)
#good
summary(larosurvfinalmod)

## peai
peaisurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, peaidata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
peaisurvmod1dharma <- simulateResiduals(peaisurvmod1)
plot(peaisurvmod1dharma)
#good
vif(peaisurvmod1)
summary(peaisurvmod1)

#Nothing signif - different from orig model, PC1:NA

peaisurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), peaidata)
peaisurvfinalmoddharma <- simulateResiduals(peaisurvfinalmod)
plot(peaisurvfinalmoddharma)
#good
summary(peaisurvfinalmod)

## plde
pldesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, pldedata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
pldesurvmod1dharma <- simulateResiduals(pldesurvmod1)
plot(pldesurvmod1dharma)
#good
vif(pldesurvmod1)
#VIF OF 8 for pH - look into this later
summary(pldesurvmod1)

#Nothing signif - different from orig model, PC1:NA

pldesurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), pldedata)
pldesurvfinalmoddharma <- simulateResiduals(pldesurvfinalmod)
plot(pldesurvfinalmoddharma)
#good
summary(pldesurvfinalmod)

## trcy
trcysurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, trcydata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
trcysurvmod1dharma <- simulateResiduals(trcysurvmod1)
plot(trcysurvmod1dharma)
#good
vif(trcysurvmod1)
summary(trcysurvmod1)

#Simplifying -- watering DOES interact with cover
# water NO interaction with neighbour abundance
# cover DOES interact with neighbour abundance
#different from the original model

trcysurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Cover:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trcydata)
trcysurvfinalmoddharma <- simulateResiduals(trcysurvfinalmod)
plot(trcysurvfinalmoddharma)
#good
summary(trcysurvfinalmod)

## tror
trorsurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, trordata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
trorsurvmod1dharma <- simulateResiduals(trorsurvmod1)
plot(trorsurvmod1dharma)
#good
vif(trorsurvmod1)
summary(trorsurvmod1)

#Simplifying -- watering DOES interact with cover
# water NO interaction with neighbour abundance
# cover NO interaction with neighbour abundance
#fundamentally the same as the original model

trorsurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trordata)
trorsurvfinalmoddharma <- simulateResiduals(trorsurvfinalmod)
plot(trorsurvfinalmoddharma)
#good
summary(trorsurvfinalmod)

## vero
verosurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                        Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot),
                      family = binomial, verodata, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
#Not converging without optimiser
verosurvmod1dharma <- simulateResiduals(verosurvmod1)
plot(verosurvmod1dharma)
#good
vif(verosurvmod1)
summary(verosurvmod1)

#Nothing signif
#fundamentally the same as the original model

verosurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), verodata)
verosurvfinalmoddharma <- simulateResiduals(verosurvfinalmod)
plot(verosurvfinalmoddharma)
#good
summary(verosurvfinalmod)

##### New seed production models ####

#ARCA
arcaseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedarca)
arcaseedmod1dharma <- simulateResiduals(arcaseedmod1)
plot(arcaseedmod1dharma)
summary(arcaseedmod1)
#simplifying - signif Treatment:NA
#same as original
arcaseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                              Treatment:std_logp1_totalabund + (1|Site/Plot), 
                            family = nbinom2, seedarca)
arcaseedfinalmoddharma <- simulateResiduals(arcaseedfinalmod)
plot(arcaseedfinalmoddharma)
#not good
summary(arcaseedfinalmod)

#hygl
#won't converge with this new model
hyglseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedhygl)
hyglseedmod1dharma <- simulateResiduals(hyglseedmod1)
plot(hyglseedmod1dharma)
summary(hyglseedmod1)
#will converge using glmer.nb but drops TreatmentDry:CoverSun
#BECAUSE there is a complete separation issue - no dry plants in the shade produced viable seeds
#because there were no dry plants in the shade - they all died
#same estimate values between models excpt for TreatmentDry estimate
#glmmTMB optimiser stops and we don't get probabilities
#glmer.nb just drops that row from the output.
#Cover:Treatment not significant anyway so dropping it from final model
hyglseedmod2 <- glmer.nb(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                           Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), seedhygl)
hyglseedmod2dharma <- simulateResiduals(hyglseedmod2)
plot(hyglseedmod2dharma)
summary(hyglseedmod2)

# ggplot(seedhygl, aes(x=Cover, y= No_viable_seeds_grouped, colour = Treatment)) +
#   geom_jitter(alpha=0.5)+
#   theme_classic()

#simplifying - nothing signif
hyglseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedhygl)
hyglseedfinalmoddharma <- simulateResiduals(hyglseedfinalmod)
plot(hyglseedfinalmoddharma)
#good
summary(hyglseedfinalmod)
#same as original

#peai
peaiseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedpeai)
peaiseedmod1dharma <- simulateResiduals(peaiseedmod1)
plot(peaiseedmod1dharma)
summary(peaiseedmod1)

#simplifying - nothing signif, different from original
peaiseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedpeai)
peaiseedfinalmoddharma <- simulateResiduals(peaiseedfinalmod)
plot(peaiseedfinalmoddharma)
#good
summary(peaiseedfinalmod)

#laro
laroseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedlaro)
laroseedmod1dharma <- simulateResiduals(laroseedmod1)
plot(laroseedmod1dharma)
summary(laroseedmod1)

#simplifying - nothing signif (Cover:NA 0.054!), different from original
laroseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedlaro)
laroseedfinalmoddharma <- simulateResiduals(laroseedfinalmod)
plot(laroseedfinalmoddharma)
#good
summary(laroseedfinalmod)

#plde
pldeseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedplde)
pldeseedmod1dharma <- simulateResiduals(pldeseedmod1)
plot(pldeseedmod1dharma)
summary(pldeseedmod1)

#simplifying - nothing signif, same as original
pldeseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedplde)
pldeseedfinalmoddharma <- simulateResiduals(pldeseedfinalmod)
plot(pldeseedfinalmoddharma)
#good
summary(pldeseedfinalmod)

#trcy
trcyseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedtrcy)
trcyseedmod1dharma <- simulateResiduals(trcyseedmod1)
plot(trcyseedmod1dharma)
summary(trcyseedmod1)

#simplifying - signif Cover:NA (and Treatment:Cover 0.052!), differnt from original
trcyseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                              Cover:std_logp1_totalabund + (1|Site/Plot), 
                            family = nbinom2, seedtrcy)
trcyseedfinalmoddharma <- simulateResiduals(trcyseedfinalmod)
plot(trcyseedfinalmoddharma)
#good
summary(trcyseedfinalmod)

#tror
trorseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedtror)
trorseedmod1dharma <- simulateResiduals(trorseedmod1)
plot(trorseedmod1dharma)
summary(trorseedmod1)

#simplifying - signif Treatment:Cover
#same as original
trorseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                              Treatment:Cover + (1|Site/Plot), 
                            family = nbinom2, seedtror)
trorseedfinalmoddharma <- simulateResiduals(trorseedfinalmod)
plot(trorseedfinalmoddharma)
#good
summary(trorseedfinalmod)

#vero
veroseedmod1 <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                          Treatment:Cover + Treatment:std_logp1_totalabund + Cover:std_logp1_totalabund + (1|Site/Plot), 
                        family = nbinom2, seedvero)
veroseedmod1dharma <- simulateResiduals(veroseedmod1)
plot(veroseedmod1dharma)
summary(veroseedmod1)

#simplifying - nothing signif, same as original
veroseedfinalmod <- glmmTMB(No_viable_seeds_grouped ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 + (1|Site/Plot), 
                            family = nbinom2, seedvero)
veroseedfinalmoddharma <- simulateResiduals(veroseedfinalmod)
plot(veroseedfinalmoddharma)
#good
summary(veroseedfinalmod)

##### New population growth models ####
#original: arcalambdamod1 <- lmer(log_lambda ~ Treatment + std_PC1 + std_PC2 + Neighbours01 +
#Treatment:std_PC1 + Treatment:Neighbours01 + std_PC1:Neighbours01 + (1|Site/Plot), lambdaarca)

## arca
arcalambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdamod1dharma <- simulateResiduals(arcalambdamod1)
plot(arcalambdamod1dharma)
#not good
vif(arcalambdamod1)
summary(arcalambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
arcalambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdaarca)
arcalambdafinalmoddharma <- simulateResiduals(arcalambdafinalmod)
plot(arcalambdafinalmoddharma)
#okay
summary(arcalambdafinalmod)

## HYGL
hygllambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdamod1dharma <- simulateResiduals(hygllambdamod1)
plot(hygllambdamod1dharma)
#okay
vif(hygllambdamod1)
summary(hygllambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
hygllambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdahygl)
hygllambdafinalmoddharma <- simulateResiduals(hygllambdafinalmod)
plot(hygllambdafinalmoddharma)
#okay
summary(hygllambdafinalmod)

## laro
larolambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdalaro)
larolambdamod1dharma <- simulateResiduals(larolambdamod1)
plot(larolambdamod1dharma)
#okay
vif(larolambdamod1)
summary(larolambdamod1)
#Model simplification step - signif Cover:Neighbours, same as original
larolambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdalaro)
larolambdafinalmoddharma <- simulateResiduals(larolambdafinalmod)
plot(larolambdafinalmoddharma)
#okay
summary(larolambdafinalmod)

## peai
peailambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdapeai)
peailambdamod1dharma <- simulateResiduals(peailambdamod1)
plot(peailambdamod1dharma)
#okay
vif(peailambdamod1)
summary(peailambdamod1)
#Model simplification step - signif Treatment:Cover, Cover:Neighbours
#same as original
peailambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Treatment:Cover + Cover:Neighbours01 + (1|Site/Plot), lambdapeai)
peailambdafinalmoddharma <- simulateResiduals(peailambdafinalmod)
plot(peailambdafinalmoddharma)
#okay
summary(peailambdafinalmod)

## plde
pldelambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdamod1dharma <- simulateResiduals(pldelambdamod1)
plot(pldelambdamod1dharma)
#okay
vif(pldelambdamod1)
summary(pldelambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
pldelambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdaplde)
pldelambdafinalmoddharma <- simulateResiduals(pldelambdafinalmod)
plot(pldelambdafinalmoddharma)
#okay
summary(pldelambdafinalmod)

## trcy
trcylambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdamod1dharma <- simulateResiduals(trcylambdamod1)
plot(trcylambdamod1dharma)
#okay
vif(trcylambdamod1)
summary(trcylambdamod1)
#Model simplification step - signif Treatment:Neighbours and Cover:Neighbours
#same as original
trcylambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdatrcy)
trcylambdafinalmoddharma <- simulateResiduals(trcylambdafinalmod)
plot(trcylambdafinalmoddharma)
#okay
summary(trcylambdafinalmod)

## tror
trorlambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdatror)
trorlambdamod1dharma <- simulateResiduals(trorlambdamod1)
plot(trorlambdamod1dharma)
#okay
vif(trorlambdamod1)
summary(trorlambdamod1)
#Model simplification step - No interactions significant, removing all
#same as original
trorlambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + (1|Site/Plot), lambdatror)
trorlambdafinalmoddharma <- simulateResiduals(trorlambdafinalmod)
plot(trorlambdafinalmoddharma)
#okay
summary(trorlambdafinalmod)

## vero
verolambdamod1 <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 +
                         Treatment:Cover + Treatment:Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdavero)
verolambdamod1dharma <- simulateResiduals(verolambdamod1)
plot(verolambdamod1dharma)
#okay
vif(verolambdamod1)
summary(verolambdamod1)
#Model simplification step - signif Cover:Neighbours 
#different from original
verolambdafinalmod <- lmer(log_lambda ~ Treatment + Cover + std_pH + std_N + Neighbours01 + Cover:Neighbours01 + (1|Site/Plot), lambdavero)
verolambdafinalmoddharma <- simulateResiduals(verolambdafinalmod)
plot(verolambdafinalmoddharma)
#okay
summary(verolambdafinalmod)


#### New extracting marginal and conditional r squared values for table ####
#Extracting values for theoretical marginal R squared
germ_model_list <- list(arcagermfinalmod, hyglgermfinalmod, larogermfinalmod, peaigermfinalmod, pldegermfinalmod, polegermfinalmod, trcygermfinalmod, trorgermfinalmod, verogermfinalmod)
rsquared = lapply(1:length(germ_model_list), function(x) {
  as.data.frame(r.squaredGLMM(germ_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
germ_rsquared_table <- do.call("rbind", rsquared)

#Rename species
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_rsquared_table <- within(germ_rsquared_table, Species[Species == '9'] <- 'Goodenia rosea')

## Table output GERMINATION ####
## Extracting values for all in a loop
germ_model_list <- list(arcagermfinalmod, hyglgermfinalmod, larogermfinalmod, peaigermfinalmod, pldegermfinalmod, polegermfinalmod, trcygermfinalmod, trorgermfinalmod, verogermfinalmod)
effects = lapply(1:length(germ_model_list), function(x) {
  as.data.frame(coef(summary(germ_model_list[[x]]))) %>% mutate(Species=paste0(x))})
germ_effects_table <- do.call("rbind", effects)

#make rownames a column
germ_effects_table <- cbind(Effect = rownames(germ_effects_table), germ_effects_table)
rownames(germ_effects_table) <- NULL

#Rename species
germ_effects_table <- within(germ_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_effects_table <- within(germ_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_effects_table <- within(germ_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_effects_table <- within(germ_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_effects_table <- within(germ_effects_table, Species[Species == '5'] <- 'Plantago debilis')
germ_effects_table <- within(germ_effects_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_effects_table <- within(germ_effects_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_effects_table <- within(germ_effects_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_effects_table <- within(germ_effects_table, Species[Species == '9'] <- 'Goodenia rosea')

#Renaming effects since loop adding values to ends
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, '(Intercept)')] <- 'Intercept'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'std_pH')] <- 'pH'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'std_N')] <- 'Nitrogen'
germ_effects_table$Effect[startsWith(germ_effects_table$Effect, 'CoverShade')] <- 'Cover (shade)'

germ_effects_table <- within(germ_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
germ_effects_table <- within(germ_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
germ_effects_table <- within(germ_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
germ_effects_table <- within(germ_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
germ_effects_table <- within(germ_effects_table, Species[Species == '5'] <- 'Plantago debilis')
germ_effects_table <- within(germ_effects_table, Species[Species == '6'] <- 'Podolepis lessonii')
germ_effects_table <- within(germ_effects_table, Species[Species == '7'] <- 'Trachymene cyanopetala')
germ_effects_table <- within(germ_effects_table, Species[Species == '8'] <- 'Trachymene ornata')
germ_effects_table <- within(germ_effects_table, Species[Species == '9'] <- 'Goodenia rosea')

#Renaming columns
germ_effects_table <- germ_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')

#Making column with Estimate (+/- SE) and p value asterisks all combined
#Add column for asterisks based on below function
germ_effects_table <- germ_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
germ_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", germ_effects_table$Estimate, germ_effects_table$SE, germ_effects_table$p_asterisks)

#Join with r squared values
#germ_effects_table <- left_join(germ_effects_table, germ_rsquared_table)

germ_effects_kbl <- germ_effects_table %>% select(Species, Effect, collated)
germ_effects_kbl <- germ_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)

#Plotting with kableR
germ_effects_kbl %>% kbl(align = 'lccccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.13/0.82"=1, "0.25/0.91"=1, "0.29/0.84"=1, "0.26/0.89"=1, "0.23/0.91"=1, "0.28/0.96"=1, "0.12/0.92"=1, "0.34/0.90"=1, "0.35/0.93"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Emergence" = 1, "n=190"=1, "n=189"=1, "n=192"=1, "n=192"=1, "n=92"=1, "n=185"=1, "n=192"=1, "n=191"=1, "n=191"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:10, width = 4)


## SURVIVAL #### 

### PROBLEM TO FIX - Estimate -23.7 for hygl Dry:Cover shade.
### r squared values
surv_model_list <- list(arcasurvfinalmod, hyglsurvfinalmod, larosurvfinalmod, peaisurvfinalmod, pldesurvfinalmod, trcysurvfinalmod, trorsurvfinalmod, verosurvfinalmod)
rsquared = lapply(1:length(surv_model_list), function(x) {
  as.data.frame(r.squaredGLMM(surv_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
surv_rsquared_table <- do.call("rbind", rsquared)

#Rename species
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '7'] <- 'Trachymene ornata')
surv_rsquared_table <- within(surv_rsquared_table, Species[Species == '8'] <- 'Goodenia rosea')

surv_model_list <- list(arcasurvfinalmod, hyglsurvfinalmod, larosurvfinalmod, peaisurvfinalmod, pldesurvfinalmod, trcysurvfinalmod, trorsurvfinalmod, verosurvfinalmod)
effects = lapply(1:length(surv_model_list), function(x) {
  as.data.frame(coef(summary(surv_model_list[[x]]))) %>% mutate(Species=paste0(x))})
surv_effects_table <- do.call("rbind", effects)

surv_effects_table <- cbind(Effect = rownames(surv_effects_table), surv_effects_table)
rownames(surv_effects_table) <- NULL

surv_effects_table$Effect[startsWith(surv_effects_table$Effect, '(Intercept)')] <- 'Intercept'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'CoverShade:std_logp1_totalabund')] <- 'Cover (shade):Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry:CoverShade')] <- 'Dry:Cover (shade)'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet:CoverShade')] <- 'Wet:Cover (shade)'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry:I(std_logp1_totalabund^2)')] <- 'Dry:Neighbour abundance^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet:I(std_logp1_totalabund^2)')] <- 'Wet:Neighbour abundance^2'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry:std_logp1_totalabund')] <- 'Dry:Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet:std_logp1_totalabund')] <- 'Wet:Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'Dodder01')] <- 'Dodder (present)'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_logp1_totalabund')] <- 'Neighbour abundance'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'CoverShade')] <- 'Cover (shade)'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_pH')] <- 'pH'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'std_N')] <- 'Nitrogen'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
surv_effects_table$Effect[startsWith(surv_effects_table$Effect, 'TreatmentWet')] <- 'Wet'

surv_effects_table <- within(surv_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
surv_effects_table <- within(surv_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
surv_effects_table <- within(surv_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
surv_effects_table <- within(surv_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
surv_effects_table <- within(surv_effects_table, Species[Species == '5'] <- 'Plantago debilis')
surv_effects_table <- within(surv_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
surv_effects_table <- within(surv_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
surv_effects_table <- within(surv_effects_table, Species[Species == '8'] <- 'Goodenia rosea')

surv_effects_table <- surv_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')
surv_effects_table <- surv_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
surv_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", surv_effects_table$Estimate, surv_effects_table$SE, surv_effects_table$p_asterisks)

surv_effects_kbl <- surv_effects_table %>% select(Species, Effect, collated)

#Not sure why Dry and Wet are duplicating, extensively troubleshooted and not sure, but grouping by effects works well
surv_effects_kbl <- surv_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)

#Plotting with kableR
surv_effects_kbl %>% kbl(align = 'lcccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.44/0.48"=1, "0.94/0.94"=1, "0.06/0.06"=1, "0.10/0.14"=1, "0.28/0.28"=1, "0.23/0.23"=1, "0.21/0.33"=1, "0.20/0.39"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Survival" = 1, "n=149"=1, "n=100"=1, "n=156"=1, "n=142"=1, "n=65"=1, "n=154"=1, "n=145"=1, "n=134"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

## FECUNDITY #### 
seed_model_list <- list(arcaseedfinalmod, hyglseedfinalmod, laroseedfinalmod, peaiseedfinalmod, pldeseedfinalmod, trcyseedfinalmod, trorseedfinalmod, veroseedfinalmod)
rsquared = lapply(1:length(seed_model_list), function(x) {
  as.data.frame(r.squaredGLMM(seed_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
seed_rsquared_table <- do.call("rbind", rsquared)

#Rename species
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '7'] <- 'Trachymene ornata')
seed_rsquared_table <- within(seed_rsquared_table, Species[Species == '8'] <- 'Goodenia rosea')

seed_model_list <- list(arcaseedfinalmod, hyglseedfinalmod, laroseedfinalmod, peaiseedfinalmod, pldeseedfinalmod, trcyseedfinalmod, trorseedfinalmod, veroseedfinalmod)
#glmmTMB model coefs need to be extracted slightly differently
effects = lapply(1:length(seed_model_list), function(x) {
  as.data.frame(coef(summary(seed_model_list[[x]]))[["cond"]]) %>% mutate(Species=paste0(x))})
seed_effects_table <- do.call("rbind", effects)

seed_effects_table <- cbind(Effect = rownames(seed_effects_table), seed_effects_table)
rownames(seed_effects_table) <- NULL

seed_effects_table$Effect[startsWith(seed_effects_table$Effect, '(Intercept)')] <- 'Intercept'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'CoverShade:std_logp1_totalabund')] <- 'Cover (shade):Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry:CoverShade')] <- 'Dry:Cover (shade)'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet:CoverShade')] <- 'Wet:Cover (shade)'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry:std_logp1_totalabund')] <- 'Dry:Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet:std_logp1_totalabund')] <- 'Wet:Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'Dodder01')] <- 'Dodder (present)'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_logp1_totalabund')] <- 'Neighbour abundance'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_pH')] <- 'pH'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'std_N')] <- 'Nitrogen'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'CoverShade')] <- 'Cover (shade)'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
seed_effects_table$Effect[startsWith(seed_effects_table$Effect, 'TreatmentWet')] <- 'Wet'

seed_effects_table <- within(seed_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
seed_effects_table <- within(seed_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
seed_effects_table <- within(seed_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
seed_effects_table <- within(seed_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
seed_effects_table <- within(seed_effects_table, Species[Species == '5'] <- 'Plantago debilis')
seed_effects_table <- within(seed_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
seed_effects_table <- within(seed_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
seed_effects_table <- within(seed_effects_table, Species[Species == '8'] <- 'Goodenia rosea')

seed_effects_table <- seed_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|z|)')
seed_effects_table <- seed_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                            p_value <0.001~"***",
                                                                            p_value <0.01~"**",
                                                                            p_value <0.05~"*"))
seed_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", seed_effects_table$Estimate, seed_effects_table$SE, seed_effects_table$p_asterisks)

seed_effects_kbl <- seed_effects_table %>% select(Species, Effect, collated)

seed_effects_kbl <- seed_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)

#Plotting with kableR
seed_effects_kbl %>% kbl(align = 'lcccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.46/0.67"=1, "0.26/0.26"=1, "0.23/0.23"=1, "0.46/0.48"=1, "0.49/0.54"=1, "0.21/0.22"=1, "0.32/0.32"=1, "0.16/0.44"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Seed production" = 1, "n=55"=1, "n=41"=1, "n=79"=1, "n=76"=1, "n=35"=1, "n=112"=1, "n=76"=1, "n=93"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

## LAMBDA ####  
#UPDATE COVERSUN TO COVERSHADE*
lambda_model_list <- list(arcalambdafinalmod, hygllambdafinalmod, larolambdafinalmod, peailambdafinalmod, pldelambdafinalmod, trcylambdafinalmod, trorlambdafinalmod, verolambdafinalmod)
rsquared = lapply(1:length(lambda_model_list), function(x) {
  as.data.frame(r.squaredGLMM(lambda_model_list[[x]]))[1,] %>% mutate(Species=paste0(x))})
lambda_rsquared_table <- do.call("rbind", rsquared)

#Rename species
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '1'] <- 'Arctotheca calendula')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '3'] <- 'Lawrencella rosea')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '4'] <- 'Pentameris airoides')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '5'] <- 'Plantago debilis')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '7'] <- 'Trachymene ornata')
lambda_rsquared_table <- within(lambda_rsquared_table, Species[Species == '8'] <- 'Goodenia rosea')

lambda_model_list <- list(arcalambdafinalmod, hygllambdafinalmod, larolambdafinalmod, peailambdafinalmod, pldelambdafinalmod, trcylambdafinalmod, trorlambdafinalmod, verolambdafinalmod)
effects = lapply(1:length(lambda_model_list), function(x) {
  as.data.frame(coef(summary(lambda_model_list[[x]]))) %>% mutate(Species=paste0(x))})
lambda_effects_table <- do.call("rbind", effects)

lambda_effects_table <- cbind(Effect = rownames(lambda_effects_table), lambda_effects_table)
rownames(lambda_effects_table) <- NULL

lambda_effects_table <- within(lambda_effects_table, Species[Species == '1'] <- 'Arctotheca calendula')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '2'] <- 'Hyalosperma glutinosum')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '3'] <- 'Lawrencella rosea')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '4'] <- 'Pentameris airoides')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '5'] <- 'Plantago debilis')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '6'] <- 'Trachymene cyanopetala')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '7'] <- 'Trachymene ornata')
lambda_effects_table <- within(lambda_effects_table, Species[Species == '8'] <- 'Goodenia rosea')

lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, '(Intercept)')] <- 'Intercept'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'CoverShade:Neighbours01Neighbours1')] <- 'Cover (shade):Neighbours (present)'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'Neighbours01Neighbours1')] <- 'Neighbours (present)'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry:CoverShade')] <- 'Dry:Cover (shade)'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet:CoverShade')] <- 'Wet:Cover (shade)'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry:Neighbours01Neighbours1')] <- 'Dry:Neighbours (present)'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet:Neighbours01Neighbours1')] <- 'Wet:Neighbours (present)'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentDry')] <- 'Dry'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'TreatmentWet')] <- 'Wet'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_pH')] <- 'pH'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'std_N')] <- 'Nitrogen'
lambda_effects_table$Effect[startsWith(lambda_effects_table$Effect, 'CoverShade')] <- 'Cover (shade)'


lambda_effects_table <- lambda_effects_table %>% select(Species, Effect, Estimate, 'SE' = 'Std. Error', 'p_value' = 'Pr(>|t|)')
lambda_effects_table <- lambda_effects_table %>% mutate(p_asterisks = case_when(p_value >=0.05~"",
                                                                                p_value <0.001~"***",
                                                                                p_value <0.01~"**",
                                                                                p_value <0.05~"*"))
lambda_effects_table$collated <- sprintf("%1.1f ± %1.1f%s", lambda_effects_table$Estimate, lambda_effects_table$SE, lambda_effects_table$p_asterisks)
lambda_effects_kbl <- lambda_effects_table %>% select(Species, Effect, collated)
lambda_effects_kbl <- lambda_effects_kbl %>% group_by(Species, Effect) %>%
  pivot_wider(names_from = Species, values_from = collated)
#Plotting with kableR
lambda_effects_kbl %>% kbl(align = 'lcccccccc') %>%
  kable_classic(full_width = T, html_font = "Times", font_size = 12) %>%
  add_header_above(c("R^2 m/c" = 1, "0.22/0.22"=1, "0.27/0.49"=1, "0.20/0.60"=1, "0.36/0.52"=1, "0.39/0.81"=1, "0.26/0.75"=1, "0.12/0.38"=1, "0.26/0.49"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  add_header_above(c("Population growth" = 1, "n=48"=1, "n=48"=1, "n=48"=1, "n=48"=1, "n=24"=1, "n=48"=1, "n=48"=1, "n=48"=1), align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"), italic = T, background = "lightgrey") %>%
  row_spec(0, italic = T) %>%
  column_spec(1, italic = F) %>%
  column_spec(1:9, width = 4)

#### Plotting germination in open v shade ####
## Boxplot below - not using
##Plotting raw data as percent_germ with boxplots in open vs shade.
# boxplot(percent_germ ~ Cover, pch=19, ylab=NA, xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, ylim=c(0,1), hygldata)
# stripchart(percent_germ ~ Cover, hygldata, pch = 19, method = "jitter", col=alpha("grey60", 0.4), vertical = TRUE, cex = 2, add = TRUE)
# text(x = 1.5, y = 0.7, "*", cex = 10, col = "red")

##Plotting from model germination, with percent_germ
hyglgermfinalmod <- glmer(cbind(total_germ, total_no_germ) ~ Cover + std_N + std_pH + (1|Site/Plot/rowID), 
                          family = binomial, hygldata)
hygl_germ_pred<-glmm.predict(mod=hyglgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                              se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
#add cover name
hygl_germ_pred$Cover <- c('Shade', 'Sun')

## Plot them
ggplot(hygl_germ_pred, aes(x = Cover, y = y))+
  geom_point(position = position_dodge(0.8), cex=2.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.3), position = position_dodge(0.8), cex=1)+
  ylab("Probability of germination")+
  xlab("Cover")+
  theme_classic()+
  my_theme+
  theme(axis.ticks.x = element_blank())
## Want to add in the raw datapoints underneath, percent_germ
at.x <- seq(1,by=.5,length.out=2) # set here the X-axis positions
ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, cex = 2)+
  geom_point(data=hygl_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.1), cex=1)+
 # geom_text(data=growth_emmeans, aes(x=Cover, y=10, label=c('a')), size =6)+
  ylab("Probability of germination")+
  #ylim(0, 1)+
  theme_classic()+
  theme(axis.text=element_text(size=16), 
        axis.title=element_text(size=16),
       aspect.ratio=1)
       # axis.title.x=element_blank(),
        #axis.text.x=element_blank())

## survival is going to be plotted for ambient plots with average
# amount of neighbours, absence of Dodder, avreage pH and N
#library(scales)
#need this to make y axis go up by 0.1 for survival


## trying for survival
##compare to boxplot
boxplot(surv_to_produce_seeds ~ Cover, pch=19, ylab=NA, xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, ylim=c(0,1.1), trcydata)
stripchart(surv_to_produce_seeds ~ Cover, trcydata, pch = 19, method = "jitter", col=alpha("grey60", 0.2), vertical = TRUE, cex = 2, add = TRUE)
#need to jitter vertically a bit

trcysurvfinalmod <- glmer(surv_to_produce_seeds ~ Treatment + Cover + std_pH + std_N + std_logp1_totalabund + Dodder01 +
                            Treatment:Cover + Cover:std_logp1_totalabund + (1|Site/Plot),
                          family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), trcydata)
summary(trcysurvfinalmod)
#holding everything at its mean
## you can make multiple predictions at once - here I made the predictions for sun =1 and sun=0
#Ambient watering, no dodder, mean pH and mean N. Estimated for sun and shade
trcy_cover_pred<-glmm.predict(mod=trcysurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1*0)*0), 
                           se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
#add cover name
trcy_cover_pred$Cover <- c('open', 'shade')

## Plot them
ggplot(trcy_cover_pred, aes(x = Cover, y = y))+
  geom_point(position = position_dodge(0.8), cex=2.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.3), position = position_dodge(0.8), cex=1)+
  ylab("Probability of survival")+
  xlab("Cover")+
  theme_classic()+
  my_theme+
  theme(axis.ticks.x = element_blank())
#could plot survival datapoints onto this as well

#seed proudction
#type='response'


#### Adding correlation tests between vital rates ####
#### Build dataframe with plot-level vital rates ####
#Need a dataframe with one value for each rate
#vitaldata (germ and survival)
#seedmodeldata (produced at least one seed)
#popdata (lambda) --- this dataframe has data for all rates split by nbhs

#popgrowthrate data has germination by plot, surv and seed production by neighbours by plot
# so does popdata, with lambda
### need to go back to the code in data prep, do same thing
# but for all data, not split by neighbours

meanvitalrates <- popgrowthratedata
#without pole
meanvitalrates <- meanvitalrates %>% filter(!(Species=="POLE"))
#already have plot_germ, need plot_surv, plot_fec and plot_lambda

##### plot_surv
plot_survdata <- vitaldata %>% filter(!(is.na(surv_to_produce_seeds)))
plot_surv_counts <- plot_survdata %>% group_by(Species, Site, Plot, surv_to_produce_seeds) %>% count()
Species <- rep(c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO"), each = 3, times = 8)
Site <- rep(c("1", "2", "3", "4", "5", "6", "7", "8"), each = 9, times = 3)
Plot <- rep(c("A", "B", "C"), times = 72)
plot_surv_prop <- cbind(Species, Site, Plot)
plot_surv_prop <- data.frame(plot_surv_prop)
plot_surv_prop <- plot_surv_prop %>% unite("idforjoining", Plot:Site, sep = ":", remove = "false")
plot_surv_counts$surv_to_produce_seeds <- as.factor(plot_surv_counts$surv_to_produce_seeds)
plot_surv_prop <- left_join(plot_surv_prop, plot_surv_counts)
plot_surv_prop <- within(plot_surv_prop, n[is.na(n)] <- 0)

#
plot_surv <- plot_surv_prop %>% group_by(Species, Site, Plot)  %>% 
  summarise(plot_surv = ifelse(n[surv_to_produce_seeds==1]==0 & n[surv_to_produce_seeds==0]==0, NA, n[surv_to_produce_seeds==1]/(n[surv_to_produce_seeds==1] + n[surv_to_produce_seeds==0])),
            total_n = sum(n))
## Merge with plot-level data
plot_surv <- plot_surv %>% select(Species, Site, Plot, plot_surv)
meanvitalrates <- left_join(meanvitalrates, plot_surv)

#Calculating mean viable seed production of logged seed values (exponentiated)
allplotfec <- viable_plot %>% mutate(log_seeds = log(No_viable_seeds_grouped)) %>%
  group_by(Species, Site, Plot) %>%
  summarise(plot_fec = exp(mean(log_seeds)))
##this is only where there were seeds produced. Need to factor in 0s?
#I don't think so
#Merge with plot-level data
meanvitalrates <- left_join(meanvitalrates, allplotfec)

#### Need plot_lambda (not in presence/absence of neighbours)
### Assigning all survival and fecundity NA value to 0 for lambda calculations
#(no plants germinated, none survived, none produced seeds)
one_poplongdata <- meanvitalrates
one_poplongdata <- within(one_poplongdata, plot_surv[is.na(plot_surv)] <- 0)
one_poplongdata <- within(one_poplongdata, plot_fec[is.na(plot_fec)] <- 0)

#Calculate population growth rates
# Per capita growth rate of a given population  i = seed survival*(1-germination)+number of viable seeds produced per germinant*germination
# prob of germination / germination fraction. # germinated / total germination, currently I have this as percent_germ which is a proportion (not percentage, despite the name)
one_poplongdata <- one_poplongdata %>% mutate(plot_lambda = seed_survival*(1-plot_germ)+plot_fec*plot_surv*plot_germ)
#Still get NAs where seed wasn't sown, working well

#Merge this back with plot-levle data
one_poplongdata <- one_poplongdata %>% select(Species, Site, Plot, plot_lambda)
meanvitalrates <- left_join(meanvitalrates, one_poplongdata)

#Simplify main dataframe
plot_vitalrates <- meanvitalrates %>% select(Species, Site, Plot, seed_survival, plot_germ, plot_surv, plot_fec, plot_lambda)

# I think we need to keep our zero values for fecundity and survival, instead of NAs
#Yes! Need zero values, these are already factored into popdata too
plot_zero <- plot_vitalrates
plot_zero <- within(plot_zero, plot_surv[is.na(plot_surv)] <- 0)
plot_zero <- within(plot_zero, plot_fec[is.na(plot_fec)] <- 0)

## Change the column names, for plotting
colnames(plot_zero) <- c("Species", "Site", "Plot", "seed_survival", "emergence", "survival", "seed production", "population growth rate")

#Split the dataframe by species
plotarca <- plot_zero %>% filter(Species == "ARCA")
plothygl <- plot_zero %>% filter(Species == "HYGL")
plotlaro <- plot_zero %>% filter(Species == "LARO")
plotpeai <- plot_zero %>% filter(Species == "PEAI")
plotplde <- plot_zero %>% filter(Species == "PLDE")
plottrcy <- plot_zero %>% filter(Species == "TRCY")
plottror <- plot_zero %>% filter(Species == "TROR")
plotvero <- plot_zero %>% filter(Species == "VERO")

#Need just the rates for plotting
corrarca <- plotarca %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrhygl <- plothygl %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrlaro <- plotlaro %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrpeai <- plotpeai %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrplde <- plotplde %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrtrcy <- plottrcy %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrtror <- plottror %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)
corrvero <- plotvero %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

#overall (all species)
justrates <- plot_zero %>% ungroup() %>% select(emergence, survival, `seed production`, `population growth rate`)

#### Correlation matrices ####

#Two different correlation plots, all per species
# regardless of environment: all data
# unmanipulated plots only: unthinned, ambient

#### ALL PLOTS ##
library(RColorBrewer)

#overall
#Need use= 'complete' to omit NAs
M = cor(justrates, use = 'complete')
#diag = FALSE removed the 1:1 correlations
corrplot(M, method = 'square', type = 'upper', diag = FALSE)
## order = 'original' keeps the order in the data frame
corrplot(M, method = 'circle', type = 'upper', order = 'original')
corrplot(M, method = 'number', type = 'upper')
corrplot.mixed(M, type = 'lower')
corrplot(M, method = 'ellipse', type = 'upper', addCoef.col='black', diag = FALSE)

#my favourite combination:
corrplot(M, method = 'square', type = 'lower', addCoef.col='grey30', 
         diag = FALSE, order='original', tl.col="black", cl.pos='n',
         p.mat = p.mat,
         insig = "label_sig",
         pch.col = "red",
         sig.level = c(.001, .01, .05))

corrplot(M, method = 'square', type = 'lower', addCoef.col='grey30', 
         diag = FALSE, order='original', tl.col="black", cl.pos='n',
         p.mat = p.mat,
         pch.col = "red",
         sig.level = c(.001, .01, .05))

## leave blank on non-significant coefficient
## add significant correlation coefficients
corrplot(M, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='grey30')

## leave blank on non-significant coefficient
## add all correlation coefficients
corrplot(M, p.mat = p.mat, method = 'circle', type = 'lower', insig='blank',
         order = 'original', diag = FALSE)$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))

# Visualize confidence interval
corrplot(M, lowCI = testRes$lowCI, uppCI = testRes$uppCI, order = 'hclust',
         tl.pos = 'd', rect.col = 'navy', plotC = 'rect', cl.pos = 'n')

#labels = NULL
#cl.pos='n' means don't draw legend
trace(corrplot, edit=TRUE)
#M<-cor(mtcars)
#res1 <- cor.mtest(mtcars, conf.level = .95)
corrplot(arca_all,
         method="square",
         type="lower",
         p.mat = p.mat,
         insig = "label_sig",
         sig.level = c(.001, .01, .05),
         pch.cex = 0.8,
         pch.col = "red",
         tl.col="black",
         tl.cex=1,
         addCoef.col = "black",
         tl.pos="n",
         diag = FALSE)
#Look into adding p values to the plot **
#How does this work statistically? E.g. accounting for my design? Blocks?
#Can leave non-signficant correlation coefficients blank

#https://github.com/caijun/ggcorrplot2
# Use corr.test() from psych package to calculate the correlation matrix and 
# corresponding p value matrix without adjustment.
library(psych)
ct <- corr.test(M, adjust = "none")
corr <- ct$r
p.mat <- ct$p
#(sig.lvl = 0.05) are indicated by X by default.
ggcorrplot.mixed(M, upper = "ellipse", lower = "number", p.mat = p.mat)

#https://stackoverflow.com/questions/63227830/r-corrplot-plot-correlation-coefficients-along-with-significance-stars
##Above article shows how to display both significance stars and R valus

#Change the place_points function within the corrplot function. To do so, run:
# 
# trace(corrplot, edit=TRUE)
# Then replace on line 443
# 
# place_points = function(sig.locs, point) {
#   text(pos.pNew[, 1][sig.locs], pos.pNew[, 2][sig.locs], 
#        labels = point, col = pch.col, cex = pch.cex, 
#        lwd = 2)
#   with:
#     
#     # adjust text(X,Y ...) according to your needs, here +0.25 is added to the Y-position    
#     place_points = function(sig.locs, point) {
#       text(pos.pNew[, 1][sig.locs], (pos.pNew[, 2][sig.locs])+0.25, 
#            labels = point, col = pch.col, cex = pch.cex, 
#            lwd = 2)
### NOTE All that does is add +0.25 to height


p.mat <- cor.mtest(justrates, conf.level = .95)

corrplot(M, method = 'square', p.mat = p.mat, type = 'upper', diag = FALSE)
corrplot(M, method = 'square', p.mat = p.mat, type = 'lower', diag = FALSE, insig='p-value')
#insig='p-value' gives the p-value of the insignificant ones
corrplot(M, method = 'circle', p.mat = p.mat, sig.level = c(0.001, 0.01, 0.05), pch.col = 'grey20', pch.cex=0.9, insig = 'label_sig', type = 'upper', diag = FALSE)

corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=10, name="PRGn"))
#PuOr, #PRGn, #PiYG, #BrBG, #RdBu
#n = 4

#Fecundity and lambda are highly correlated.... 0.68.
modfeclambda <- lmer(plot_lambda ~ plot_fec + (1|Site/Plot), plot_zero)
summary(modfeclambda)
r.squaredGLMM(modfeclambda)
#0.39 R squared though. hm.
#0.40 when we use 0 values instead of NAs for surv and fecundity.

#par(mfrow=c(2,4))
#rows, columns

library(ggcorrplot2)

dev.off()
pdf("Output/Figures/panel_all_corr.pdf", width = 21, height = 12)
par(mfrow=c(2,4), oma =c(0,14,14,0))

arca_all = cor(corrarca, use = 'complete')
corrplot(arca_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#hygl
hygl_all = cor(corrhygl, use = 'complete')
corrplot(hygl_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#laro
laro_all = cor(corrlaro, use = 'complete')
corrplot(laro_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#peai
peai_all = cor(corrpeai, use = 'complete')
corrplot(peai_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#plde
plde_all = cor(corrplde, use = 'complete')
corrplot(plde_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#trcy
trcy_all = cor(corrtrcy, use = 'complete')
corrplot(trcy_all, method = 'square', type = 'lower', diag = FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#tror
tror_all = cor(corrtror, use = 'complete')
corrplot(tror_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
#vero
vero_all = cor(corrvero, use = 'complete')
corrplot(vero_all, method = 'square', type = 'lower', diag=FALSE, p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='n', cl.pos='n', number.cex=3, tl.cex=2)
dev.off()

### tryinga different way
dev.off()
pdf("Output/Figures/panel_all_corr2.pdf", width = 21, height = 12)
par(mfrow=c(2,4), oma =c(1,1,1,1))

arca_all = cor(corrarca, use = 'complete')
corrplot(arca_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
t#hygl
hygl_all = cor(corrhygl, use = 'complete')
corrplot(hygl_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#laro
laro_all = cor(corrlaro, use = 'complete')
corrplot(laro_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#peai
peai_all = cor(corrpeai, use = 'complete')
corrplot(peai_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#plde
plde_all = cor(corrplde, use = 'complete')
corrplot(plde_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#trcy
trcy_all = cor(corrtrcy, use = 'complete')
corrplot(trcy_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#tror
tror_all = cor(corrtror, use = 'complete')
corrplot(tror_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
#vero
vero_all = cor(corrvero, use = 'complete')
corrplot(vero_all, method = 'square', type = 'lower', p.mat = p.mat, insig = 'blank', sig.level = 0.05, order = 'original', tl.col="black", addCoef.col='black', tl.pos='d', cl.pos='n', number.cex=3, tl.cex=0.01)
dev.off()



##Creating legend to add in manually
corrplot(vero_all, method = 'circle', type = 'lower', tl.col="black")

ggplot(corrtror, aes(x=plot_fec, y=plot_lambda))+
  geom_point()+
  geom_smooth(method='lm')+
  theme_classic()
#check this - highlight correlated except for the one outlier with high fecundity
modfeclambda <- lm(plot_lambda ~ plot_fec, corrtror)
summary(modfeclambda)

### Plotting together
p1 <- corrplot(M, method = 'circle', type = 'upper', addCoef.col='black', diag = FALSE, order='original', tl.col="black")
p2 <- corrplot(M, method = 'number', type = 'upper', diag = FALSE, order='original', tl.col="black")

p1 <- ggcorrplot(M, method = 'circle', type = 'upper')
p2 <- ggcorrplot(M, method = 'number', type = 'upper')


library(cowplot)
prow <- plot_grid(p1 + theme(legend.position = "none"),
                  p2 + theme(legend.position = "none"),
                  rel_widths = c(1, 1), nrow = 1, align = 'hv',
                  labels = c("(a)", "(b)"), label_x = 0, label_y = 1)

plot_grid(p1, p2)
# Extract the legend from the first corrgram
#legend <- get_legend(p1)
# Add the legend to the bottom of the plot row we made earlier.
p <- cowplot::plot_grid(prow, ncol = 1)
p

ggcorrplot(corr, type = "lower")

### UNMANIPULATED PLOTS ONLY ##

#Filtering to plots with neighbours and ambient watering
nbh_justrates <- popdata %>% filter(Neighbours01 == "Neighbours1") %>% ungroup() %>% 
  select(plot_germ, plot_survival, plot_fecundity, lambda)
#still need to do ambient watering **
test <- vitaldata %>% filter()
M_nbh = cor(nbh_justrates, use = 'complete')
corrplot(M_nbh, method = 'circle', type = 'upper')
corrplot(M_nbh, method = 'number', type = 'upper')



