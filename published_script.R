#everything needed to make figures for JoE paper
# Resubmission 2024
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
library(scales)
library(cowplot)


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
  geom_jitter(alpha = 0.3, width=0.1, height=0.2)+
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



#### Figure 3 - rates ~ cover ####
dev.off()
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(2,10,1,1))
### ARCA ####
coef(arcagermfinalmod)
arca_germ_pred<-glmm.predict(mod=arcagermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
arca_germ_pred$Cover <- c('Shade', 'Sun')
a <- ggplot()+
  geom_jitter(data=arcadata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=arca_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(arcasurvfinalmod)
arca_surv_pred<-glmm.predict(mod=arcasurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
arca_surv_pred$Cover <- c('Shade', 'Sun')
b <- ggplot()+
  geom_jitter(data=arcadata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=arca_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
#arca
#CAN'T go below 0, the estimate from model. troubleshooting below
# coef(arcaseedfinalmod)
# arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
#                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=TRUE)
# arca_seed_pred$Cover <- c('Shade', 'Sun')
# arca_seed_pred$lower <- c(0.165081, 0.3698253)
# #Get this lower value from log_link=TRUE, but it does fix the 0 issue
# #or I could manually set it to zero?
# ggplot()+
#   geom_jitter(data=seedarca, aes(x = Cover, y = log(No_viable_seeds_grouped+1)), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
#   geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
#   geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
#   theme_classic()+
#   theme(axis.text=element_text(size=24), 
#         axis.title=element_text(size=24),
#         axis.title.y=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank())
#Using log10 y axis but still below 0 so adding 1
#but adding 1 changes the error bars' position compared to plot above
# coef(arcaseedfinalmod)
# arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
#                              se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
# arca_seed_pred$Cover <- c('Shade', 'Sun')
# arca_seed_pred$lower <- arca_seed_pred$lower+1
# arca_seed_pred$upper <- arca_seed_pred$upper+1
# #add one to the mean??
# ggplot()+
#   geom_jitter(data=seedarca, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
#   geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
#   geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
#   scale_y_continuous(trans="log10")+
#   theme_classic()+
#   theme(axis.text=element_text(size=24),
#         axis.title=element_text(size=24),
#         axis.title.y=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank())

#plotting on a log scale, so need seeds+1 and +1 to negative lower limits
#exponentiating after logging, log_link=TRUE IS exponentiating it for us
coef(arcaseedfinalmod)
arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
arca_seed_pred$Cover <- c('Shade', 'Sun')
#added one to the negative lower limit
arca_seed_pred$lower <- c(1.165081, 1.447482)
c <- ggplot()+
  geom_jitter(data=seedarca, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

#in the absence of neighbours
coef(arcalambdafinalmod)
arca_lambda_pred<-glmm.predict(mod=arcalambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
arca_lambda_pred$Cover <- c('Shade', 'Sun')
d <- ggplot()+
  geom_jitter(data=lambdaarca, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=arca_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### hygl ####
coef(hyglgermfinalmod)
hygl_germ_pred<-glmm.predict(mod=hyglgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_germ_pred$Cover <- c('Shade', 'Sun')
e <- ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=hygl_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hyglsurvfinalmod)
hygl_surv_pred<-glmm.predict(mod=hyglsurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), 0*0, 0*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_surv_pred$Cover <- c('Shade', 'Sun')
f <- ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=hygl_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hyglseedfinalmod)
hygl_seed_pred<-glmm.predict(mod=hyglseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
hygl_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
g <- ggplot()+
  geom_jitter(data=seedhygl, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=hygl_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hygllambdafinalmod)
hygl_lambda_pred<-glmm.predict(mod=hygllambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
hygl_lambda_pred$Cover <- c('Shade', 'Sun')
#negative values here are fine! Population growth rate
h <- ggplot()+
  geom_jitter(data=lambdahygl, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=hygl_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### laro ####
coef(larogermfinalmod)
laro_germ_pred<-glmm.predict(mod=larogermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_germ_pred$Cover <- c('Shade', 'Sun')
i <- ggplot()+
  geom_jitter(data=larodata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=laro_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(larosurvfinalmod)
laro_surv_pred<-glmm.predict(mod=larosurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_surv_pred$Cover <- c('Shade', 'Sun')
j <- ggplot()+
  geom_jitter(data=larodata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=laro_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(laroseedfinalmod)
laro_seed_pred<-glmm.predict(mod=laroseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
laro_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
k <- ggplot()+
  geom_jitter(data=seedlaro, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=laro_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(larolambdafinalmod)
laro_lambda_pred<-glmm.predict(mod=larolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
laro_lambda_pred$Cover <- c('Shade', 'Sun')
l <- ggplot()+
  geom_jitter(data=lambdalaro, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=laro_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### peai ####
coef(peaigermfinalmod)
peai_germ_pred<-glmm.predict(mod=peaigermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_germ_pred$Cover <- c('Shade', 'Sun')
m <- ggplot()+
  geom_jitter(data=peaidata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=peai_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peaisurvfinalmod)
peai_surv_pred<-glmm.predict(mod=peaisurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_surv_pred$Cover <- c('Shade', 'Sun')
n <- ggplot()+
  geom_jitter(data=peaidata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=peai_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peaiseedfinalmod)
peai_seed_pred<-glmm.predict(mod=peaiseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
peai_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
o <- ggplot()+
  geom_jitter(data=seedpeai, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=peai_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  annotate("text", x = 1.5, y= 120, label = "*", family = "", size = 20, colour = "red")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
coef(peailambdafinalmod)
peai_lambda_pred<-glmm.predict(mod=peailambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
peai_lambda_pred$Cover <- c('Shade', 'Sun')
p <- ggplot()+
  geom_jitter(data=lambdapeai, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=peai_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### plde ####
coef(pldegermfinalmod)
plde_germ_pred<-glmm.predict(mod=pldegermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_germ_pred$Cover <- c('Shade', 'Sun')
q <- ggplot()+
  geom_jitter(data=pldedata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=plde_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldesurvfinalmod)
plde_surv_pred<-glmm.predict(mod=pldesurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_surv_pred$Cover <- c('Shade', 'Sun')
r <- ggplot()+
  geom_jitter(data=pldedata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=plde_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldeseedfinalmod)
plde_seed_pred<-glmm.predict(mod=pldeseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plde_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
s <- ggplot()+
  geom_jitter(data=seedplde, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=plde_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  annotate("text", x = 1.5, y= 20, label = "*", family = "", size = 20, colour = "red")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldelambdafinalmod)
plde_lambda_pred<-glmm.predict(mod=pldelambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plde_lambda_pred$Cover <- c('Shade', 'Sun')
t <- ggplot()+
  geom_jitter(data=lambdaplde, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=plde_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### trcy ####
coef(trcygermfinalmod)
trcy_germ_pred<-glmm.predict(mod=trcygermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_germ_pred$Cover <- c('Shade', 'Sun')
u <- ggplot()+
  geom_jitter(data=trcydata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=trcy_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcysurvfinalmod)
trcy_surv_pred<-glmm.predict(mod=trcysurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_surv_pred$Cover <- c('Shade', 'Sun')
v <- ggplot()+
  geom_jitter(data=trcydata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=trcy_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcyseedfinalmod)
trcy_seed_pred<-glmm.predict(mod=trcyseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, c(1,0)*0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
trcy_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
w <- ggplot()+
  geom_jitter(data=seedtrcy, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=trcy_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcylambdafinalmod)
trcy_lambda_pred<-glmm.predict(mod=trcylambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0*0, 0*0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
trcy_lambda_pred$Cover <- c('Shade', 'Sun')
x <- ggplot()+
  geom_jitter(data=lambdatrcy, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=trcy_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### tror ####
coef(trorgermfinalmod)
tror_germ_pred<-glmm.predict(mod=trorgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_germ_pred$Cover <- c('Shade', 'Sun')
y <- ggplot()+
  geom_jitter(data=trordata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=tror_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorsurvfinalmod)
tror_surv_pred<-glmm.predict(mod=trorsurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0)),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_surv_pred$Cover <- c('Shade', 'Sun')
z <- ggplot()+
  geom_jitter(data=trordata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=tror_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorseedfinalmod)
tror_seed_pred<-glmm.predict(mod=trorseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0)),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
tror_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
ab <- ggplot()+
  geom_jitter(data=seedtror, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=tror_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorlambdafinalmod)
tror_lambda_pred<-glmm.predict(mod=trorlambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
tror_lambda_pred$Cover <- c('Shade', 'Sun')
bc <- ggplot()+
  geom_jitter(data=lambdatror, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=tror_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### vero ####
coef(verogermfinalmod)
vero_germ_pred<-glmm.predict(mod=verogermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_germ_pred$Cover <- c('Shade', 'Sun')
cd <- ggplot()+
  geom_jitter(data=verodata, aes(x = Cover, y = percent_germ), alpha = 0.2, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=vero_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  geom_text(aes(label = "*"), size = 20, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(verosurvfinalmod)
vero_surv_pred<-glmm.predict(mod=verosurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_surv_pred$Cover <- c('Shade', 'Sun')
de <- ggplot()+
  geom_jitter(data=verodata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.2, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=vero_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(veroseedfinalmod)
vero_seed_pred<-glmm.predict(mod=veroseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
vero_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
ef <- ggplot()+
  geom_jitter(data=seedvero, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.3, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=vero_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(verolambdafinalmod)
vero_lambda_pred<-glmm.predict(mod=verolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
vero_lambda_pred$Cover <- c('Shade', 'Sun')
fg <- ggplot()+
  geom_jitter(data=lambdavero, aes(x = Cover, y = log_lambda), alpha = 0.3, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=vero_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())
#### plotting cover_panel from above ####
pdf("Output/Figures/panel_cover_2024.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(2,10,1,1))

plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,ab,bc,cd,de,ef,fg, align="hv", ncol=4)
dev.off()
#hjust=-5.5, 

pdf("Output/Figures/panel_cover_labels.pdf", width=21, height=21)
plot_row <- plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,ab,bc,cd)
title <- ggdraw() + 
  draw_label(
    "Probability of emergence",
    # fontface = 'bold',
    x = 0,
    hjust = 0,
    angle = 90
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(
  title, plot_row,
  ncol = 4, align="hv"
)
dev.off()
### manually adding axis and species labels in Adobe Illustrator
