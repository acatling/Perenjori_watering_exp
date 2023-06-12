# Alexandra Catling
# Code needed to generate figures and tables in the paper entitled:
# Variation in demographic responses to competition and abiotic conditions in an annual plant community
# By Catling, Mayfield and Dwyer
# 2023

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
library(dplyr)
library(sf)
library(ggfortify)
library(exactLTRE)
library(popbio)

#### Running below models needed for figures and tables ####
## Germination ####
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

## Survival ####
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
#problem here - only using this species for germination data, low germination
# polesurvmod1 <- glmer(surv_to_produce_seeds ~ Treatment + std_PC1 + std_PC2 + std_logp1_totalabund + Dodder01 +
#                         Treatment:std_PC1 + Treatment:std_logp1_totalabund + std_PC1:std_logp1_totalabund + (1|Site/Plot), 
#                       family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), poledata)

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

## Fecundity ####
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

## Lambda ####
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

#### FIGURE 1 - Map study site ####
# Figure created in adobe illustrator using following figure graphics

#Read in the SA2 shapefile downloaded from the ABS
#Data from ABS localities
# https://data.gov.au/dataset/ds-wa-d2dc22c6-0840-448c-819f-b6fb21411517/details?q=

ausplotdata <- read_sf("Data_and_code_from_others/SA2_2016_AUST.shp")
#filter the Australian SA2 shapefile for only WA
waplotdata <- ausplotdata %>% filter(STE_NAME16 == "Western Australia")
#import a shapefile of state boundaries
aus_state_data <- read_sf("Data_and_code_from_others/STE_2016_AUST.shp")
#make map of just state boundaries
ggplot()+
  geom_sf(data = aus_state_data, fill = "white")+
  theme_classic()
#make a new dataset with West Perenjori Nature Reserve coordinations
wa_cities2 <- tribble(
  ~city, ~lat, ~long, 
  "Perenjori", -29.443172, 116.288301)
#convert columns to geometry column with the st_as_sf() function. Google Maps uses the coordinate reference system 4326 (the GPS system).
wa_cities_geometry2 <- wa_cities2 %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326)

#Using ggrepel package to try and offset city labels and give them points
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

#### FIGURE 2 - Demographic responses to PC1 (shade, litter, soil) ####
dev.off()
pdf("Output/Figures/panel_PC1_int.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(12, 20, 5, 1), mar =c(2,10,1,1))
#ARCA
x_to_plot<-seq.func(arcadata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcagermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival -signif. int.
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01==1, "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcasurvfinalmod
x_to_plot_no_nbh <- seq.func(arcadata$std_PC1[arcadata$Neighbours01=='0'])
x_to_plot_nbh <- seq.func(arcadata$std_PC1[arcadata$Neighbours01=='1'])
#Without nbhs - blue
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, -0.79, 0, 0*x_to_plot, 0*x_to_plot, x_to_plot*-0.79), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#With nbhs - red
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0.97*x_to_plot, 0*x_to_plot, x_to_plot*0.97), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Seed production
x_to_plot<-seq.func(seedarca$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedarca)
model <- arcaseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdaarca$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaarca)
model <- arcalambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#HYGL
x_to_plot<-seq.func(hygldata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglgermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedhygl$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedhygl)
model <- hyglseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdahygl$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdahygl)
model <- hygllambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#LARO
x_to_plot<-seq.func(larodata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, larodata)
model <- larogermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, larodata)
model <- larosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production - signif. int.
x_to_plot<-seq.func(seedlaro$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha(ifelse(Neighbours01==1, "#CC79A7", "#0072B2"), 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedlaro)
model <- laroseedfinalmod
x_to_plot_no_nbh <- seq.func(seedlaro$std_PC1[seedlaro$Neighbours01=='0'])
x_to_plot_nbh <- seq.func(seedlaro$std_PC1[seedlaro$Neighbours01=='1'])
#Without nbhs - blue
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, -0.79, 0, x_to_plot*-0.79), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#With nbhs - red
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0.97, 0, x_to_plot*0.97), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Population growth - signif. int.
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdalaro)
model <- larolambdafinalmod
x_to_plot_no_nbh <- seq.func(lambdalaro$std_PC1[lambdalaro$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdalaro$std_PC1[lambdalaro$Neighbours01=='Neighbours1'])
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#PEAI
x_to_plot<-seq.func(peaidata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, peaidata)
model <- peaigermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival -signif. int.
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01==1, "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, peaidata)
model <- peaisurvfinalmod
x_to_plot_no_nbh <- seq.func(peaidata$std_PC1[peaidata$Neighbours01=='0'])
x_to_plot_nbh <- seq.func(peaidata$std_PC1[peaidata$Neighbours01=='1'])
#Without nbhs - blue
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, -0.79, 0, x_to_plot*-0.79), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#With nbhs - red
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0.97, 0, x_to_plot*0.97), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Seed production
x_to_plot<-seq.func(seedpeai$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedpeai)
model <- peaiseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 1.35,y = 220,"*", cex = 10, col = "red")
#Population growth - signif. int.
x_to_plot_no_nbh <- seq.func(lambdapeai$std_PC1[lambdapeai$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdapeai$std_PC1[lambdapeai$Neighbours01=='Neighbours1'])
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdapeai)
model <- peailambdafinalmod
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot^2, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, x_to_plot^2, 0*x_to_plot, 0*x_to_plot, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#PLDE
x_to_plot<-seq.func(pldedata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, pldedata)
model <- pldegermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, pldedata)
model <- pldesurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot^2, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedplde$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedplde)
model <- pldeseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdaplde$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaplde)
model <- pldelambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#TRCY
x_to_plot<-seq.func(trcydata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcygermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcysurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0^2, 0*x_to_plot, 0*x_to_plot, 0*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedtrcy$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedtrcy)
model <- trcyseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 1.35,y = 42,"*", cex = 10, col = "red")
#Population growth - signif. int.
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatrcy)
model <- trcylambdafinalmod
x_to_plot_no_nbh <- seq.func(lambdatrcy$std_PC1[lambdatrcy$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdatrcy$std_PC1[lambdatrcy$Neighbours01=='Neighbours1'])
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0*0, 0*0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, 0*1, 0*1, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#TROR
x_to_plot<-seq.func(trordata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorgermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedtror$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedtror)
model <- trorseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdatror$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatror)
model <- trorlambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#VERO
x_to_plot<-seq.func(verodata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verogermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedvero$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedvero)
model <- veroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot^2), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth - signif. int.
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdavero)
model <- verolambdafinalmod
x_to_plot_no_nbh <- seq.func(lambdavero$std_PC1[lambdavero$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdavero$std_PC1[lambdavero$Neighbours01=='Neighbours1'])
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, 0*x_to_plot, 0*x_to_plot, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

###Overall text
##x labels
mtext("PC1", adj = 0.15, side = 1, line = 3, cex = 3,outer = TRUE)
mtext("PC1", adj = 0.4, side = 1, line = 3, cex = 3,outer = TRUE)
mtext("PC1", adj = 0.66, side = 1, line = 3, cex = 3, outer = TRUE)
mtext("PC1", adj = 0.93, side = 1, line = 3, cex = 3, outer = TRUE)
##y labels
mtext("Probability of emergence", side = 2, cex = 3, outer=TRUE, line=-5)
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line=-40)
mtext("Number of viable seeds produced (+1)", side = 2, cex = 3, outer=TRUE, line=-74)
mtext("Population growth rate (log)", side = 2, cex = 3, outer=TRUE, line=-108)
##main labels
mtext("Emergence", outer=TRUE, adj=0.1,side = 3, cex = 3)
mtext("Survival", outer=TRUE, adj=0.39,side = 3, cex = 3)
mtext("Seed production", outer=TRUE, adj = 0.68, side = 3, cex = 3)
mtext("Population growth", outer=TRUE, adj=1, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 4, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 27, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 34, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 34, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 50, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 56, side = 3, cex = 2.5, outer = TRUE)
reset()
legend("bottom", title=NULL, horiz=T, legend=c("No neighbours", "Neighbours"),
       col=c("#0072B2", "#CC79A7"), pch=19, cex=3, bty="n")
dev.off()

#### FIGURE 3 - Demographic responses to neighbours ####
dev.off()
pdf("Output/Figures/panel_NA.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(5, 20, 5, 1), mar =c(3,10,1,1))
#ARCA
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab = NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, arcadata)
model <- arcasurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedarca$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedarca)
model <- arcaseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab=NA, xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdaarca)
stripchart(log_lambda ~ Neighbours01, lambdaarca, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex = 2, add = TRUE)
text(x = 1.5, y = 2.55, "*", cex = 10, col = "red")

#Adding HYGL
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, hygldata)
model <- hyglsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedhygl$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedhygl)
model <- hyglseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdahygl)
stripchart(log_lambda ~ Neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding LARO
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, larodata)
model <- larosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedlaro$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedlaro)
model <- laroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdalaro)
stripchart(log_lambda ~ Neighbours01, lambdalaro, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding PEAI
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, peaidata)
model <- peaisurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedpeai$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedpeai)
model <- peaiseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdapeai)
stripchart(log_lambda ~ Neighbours01, lambdapeai, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding PLDE
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, pldedata)
model <- pldesurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0^2, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedplde$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedplde)
model <- pldeseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 2.3, y = 50, "*", cex = 10, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdaplde)
stripchart(log_lambda ~ Neighbours01, lambdaplde, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)
text(x = 1.5, y = 2, "*", cex = 10, col = "red")

#Adding TRCY
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, trcydata)
model <- trcysurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, x_to_plot^2, 0*0, 0*0, 0*x_to_plot^2, 0*x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedtrcy$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedtrcy)
model <- trcyseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdatrcy)
stripchart(log_lambda ~ Neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding TROR
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, trordata)
model <- trorsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedtror$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedtror)
model <- trorseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdatror)
stripchart(log_lambda ~ Neighbours01, lambdatror, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding VERO
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, verodata)
model <- verosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedvero$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedvero)
model <- veroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab="", names = c("", ""), col = "white", cex.axis = 2.5, lambdavero)
stripchart(log_lambda ~ Neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)
mtext(side=1, "Absent", adj=0.1, line=2.5, cex =2.5)
mtext(side=1, "Present", adj=0.9, line=2.5, cex =2.5)

###Overall text
##x labels
mtext("Neighbour abundance", adj = 0.09, side = 1, line=3, cex = 3,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, line=3, cex = 3, outer = TRUE)
mtext("Neighbour presence", adj = 0.98, side = 1, line=3, cex = 3, outer = TRUE)
##y labels
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line=-6)
mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 3, outer=TRUE, line=-50.5)
mtext("Population growth rate (log)", side = 2, cex = 3, outer=TRUE, line=-97.5)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 3)
mtext("Seed production", outer=TRUE, adj = 0.55, side = 3, cex = 3)
mtext("Population growth", outer=TRUE, adj=0.96, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 29, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 37, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 35, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 52, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 59, side = 3, cex = 2.5, outer = TRUE)

dev.off()
#### SUPP. FIGURE S3 - PCA abiotic variables ####
#Selecting variables of interest for PCA
soilpcadata <- soildata %>% select(Site, Plot, pH, log_NH4N, log_NO3N, log_P, log_K)
abioticpcadata <- merge(soilpcadata, canopydatatrim)
abioticpcadata <- abioticpcadata %>% unite("plotid", Site:Plot, remove = "false")
#Removing Site and Plot
abioticpcatrim <- abioticpcadata %>% select(-c(Site, Plot))

#This assigns the rowsnames so that the datapoints will come up as the site/plot names on the plot
rownames(abioticpcatrim) <- abioticpcatrim$plotid
#This drops the first column (old rownames/plotid)
abioticpcatrim <- as.data.frame(abioticpcatrim[,-1])
### Renaming columns for nicer plot
colnames(abioticpcatrim) <- c("pH", "NH4N", "NO3N", "P", "K", "canopy closure", "litter cover")
soil_pca <- princomp(abioticpcatrim, cor = TRUE)

dev.off()
pdf("Output/Figures/pca-abiotic.pdf")
par(pty="s")
autoplot(soil_pca, label = TRUE, shape = TRUE,
         loadings = TRUE, loadings.colour = 'slateblue', 
         loadings.label.repel = TRUE, loadings.label.size = 5,
         loadings.label = TRUE, loadings.label.colour = 'slateblue')+
  xlab("PC1 (54.1%)")+
  ylab("PC2 (18.5%)")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16))
dev.off()
#Summary and loadings info
summary(soil_pca)
loadings(soil_pca)

#### SUPP. TABLE S1 - p-values quadratic fixed effects ####
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

#No neighbour abundance analysis for germination

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
#Note that HYGL doesn't converge here, not quadratic (plotted to check)
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
#Struggling to save this using save_kable and as_image()
#Can copy it from the Viewer using copy to clipboard, maintain aspect ratio, first value 500
#### SUPP. TABLE S2 - below model outputs table ####
#### Extracting marginal and conditional r squared values for supp. table s2
### Need to run models first
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

## Germination ####
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

## Survival #### 
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

## Fecundity #### 
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

## Lambda ####  
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





#### SUPP. FIGURE S4 - interactions PC1 and watering ####
### Survival, seed production and pop growth ~ PC1 and watering
##Data frame for wet and dry colours
colours <- data.frame(Treatment=c("Dry","Ambient","Wet"),
                      colour = c("#CC79A7","grey60","#0072B2"))
dev.off()
pdf("Output/Figures/supp_PC1_watering2.pdf", width=20, height=18)
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

#### SUPP. FIGURE S5 - LTRE bar chart ####
ltre_sp <- c("ARCA", "PLDE", "VERO")
ltre_sp <- as.data.frame(ltre_sp)
#Add in contribution values
ltre_sp$a12 <- c(0.31, 0.73, 0.26)
ltre_sp$a21 <- c(0.52, 0.16, 0.16)
ltre_sp$a22 <- c(0.17, 0.11, 0.58)

#Need long form
ltre_sp <- ltre_sp %>% pivot_longer(cols = c(a12, a21, a22), names_to = "transition", values_to = "contribution")
str(ltre_sp)
## Make barchart
## Actually, vero has an interaction so no signif main effect of neighbours on lambda!
ltre_sp_no_vero <- ltre_sp %>% filter(!(ltre_sp=="VERO"))
dev.off()
pdf("Output/Figures/supp_8_ltre.pdf", width=8, height=6)
ggplot(ltre_sp_no_vero, aes(x=transition, y=contribution))+
  geom_bar(aes(fill=ltre_sp), position = "dodge2", stat = "identity")+
  ylab("Contribution to difference in population growth rates")+
  xlab("Transition")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), labels=c(expression(italic('A. calendula')), expression(italic('P. debilis'))))+
  theme_classic()+
  guides(fill = guide_legend(title = "Species", label.hjust = 0, title.hjust=0.5)) +
  theme(axis.text = element_text(size=14),
        legend.text = element_text(size=14), legend.title = element_text(size=14),
        axis.title=element_text(size=14))
dev.off()

##### Decomposing differences in lambda in presence and absence of neighbours ####
## for species that had a signif affect of nbhs on lambda, so for ARCA and PLDE
#### Matrix elements in order of a11, a21, a12, a22 where:

#a: dormancy rate of seed bank seeds
#z: seed survival rate (maintaining viability)
#f: number of viable seeds produced per adult
#g: rate of seedling emergence
#s: plant survivorship to reproduction

#a11 = az
#a12 = afz
#a21 = (1-a)zgs
#a22 = (1-a)zgsf

#Rates common to both presence and absence of neighbours: 
# dormancy rate (a), emergence rate (g), seed survival (z)
## z - seed survival rate from Isaac Towers' work
#species_level_seed_fill data frame

## Need to calculate average germination rates and 
#fecundity and survival rates in the presence or absence of neighbours

sp_ltre <- popdata %>% group_by(Species, Neighbours01) %>%
  summarise(avg_emerg = mean(plot_germ, na.rm=T),
            avg_surv = mean(plot_survival, na.rm=T),
            avg_fecundity = mean(plot_fecundity, na.rm=T),
            lambda = mean(lambda, na.rm=T))

### ARCA ###
# ARCA - Neighbours0 - absence of neighbours
#a: 0.76 (average from other species in Maia Raymundo's 2016 Perenjori dataset)
#z: 0.91
#g: 0.24
#f: 10.17
#s: 0.48
# ARCA - Neighbours1 - presence of neighbours
#a: 0.76 (average from other species in Maia Raymundo's 2016 Perenjori dataset)
#z: 0.91
#g: 0.24
#f: 6.23
#s: 0.20

#Build matrices using rates for absence or presence of neighbours
#A11, A21, A12, A22
#data=c(az, (1-a)zgs, afz, (1-a)zgsf)
arca_no_nbh<- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.48, 0.76*10.17*0.91, (1-0.76)*0.91*0.24*0.48*10.17), nrow=2, ncol=2)
arca_nbh <- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.20, 0.76*6.23*0.91, (1-0.76)*0.91*0.24*0.20*6.23), nrow=2, ncol=2)
#Calculate eigenvalues for each matrix (lambda, largest eigenvalue)
eigen(arca_no_nbh) #0.95
eigen(arca_nbh) #0.76

##Using exactLTRE to perform fixed directional LTRE analysis with no nbh as reference matrix
#Calculating contributions to the difference in lambda due to vital rates
#calcuated by setting all other parameters to their mean

cont_diff <- exactLTRE(list(arca_no_nbh, arca_nbh), method = 'fixed', fixed.directional = TRUE)
# matrix of contributions
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# sum of contributions
sum(cont_diff$epsilons) 
# true difference in lambda  #-0.191
lamDiff(list(arca_no_nbh,arca_nbh)) 

#approximate LTRE approach and popbio LTRE give the same values

### PLDE ###
# PLDE - Neighbours0 - absence of neighbours
#a: 0.76 (average from other species in Maia's 2016 Perenjori dataset)
#z: 0.98
#g: 0.33
#f: 8.02
#s: 0.25
# PLDE - Neighbours1 - presence of neighbours
#a: 0.76 (average from other species in Maia's 2016 Perenjori dataset)
#z: 0.98
#g: 0.33
#f: 2.55
#s: 0.20

#A11, A21, A12, A22
#data=c(az, (1-a)zgs, afz, (1-a)zgsf)
plde_no_nbh<- matrix(data=c(0.76*0.98, (1-0.76)*0.98*0.33*0.25, 0.76*8.02*0.98, (1-0.76)*0.98*0.33*0.25*8.02), nrow=2, ncol=2)
plde_nbh <- matrix(data=c(0.76*0.98, (1-0.76)*0.98*0.33*0.20, 0.76*2.55*0.98, (1-0.76)*0.98*0.33*0.20*2.55), nrow=2, ncol=2)
eigen(plde_no_nbh) #0.90
eigen(plde_nbh) #0.78
# true difference in lambda  #-0.116
lamDiff(list(plde_no_nbh,plde_nbh)) 

cont_diff <- exactLTRE(list(plde_no_nbh, plde_nbh), method = 'fixed')
# matrix of contributions
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# sum of contributions
sum(cont_diff$epsilons)

LTRE(plde_nbh, plde_no_nbh)
approximateLTRE(list(plde_no_nbh,plde_nbh), method='fixed') # contributions to the difference in lambda

