## Full model for WA data
#July - October 2021

library(tidyverse)  
library(lme4)
library(ggplot2)
library(lmerTest)

#Data imported from data preparation sheet
source("data_preparation.R")
#dataall has everything (germination, survival, seed production and neighbour info) combined
#datanotscaled is same dataset but without standardised predictors

################################### MODELS ###################################
####Survival models (ProducedSeeds as proxy) #####
#Survival to seed production, yes or no
#Note that Control watering treatment is the reference
#Many died before the watering treatment began? Need to account for this!!!!**
##ARCA
#Only one case of one intraspecific neighbour
arcasurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, arcamortality)
summary(arcasurv1)
arcasurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, arcamortality)
summary(arcasurv2)
#Warning messages
arcasurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, arcamortality)
summary(arcasurv3)
arcasurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 +  
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, arcamortality)
#Doesn't converge
arcasurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, arcamortality)
summary(arcasurv5)
arcasurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, arcamortality)
summary(arcasurv6)
#Warning messages
AIC(arcasurv1, arcasurv2, arcasurv3, arcasurv5, arcasurv6)
#arcasurv5 best fit

###hygl
hyglsurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, hyglmortality)
summary(hyglsurv1)
hyglsurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, hyglmortality)
summary(hyglsurv2)
hyglsurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, hyglmortality)
summary(hyglsurv3)
hyglsurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, hyglmortality)
summary(hyglsurv4)
hyglsurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, hyglmortality)
#Doesn't converge
hyglsurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, hyglmortality)
summary(hyglsurv6)
AIC(hyglsurv1, hyglsurv2, hyglsurv3, hyglsurv4, hyglsurv6)
###hyglsurv3 most parsimonious

#laro
#No intraspecific neighbours
larosurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, laromortality)
summary(larosurv1)
larosurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, laromortality)
summary(larosurv3)
larosurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, laromortality)
summary(larosurv5)
larosurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, laromortality)
summary(larosurv6)
AIC(larosurv1, larosurv3, larosurv5, larosurv6)
#larosurv1 most parsimonious

#peai
peaisurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, peaimortality)
summary(peaisurv1)
peaisurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, peaimortality)
summary(peaisurv2)
peaisurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, peaimortality)
summary(peaisurv3)
peaisurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, peaimortality)
summary(peaisurv4)
peaisurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, peaimortality)
summary(peaisurv5)
peaisurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, peaimortality)
summary(peaisurv6)
AIC(peaisurv1, peaisurv2, peaisurv3, peaisurv4, peaisurv5, peaisurv6)
#peaisurv1 most parsimonious

#plde
pldesurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, pldemortality)
summary(pldesurv1)
pldesurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, pldemortality)
summary(pldesurv2)
pldesurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, pldemortality)
summary(pldesurv3)
pldesurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, pldemortality)
summary(pldesurv4)
pldesurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, pldemortality)
summary(pldesurv5)
pldesurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, pldemortality)
summary(pldesurv6)
AIC(pldesurv1, pldesurv2, pldesurv3, pldesurv4, pldesurv5, pldesurv6)
#peaisurv1 most parsimonious

#pole
#Only one case of intraspecific neighbours (15)
polesurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, polemortality)
summary(polesurv1)
polesurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, polemortality)
summary(polesurv2)
polesurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, polemortality)
summary(polesurv3)
polesurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, polemortality)
summary(polesurv4)
#Warning messages
polesurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, polemortality)
#Doesn't converge
polesurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, polemortality)
summary(polesurv6)
#Rank deficient
AIC(polesurv1, polesurv2, polesurv3, polesurv4)
#polesurv1 most parsimonious

#trcy
trcysurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, trcymortality)
summary(trcysurv1)
trcysurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, trcymortality)
summary(trcysurv2)
trcysurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, trcymortality)
summary(trcysurv3)
trcysurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, trcymortality)
summary(trcysurv4)
trcysurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, trcymortality)
summary(trcysurv5)
trcysurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, trcymortality)
#Doesn't converge
AIC(trcysurv1, trcysurv2, trcysurv3, trcysurv4, trcysurv5)
#trcysurv1 most parsimonious

#tror
#Only one case of intraspecific neighbours (1)
trorsurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, trormortality)
summary(trorsurv1)
trorsurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, trormortality)
#Doesn't converge
trorsurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, trormortality)
summary(trorsurv3)
trorsurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, trormortality)
summary(trorsurv4)
#Warning messages
trorsurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, trormortality)
summary(trorsurv5)
trorsurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, trormortality)
summary(trorsurv6)
AIC(trorsurv1, trorsurv3, trorsurv4, trorsurv5, trorsurv6)
#trorsurv3 most parsimonious

#vero
verosurv1 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Total_abundance + (1|Site/Plot), family = binomial, veromortality)
summary(verosurv1)
verosurv2 <- glmer(ProducedSeeds ~ cc_percentage + Treatment + Intra_abundance + 
                     Inter_abundance + (1|Site/Plot), family = binomial, veromortality)
summary(verosurv2)
verosurv3 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                     (1|Site/Plot), family = binomial, veromortality)
summary(verosurv3)
verosurv4 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment)^2 + 
                     Intra_abundance + Inter_abundance + (1|Site/Plot), family = binomial, veromortality)
summary(verosurv4)
verosurv5 <- glmer(ProducedSeeds ~ (cc_percentage + Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, veromortality)
summary(verosurv5)
verosurv6 <- glmer(ProducedSeeds ~ cc_percentage + (Treatment + Total_abundance)^2 + 
                     (1|Site/Plot), family = binomial, veromortality)
summary(verosurv6)
AIC(verosurv1, verosurv2, verosurv3, verosurv4, verosurv5, verosurv6)
#verosurv1 most parsimonious

############ Seed production data models ######################
# If a plant survived to produce any seeds (whether viable or inviable), ProducedSeeds = 1

seedmodeldata <- dataall %>% filter(ProducedSeeds == "1")
seedarca <- seedmodeldata %>% filter(Species == "ARCA")
seedhygl <- seedmodeldata %>% filter(Species == "HYGL")
seedlaro <- seedmodeldata %>% filter(Species == "LARO")
seedpeai <- seedmodeldata %>% filter(Species == "PEAI")
seedplde <- seedmodeldata %>% filter(Species == "PLDE")
seedpole <- seedmodeldata %>% filter(Species == "POLE")
seedtrcy <- seedmodeldata %>% filter(Species == "TRCY")
seedtror <- seedmodeldata %>% filter(Species == "TROR")
seedvero <- seedmodeldata %>% filter(Species == "VERO")

ggplot(seedmodeldata, aes(x = log(No_viable_seeds_grouped+1), y = log(No_inviable_seeds_grouped+1)))+
  geom_point(aes(alpha = 0.4))+
  geom_smooth(method = "lm")+
  theme_classic()+
  facet_wrap(vars(Species))
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  corr <- lm(log(No_inviable_seeds_grouped+1) ~ log(No_viable_seeds_grouped+1), data = filter(seedmodeldata, Species == specieslist[i]))
  print(summary(corr))
}

## Model 6 is asking whether water availability interacts with neighbour abundance to affect fecundity
#glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot)

###arca
#ARCA doesn't have any intraspecific neighbours so in models total_abundance = intra_abundance
arca1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedarca)
summary(arca1)
arca3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedarca)
summary(arca3)
#Doesn't converge
arca5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedarca)
#Doesn't converge
arca6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedarca)
summary(arca6)
AIC(arca1, arca6)
#arca1 most parsimonious model

###hygl
hygl1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedhygl)
summary(hygl1)
hygl2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedhygl)
#Doesn't converge
hygl3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedhygl)
#Doesn't converge
hygl4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedhygl)
#Doesn't converge
hygl5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedhygl)
#Doesn't converge
hygl6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedhygl)
#Doesn't converge
#hygl1 best fit, only one that converged

###laro
#Has no intraspecific neighbours!
laro1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedlaro)
summary(laro1)
laro3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedlaro)
summary(laro3)
laro5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedlaro)
summary(laro5)
#Doesn't converge
laro6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedlaro)
summary(laro6)
AIC(laro1, laro2, laro3, laro4, laro6)
#laro1 best fit

###peai
peai1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedpeai)
summary(peai1)
peai2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedpeai)
summary(peai2)
peai3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedpeai)
summary(peai3)
peai4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedpeai)
summary(peai4)
peai5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedpeai)
summary(peai5)
peai6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedpeai)
summary(peai6)
AIC(peai1, peai2, peai3, peai4, peai5, peai6)
#peai2 best fit

###plde
plde1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedplde)
summary(plde1)
plde2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedplde)
summary(plde2)
plde3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedplde)
summary(plde3)
plde4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedplde)
summary(plde4)
plde5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedplde)
#Doesn't converge
plde6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedplde)
summary(plde6)
AIC(plde1, plde2, plde3, plde4, plde6)
#plde1 most parsimonious

###pole
pole1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedpole)
summary(pole1)
pole2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedpole)
summary(pole2)
pole3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedpole)
summary(pole3)
pole4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedpole)
summary(pole4)
pole5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedpole)
#Doesn't converge
pole6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedpole)
summary(pole6)
AIC(pole1, pole2, pole3, pole4, pole5, pole6)
#pole1 most parsimonious

###trcy
trcy1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedtrcy)
summary(trcy1)
trcy2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedtrcy)
summary(trcy2)
trcy3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedtrcy)
summary(trcy3)
trcy4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedtrcy)
summary(trcy4)
trcy5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedtrcy)
summary(trcy5)
trcy6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedtrcy)
summary(trcy6)
AIC(trcy1, trcy2, trcy3, trcy4, trcy5, trcy6)
#trcy1 most parsimonious

###tror
#Only one occurrence of intra_abundance (1)
tror1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedtror)
summary(tror1)
tror2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedtror)
#Doesn't converge
tror3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedtror)
summary(tror3)
tror4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedtror)
#Doesn't converge
tror5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedtror)
#Doesn't converge
tror6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedtror)
summary(tror6)
#TreatmentDry, TreatmentDry:Total_abundance
AIC(tror1, tror3, tror6)
#tror3 best fit

###vero
vero1 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                    (1|Site/Plot), seedvero)
#Doesn't converge
vero2 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedvero)
#Doesn't converge
vero3 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                    (1|Site/Plot), seedvero)
summary(vero3)
vero4 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                    Inter_abundance + (1|Site/Plot), seedvero)
summary(vero4)
vero5 <- glmer.nb(No_viable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedvero)
summary(vero5)
vero6 <- glmer.nb(No_viable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedvero)
summary(vero6)
AIC(vero3, vero4, vero5, vero6)
#vero3 most parsimonious

#### Inviable seed models ####

###arca
#ARCA doesn't have any intraspecific neighbours
arcain1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedarca)
summary(arcain1)
arcain3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedarca)
summary(arcain3)
arcain5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedarca)
summary(arcain5)
arcain6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedarca)
summary(arcain6)
AIC(arcain1, arcain3, arcain5, arcain6)
#arcain3 best model

###hygl
hyglin1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedhygl)
summary(hyglin1)
hyglin2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedhygl)
summary(hyglin2)
hyglin3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedhygl)
summary(hyglin3)
hyglin4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedhygl)
summary(hyglin4)
hyglin5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedhygl)
summary(hyglin5)
hyglin6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedhygl)
summary(hyglin6)
AIC(hyglin1, hyglin2, hyglin3, hyglin4, hyglin5, hyglin6)
#hyglin1 most parsimonious

###laro
#No intraspecific neighbours
laroin1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedlaro)
summary(laroin1)
laroin3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedlaro)
summary(laroin3)
laroin5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedlaro)
summary(laroin5)
laroin6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedlaro)
summary(laroin6)
AIC(laroin1, laroin3, laroin5, laroin6)
#laroin1 most parsimonious

###peai
peaiin1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedpeai)
summary(peaiin1)
peaiin2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedpeai)
summary(peaiin2)
peaiin3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedpeai)
#Doesn't converge
peaiin4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedpeai)
summary(peaiin4)
peaiin5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedpeai)
summary(peaiin5)
peaiin6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedpeai)
summary(peaiin6)
AIC(peaiin1, peaiin2, peaiin4, peaiin5, peaiin6)
#peaiin1 most parsimonious

###plde
pldein1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedplde)
#Doesn't converge
pldein2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedplde)
#Doesn't converge
pldein3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedplde)
summary(pldein3)
pldein4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedplde)
#Doesn't converge
pldein5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedplde)
summary(pldein5)
pldein6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedplde)
#Doesn't converge
AIC(pldein3, pldein5)
#pldein5 best fit? Warning messages. Check this. I think pldein3 only model that converges

###pole
polein1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedpole)
summary(polein1)
polein2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedpole)
summary(polein2)
polein3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedpole)
summary(polein3)
polein4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedpole)
summary(polein4)
polein5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedpole)
#Doesn't converge
polein6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedpole)
summary(polein6)
AIC(polein1, polein2, polein3, polein4, polein6)
#polein1 most parsimonious

###trcy
trcyin1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedtrcy)
summary(trcyin1)
trcyin2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedtrcy)
summary(trcyin2)
trcyin3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedtrcy)
#Doesn't converge
trcyin4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedtrcy)
summary(trcyin4)
trcyin5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedtrcy)
#Doesn't converge
trcyin6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedtrcy)
summary(trcyin6)
AIC(trcyin1, trcyin2, trcyin4, trcyin6)
#trcyin1 most parsimonious

###tror
#Only one occurence of intra_abundance
trorin1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|Site/Plot), seedtror)
summary(trorin1)
trorin2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedtror)
summary(trorin2)
trorin3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedtror)
summary(trorin3)
trorin4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedtror)
summary(trorin4)
trorin5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedtror)
#Doesn't converge
trorin6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedtror)
summary(trorin6)
AIC(trorin1, trorin2, trorin3, trorin4, trorin6)
#trorin6 best fit?

###vero
veroin1 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Total_abundance + 
                      (1|plotid), seedvero)
summary(veroin1)
veroin2 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + Treatment + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedvero)
summary(veroin2)
veroin3 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Total_abundance + 
                      (1|Site/Plot), seedvero)
summary(veroin3)
veroin4 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment)^2 + Intra_abundance +
                      Inter_abundance + (1|Site/Plot), seedvero)
summary(veroin4)
veroin5 <- glmer.nb(No_inviable_seeds_grouped ~ (cc_percentage + Treatment + Total_abundance)^2 + (1|Site/Plot), seedvero)
summary(veroin5)
veroin6 <- glmer.nb(No_inviable_seeds_grouped ~ cc_percentage + (Treatment + Total_abundance)^2 + (1|Site/Plot), seedvero)
summary(veroin6)
AIC(veroin1, veroin2, veroin3, veroin4, veroin5, veroin6)
#veroin1 most parsimonious

#### Plotting coefficients using sjPlot ####

library(sjPlot)
library(gridExtra)
# Viable seed coef plots ####
#These are the best fitting models
allmodels <- list()
allmodels[[1]] <- arca1
allmodels[[2]] <- hygl1
allmodels[[3]] <- laro1
allmodels[[4]] <- peai2
allmodels[[5]] <- plde1
allmodels[[6]] <- pole1
allmodels[[7]] <- trcy1
allmodels[[8]] <- tror3
allmodels[[9]] <- vero3

## Splitting into two panels to spread out info. 
# Panel A - Canopy, Dry, Wet, Canopy:Dry, Canopy:Wet
#Panel B - Total_abundance, Intra_abundance, Inter_abundance
coefplot1 <- plot_models(allmodels, transform = NULL, vline.color = "grey", show.legend = FALSE, 
            rm.terms = c("Total_abundance", "Intra_abundance", "Inter_abundance"),
            axis.labels=c("Canopy:Wet", "Canopy:Dry", "Wet", "Dry", "Canopy Cover"),
            axis.lim=c(-5,5), dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
coefplot2 <- plot_models(allmodels, transform = NULL, vline.color = "grey", legend.title = "Species", 
            rm.terms = c("cc_percentage", "TreatmentDry", "TreatmentWet", "cc_percentage:TreatmentDry", "cc_percentage:TreatmentWet"), 
            axis.labels=c("Interspecific abundance", "Intraspecific abundance", "Total abundance"),
            dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))

grid.arrange(coefplot1, coefplot2, ncol = 2)
#Want to reduce spacing in between terms and add in grey horizontal lines separating them

#Trying a multi panel plot
a <- plot_model(arca1, transform = NULL, title = "ARCA", vline.color = "grey") + ylab("Estimate") +theme_classic()
b <- plot_model(plde1, transform = NULL, title = "PLDE", vline.color = "grey") + ylab("Estimate")+theme_classic()
c <- plot_model(pole1, transform = NULL, title = "POLE", vline.color = "grey")+ ylab("Estimate")+theme_classic()
d <- plot_model(laro1, transform = NULL, title = "POLE", vline.color = "grey")+ ylab("Estimate")+theme_classic()
e <- plot_model(trcy1, transform = NULL, title = "TRCY", vline.color = "grey")+ ylab("Estimate")+theme_classic()
f <- plot_model(vero3, transform = NULL, title = "VERO", vline.color = "grey")+ ylab("Estimate")+theme_classic()
g <- plot_model(hygl3, transform = NULL, title = "HYGL", vline.color = "grey")+ ylab("Estimate")+theme_classic()
h <- plot_model(tror3, transform = NULL, title = "TROR", vline.color = "grey")+ ylab("Estimate")+theme_classic()
i <- plot_model(peai2, transform = NULL, title = "PEAI", vline.color = "grey")+ ylab("Estimate")+theme_classic()

plot_grid(list(a, b, c, d, e, f, g, h, i), margin = c(.1, .1, .1, .1))

#Inviable seed coef plots ####

#These are the best fitting models
#This is what it was like with trorin1 being the best fit, trickier with trorin6
inmodels <- list()
inmodels[[1]] <- arcain3
inmodels[[2]] <- hyglin1
inmodels[[3]] <- laroin1
inmodels[[4]] <- peaiin1
inmodels[[5]] <- pldein3
inmodels[[6]] <- polein1
inmodels[[7]] <- trcyin1
inmodels[[8]] <- trorin1
inmodels[[9]] <- veroin1

coefplotin1 <- plot_models(inmodels, transform = NULL, vline.color = "grey", show.legend = FALSE, 
                         rm.terms = c("Total_abundance", "cc_percentage:TreatmentWet", "cc_percentage:TreatmentDry"),
                         axis.labels=c("Wet", "Dry", "Canopy Cover"),
                         axis.lim=c(-4,4), dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
coefplotin2 <- plot_models(inmodels, transform = NULL, vline.color = "grey", legend.title = "Species", 
                         rm.terms = c("cc_percentage", "TreatmentDry", "TreatmentWet"), 
                         axis.labels=c("Canopy:Wet", "Canopy:Dry", "Total abundance"),
                         dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
grid.arrange(coefplotin1, coefplotin2, ncol = 2)

## Not quite right yet - trorin6 best model:
# inmodels <- list()
# inmodels[[1]] <- arcain3
# inmodels[[2]] <- hyglin1
# inmodels[[3]] <- laroin1
# inmodels[[4]] <- peaiin1
# inmodels[[5]] <- pldein3
# inmodels[[6]] <- polein1
# inmodels[[7]] <- trcyin1
# inmodels[[8]] <- trorin6
# inmodels[[9]] <- veroin1
# 
# coefplotin1 <- plot_models(inmodels, transform = NULL, vline.color = "grey", show.legend = FALSE, 
#                            rm.terms = c("cc_percentage:TreatmentWet", "cc_percentage:TreatmentDry", "TreatmentDry:Total_abundance", "TreatmentWet:Total_abundance"),
#                            axis.labels=c("Total abundance", "Wet", "Dry", "Canopy Cover"),
#                            axis.lim=c(-4,4), dot.size = 2, line.size = 1)+
#   ylab("Estimate")+
#   theme_classic()+
#   scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
# coefplotin2 <- plot_models(inmodels, transform = NULL, vline.color = "grey", legend.title = "Species", 
#                            rm.terms = c("cc_percentage", "TreatmentDry", "TreatmentWet", "Total_abundance"), 
#                            axis.labels=c("Canopy:Wet", "Canopy:Dry", "Dry:Abundance", "Wet:Abundance"),
#                            dot.size = 2, line.size = 1)+
#   ylab("Estimate")+
#   theme_classic()+
#   scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
# grid.arrange(coefplotin1, coefplotin2, ncol = 2)

# Survival coef plots ####
#plot_models(arcasurv5, hyglsurv3, larosurv1, peaisurv1, pldesurv1, polesurv1, trcysurv1, trorsurv3, verosurv1, transform = NULL, vline.color = "grey", legend.title = "Species")

survivalmodels <- list()
survivalmodels[[1]] <- arcasurv5
survivalmodels[[2]] <- hyglsurv3
survivalmodels[[3]] <- larosurv1
survivalmodels[[4]] <- peaisurv1
survivalmodels[[5]] <- pldesurv1
survivalmodels[[6]] <- polesurv1
survivalmodels[[7]] <- trcysurv1
survivalmodels[[8]] <- trorsurv3
survivalmodels[[9]] <- verosurv1

coefsurvplot1 <- plot_models(survivalmodels, transform = NULL, vline.color = "grey", legend.title = "Species", 
                         rm.terms = c("cc_percentage:Total_abundance", "TreatmentDry:Total_abundance", "TreatmentWet:Total_abundance"),
                         axis.labels=c("Canopy:Wet", "Canopy:Dry", "Total_abundance", "Wet", "Dry", "Canopy Cover"),
                         axis.lim=c(-5,5), dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
coefsurvplot2 <- plot_models(survivalmodels, transform = NULL, vline.color = "grey", show.legend = FALSE, 
                         rm.terms = c("cc_percentage", "TreatmentDry", "TreatmentWet", "Total_abundance", "cc_percentage:TreatmentDry", "cc_percentage:TreatmentWet"), 
                         axis.labels=c("Canopy cover", "Dry:abundance", "Wet:abundance"),
                         dot.size = 2, line.size = 1)+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
grid.arrange(coefsurvplot1, coefsurvplot2, ncol = 2)

#Trying this another way
#Still need to try panels by sp.
survalist <- list()
survalist[[1]] <- larosurv1
survalist[[2]] <- peaisurv1
survalist[[3]] <- pldesurv1
survalist[[4]] <- polesurv1
survalist[[5]] <- trcysurv1
survalist[[6]] <- verosurv1
survblist <- list()
survblist[[1]] <- hyglsurv3
survblist[[2]] <- trorsurv3
survblist[[3]] <- arcasurv5

plot_models(survalist, transform = NULL, vline.color = "grey", legend.title = "Species",
                     axis.labels=c("Total_abundance", "Wet", "Dry", "Canopy Cover"),
                     axis.lim=c(-4,4), dot.size = 2, line.size = 1)+
                      ylab("Estimate")+
                      theme_classic()+
                  scale_colour_discrete(labels = c("VERO", "TRCY", "POLE", "PLDE", "PEAI", "LARO"))
plot_models(survblist, transform = NULL, vline.color = "grey", legend.title = "Species",
                     axis.labels=c("Wet:Abundance", "Dry:Abundance", "Canopy:Abundance", "Canopy:Wet", "Canopy:Dry", "Total_abundance", "Wet", "Dry", "Canopy Cover"),
                     dot.size = 2, line.size = 1)+
                      ylab("Estimate")+
                      theme_classic()+
                      scale_colour_discrete(labels = c("ARCA", "TROR", "HYGL"))
# Why is the ARCA dry:abundance estimate so high?!

## Viable seed curve plots ####
##John's code 20th July 2021
#Need to plot relationship between fecundity and variables to see how coefficient adjusts relationship
#1,x, 0, 0: 1 = intercept on. x = x. 0, 0 = dry and wet off. So this is for the control data
#1, x, 1, 0: Dry on, wet off. 
#1, x, 0, 1: Dry off, wet on.
#Black is control plots, red is dry plots, green is wet plots

par(mfrow=c(3,3))
summary(arca1)
with(seedarca, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='# viable seeds', xlab='', main = 'ARCA'))
curve(exp(cbind(1,x,0,0,mean(seedarca$Total_abundance))%*%fixef(arca1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedarca$Total_abundance))%*%fixef(arca1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedarca$Total_abundance))%*%fixef(arca1)), add=T, col = 3)

summary(hygl1)
with(seedhygl, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab = '', xlab='', main = 'HYGL'))
curve(exp(cbind(1,x,0,0,mean(seedhygl$Total_abundance))%*%fixef(hygl1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedhygl$Total_abundance))%*%fixef(hygl1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedhygl$Total_abundance))%*%fixef(hygl1)), add=T, col = 3)

summary(laro1)
with(seedlaro, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='', main = 'LARO'))
curve(exp(cbind(1,x,0,0,mean(seedlaro$Total_abundance))%*%fixef(laro1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedlaro$Total_abundance))%*%fixef(laro1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedlaro$Total_abundance))%*%fixef(laro1)), add=T, col = 3)

summary(peai2)
with(seedpeai, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='# viable seeds', xlab='', main = 'PEAI'))
curve(exp(cbind(1,x,0,0,mean(seedpeai$Intra_abundance), mean(seedpeai$Inter_abundance))%*%fixef(peai2)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedpeai$Intra_abundance), mean(seedpeai$Inter_abundance))%*%fixef(peai2)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedpeai$Intra_abundance), mean(seedpeai$Inter_abundance))%*%fixef(peai2)), add=T, col = 3)

summary(plde1)
with(seedplde, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='', main = 'PLDE'))
curve(exp(cbind(1,x,0,0,mean(seedplde$Total_abundance))%*%fixef(plde1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedplde$Total_abundance))%*%fixef(plde1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedplde$Total_abundance))%*%fixef(plde1)), add=T, col = 3)

summary(pole1)
with(seedpole, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='', main = 'POLE'))
curve(exp(cbind(1,x,0,0,mean(seedpole$Total_abundance))%*%fixef(pole1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedpole$Total_abundance))%*%fixef(pole1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedpole$Total_abundance))%*%fixef(pole1)), add=T, col = 3)

summary(trcy1)
with(seedtrcy, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='# viable seeds', xlab='Canopy cover', main = 'TRCY'))
curve(exp(cbind(1,x,0,0,mean(seedtrcy$Total_abundance))%*%fixef(trcy1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedtrcy$Total_abundance))%*%fixef(trcy1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedtrcy$Total_abundance))%*%fixef(trcy1)), add=T, col = 3)

summary(tror3)
with(seedtror, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment),ylab='', xlab='Canopy cover', main = 'TROR'))
curve(exp(cbind(1,x,0,0,mean(seedtror$Total_abundance), x*0, x*0)%*%fixef(tror3)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedtror$Total_abundance), x*1, x*0)%*%fixef(tror3)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedtror$Total_abundance), x*0, x*1)%*%fixef(tror3)), add=T, col = 3)

summary(vero3)
with(seedvero, plot(No_viable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='Canopy cover', main = 'VERO'))
curve(exp(cbind(1,x,0,0,mean(seedvero$Total_abundance), x*0, x*0)%*%fixef(vero3)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedvero$Total_abundance), x*1, x*0)%*%fixef(vero3)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedvero$Total_abundance), x*0, x*1)%*%fixef(vero3)), add=T, col = 3)

### Inviable seed curve plots ####
par(mfrow=c(1,1))
par(mfrow=c(3,3))
summary(arcain3)
with(seedarca, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='# inviable seeds', xlab='', main = 'ARCA'))
curve(exp(cbind(1,x,0,0,mean(seedarca$Total_abundance), x*0, x*0)%*%fixef(arcain3)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedarca$Total_abundance), x*1, x*0)%*%fixef(arcain3)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedarca$Total_abundance), x*0, x*1)%*%fixef(arcain3)), add=T, col = 3)

summary(hyglin1)
with(seedhygl, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab = '', xlab='', main = 'HYGL'))
curve(exp(cbind(1,x,0,0,mean(seedhygl$Total_abundance))%*%fixef(hyglin1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedhygl$Total_abundance))%*%fixef(hyglin1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedhygl$Total_abundance))%*%fixef(hyglin1)), add=T, col = 3)

summary(laroin1)
with(seedlaro, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='', main = 'LARO'))
curve(exp(cbind(1,x,0,0,mean(seedlaro$Total_abundance))%*%fixef(laroin1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedlaro$Total_abundance))%*%fixef(laroin1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedlaro$Total_abundance))%*%fixef(laroin1)), add=T, col = 3)

summary(peaiin1)
with(seedpeai, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='# inviable seeds', xlab='', main = 'PEAI'))
curve(exp(cbind(1,x,0,0,mean(seedpeai$Total_abundance))%*%fixef(peaiin1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedpeai$Total_abundance))%*%fixef(peaiin1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedpeai$Total_abundance))%*%fixef(peaiin1)), add=T, col = 3)

summary(pldein3)
with(seedplde, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='', main = 'PLDE'))
curve(exp(cbind(1,x,0,0,mean(seedplde$Total_abundance), x*0, x*0)%*%fixef(pldein3)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedplde$Total_abundance), x*1, x*0)%*%fixef(pldein3)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedplde$Total_abundance), x*0, x*1)%*%fixef(pldein3)), add=T, col = 3)

summary(polein1)
with(seedpole, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='', main = 'POLE'))
curve(exp(cbind(1,x,0,0,mean(seedpole$Total_abundance))%*%fixef(polein1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedpole$Total_abundance))%*%fixef(polein1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedpole$Total_abundance))%*%fixef(polein1)), add=T, col = 3)

summary(trcyin1)
with(seedtrcy, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='# inviable seeds', xlab='Canopy cover', main = 'TRCY'))
curve(exp(cbind(1,x,0,0,mean(seedtrcy$Intra_abundance))%*%fixef(trcyin1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedtrcy$Intra_abundance))%*%fixef(trcyin1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedtrcy$Intra_abundance))%*%fixef(trcyin1)), add=T, col = 3)

summary(trorin1)
with(seedtror, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment),ylab='', xlab='Canopy cover', main = 'TROR'))
curve(exp(cbind(1,x,0,0,mean(seedtror$Total_abundance))%*%fixef(trorin1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedtror$Total_abundance))%*%fixef(trorin1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedtror$Total_abundance))%*%fixef(trorin1)), add=T, col = 3)

summary(veroin1)
with(seedvero, plot(No_inviable_seeds_grouped ~ cc_percentage, col=as.factor(Treatment), ylab='', xlab='Canopy cover', main = 'VERO'))
curve(exp(cbind(1,x,0,0,mean(seedvero$Total_abundance))%*%fixef(veroin1)), add=T, col = 1)
curve(exp(cbind(1,x,1,0,mean(seedvero$Total_abundance))%*%fixef(veroin1)), add=T, col = 2)
curve(exp(cbind(1,x,0,1,mean(seedvero$Total_abundance))%*%fixef(veroin1)), add=T, col = 3)


## Coefficient summary tables ####
#Viable seeds
tab_model(arca1, hygl1, laro1, peai2, plde1, pole1, trcy1, tror3, vero3, transform = NULL)
#Inviable seeds
tab_model(arcain3, hyglin1, laroin1, peaiin1, pldein3, polein1, trcyin2, trorin6, veroin1, transform = NULL)
#Survival
tab_model(arcasurv5, hyglsurv3, larosurv1, peaisurv1, pldesurv1, polesurv1, trcysurv1, trorsurv3, verosurv1, transform = NULL)

### Diversity as a predictor #####
# Work in progress!

library(MuMIn)
r.squaredGLMM(arca1)

#A way to test how much variance is being soaked up by random effect
lmmod <-  lm(No_viable_seeds_grouped ~ 1 + as.factor(plotid), data=seedarca)
summary(lmmod)






