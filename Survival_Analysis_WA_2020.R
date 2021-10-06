#### 02/11/2020 Ali Catling Perenjori Experiment 2020
## Analysis of survival data for ESA speed talk

library(tidyverse)
library(lme4)
library(ggplot2)
library(interactions)
library(emmeans)

#### FORMING DATASETS

#mortality data
#Note that this is not the most up to date -- want survival to produce viable seeds
mortalitydataraw <- read_csv("Data/mortalityofgerminated.csv")
## Watering and shade/sun treatment data
treatments <- read_csv("Data/treatments_meta.csv")
#merging them
mortalitydatamerge <- full_join(mortalitydataraw, treatments)
head(mortalitydatamerge)
#Assigning Neighbour values to C, E and T NAs
#C (competition), E (environment) and T (trait) subplots
#C has neighbours (1), E and T do not (0)
# If Neighbour == NA, then make C_E_or_T Cs 1, Es and Ts 0
#Test: 2 A 9, 1 VERO T1 has Neighbours (1) --> should stay one
mortalitydata1 <- within(mortalitydatamerge, Neighbours[C_E_or_T == 'C' & is.na(Neighbours)] <- '1')
mortalitydata2 <- within(mortalitydata1, Neighbours[C_E_or_T == 'T' & is.na(Neighbours)] <- '0')
mortalitydata <- within(mortalitydata2, Neighbours[C_E_or_T == 'E' & is.na(Neighbours)] <- '0')

#### Removing some subplots
#Removing these due to death from extraneous factors (e.g. kangaroo, me):
#Site 1 C TROR C1, Site 1 A VERO E 3, Site 1 A LARO E3, Site 1 A PLDE C1, Site 1 C ARCA C 2, Site 4 B ARCA E 3
# Site 2 C PLDE T 1, Site 2 C TROR T 1, Site 6 A TROR E 2 and Site 8 B TRCY E 2 (keeping this last one)
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "C" & Species == "TROR" & C_E_or_T == "C" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "A" & Species == "VERO" & C_E_or_T == "E" & Rep == "3"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "A" & Species == "LARO" & C_E_or_T == "E" & Rep == "3"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "A" & Species == "PLDE" & C_E_or_T == "C" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "C" & Species == "ARCA" & C_E_or_T == "C" & Rep == "2"))
mortalitydata <- mortalitydata %>% filter(!(Site == "2" & Plot == "C" & Species == "PLDE" & C_E_or_T == "T" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "2" & Plot == "C" & Species == "TROR" & C_E_or_T == "T" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "6" & Plot == "A" & Species == "TROR" & C_E_or_T == "E" & Rep == "2"))
mortalitydata <- mortalitydata %>% filter(!(Site == "4" & Plot == "B" & Species == "ARCA" & C_E_or_T == "E" & Rep == "3"))

#Removing these due to it not being (or can't tell) my focal:
# Site 2 C HYGL E 1, Site 5 A ARCA C2, Site 5 B POLE E1 and Site 5 B POLE E 2, site 6 A VERO C 3, 
#Site 7 B POLE C 1, Site 1 B POLE E 2, Site 1 A HYGL  E 1, Site 1 C VERO E 3, Site 3 B TROR C1, Site 3 C ARCA C3
mortalitydata <- mortalitydata %>% filter(!(Site == "2" & Plot == "C" & Species == "HYGL" & C_E_or_T == "E" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "5" & Plot == "A" & Species == "ARCA" & C_E_or_T == "C" & Rep == "2"))
mortalitydata <- mortalitydata %>% filter(!(Site == "5" & Plot == "B" & Species == "POLE" & C_E_or_T == "E" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "5" & Plot == "B" & Species == "POLE" & C_E_or_T == "E" & Rep == "2"))
mortalitydata <- mortalitydata %>% filter(!(Site == "6" & Plot == "A" & Species == "VERO" & C_E_or_T == "C" & Rep == "3"))
mortalitydata <- mortalitydata %>% filter(!(Site == "7" & Plot == "B" & Species == "POLE" & C_E_or_T == "C" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "B" & Species == "POLE" & C_E_or_T == "E" & Rep == "2"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "A" & Species == "HYGL" & C_E_or_T == "E" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "1" & Plot == "C" & Species == "VERO" & C_E_or_T == "E" & Rep == "3"))
mortalitydata <- mortalitydata %>% filter(!(Site == "3" & Plot == "B" & Species == "TROR" & C_E_or_T == "C" & Rep == "1"))
mortalitydata <- mortalitydata %>% filter(!(Site == "3" & Plot == "C" & Species == "ARCA" & C_E_or_T == "C" & Rep == "3"))
mortalitydata <- mortalitydata %>% filter(!(Site == "8" & Plot == "B" & Species == "POLE" & C_E_or_T == "E" & Rep == "3"))


### Updating this based on more accurate Neighbours data in June 2021, from seed_count_2020 file
#WORK IN PROGRESS
#mortalitydata %>% select(Site, Plot, Species, C_E_or_T, Rep, No_germinated, Survival)

#Check this:
#Should I assign Site 1 C TROR C 1 needs to be assigned a 0 for survival? Died 

#Filtering by Species

arcamortality <- mortalitydata %>% filter(Species == "ARCA")
hyglmortality <- mortalitydata %>% filter(Species == "HYGL")
laromortality <- mortalitydata %>% filter(Species == "LARO")
peaimortality <- mortalitydata %>% filter(Species == "PEAI")
pldemortality <- mortalitydata %>% filter(Species == "PLDE")
polemortality <- mortalitydata %>% filter(Species == "POLE")
trcymortality <- mortalitydata %>% filter(Species == "TRCY")
trormortality <- mortalitydata %>% filter(Species == "TROR")
veromortality <- mortalitydata %>% filter(Species == "VERO")

#### MODELLING DATA
# Mixed effects binomial model
#modelallsp1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
#                     family = binomial,
#                     data = mortalitydata)
#summary(modelallsp1)

#Modelling by individual species
#ARCA
arcamodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                   family = binomial,
                   data = arcamortality)
summary(arcamodel1)

arcamodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel2)

arcamodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel3)

arcamodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel4)

arcamodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel5)

#########25/03/21, testing simpler models (removing fixed effects):

arcamodel6 <- glmer(Survival ~ Neighbours + Treatment + (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel6)

AIC(arcamodel5, arcamodel6)

arcamodel7 <- glmer(Survival ~ Neighbours + Cover +  (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel7)

AIC(arcamodel7, arcamodel5)
#Model 7 seemingly a better fit

arcamodel8 <- glmer(Survival ~ Neighbours + (1|Site/Plot),
                    family = binomial,
                    data = arcamortality)
summary(arcamodel8)
AIC(arcamodel8, arcamodel5)
#Model 8 seemingly a better fit again...

############ End of testing removing fixed effects

AIC(arcamodel1, arcamodel2, arcamodel3, arcamodel4, arcamodel5)
#Models 3, 4 and 5 have similar fit. Model 5 most parsimonious. Significant Neighbours

#PLOT
#Not working yet: #arcamortality$Treatment %>% rename_if(Control, Ambient)
arcamortality$Treatment <- factor(arcamortality$Treatment, level = c("Dry", "Control", "Wet"))

#note that mod2 = cover splits represents shade and sun separately
cat_plot(arcamodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = arcamortality, main = "ARCA")

#HYGL
hyglmodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = hyglmortality)
summary(hyglmodel1)

hyglmodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = hyglmortality)
summary(hyglmodel2)

hyglmodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = hyglmortality)
summary(hyglmodel3)

hyglmodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = hyglmortality)
summary(hyglmodel4)

hyglmodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = hyglmortality)
summary(hyglmodel5)

AIC(hyglmodel1, hyglmodel2, hyglmodel3, hyglmodel4, hyglmodel5)
#Models 1-4 failing to converge.
# Models 3 and 4 have similarly lowest fits. Model 4 most parsimonious. Significant negative treatmentDry
# and positive TreatmentDry:CoverSun

#PLOT
hyglmortality$Treatment <- factor(hyglmortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(hyglmodel1, pred = Treatment, modx = Neighbours, data = hyglmortality)
cat_plot(hyglmodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = hyglmortality, main = "HYGL")

#LARO
laromodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = laromortality)
summary(laromodel1)

laromodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = laromortality)
summary(laromodel2)

laromodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = laromortality)
summary(laromodel3)

laromodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = laromortality)
summary(laromodel4)

laromodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = laromortality)
summary(laromodel5)

AIC(laromodel1, laromodel2, laromodel3, laromodel4, laromodel5)
#Model 5 most parsimonious. Significant effect of neighbours.

#PLOT
laromortality$Treatment <- factor(laromortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(laromodel1, pred = Treatment, modx = Neighbours, data = laromortality)
cat_plot(laromodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = laromortality, main = "LARO")

#PEAI
peaimodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = peaimortality)
summary(peaimodel1)

peaimodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = peaimortality)
summary(peaimodel2)

peaimodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = peaimortality)
summary(peaimodel3)

peaimodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = peaimortality)
summary(peaimodel4)

peaimodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = peaimortality)
summary(peaimodel5)

AIC(peaimodel1, peaimodel2, peaimodel3, peaimodel4, peaimodel5)
#Model 5 most parsimonious. Significant intercept and TreatmentDry

#PLOT
peaimortality$Treatment <- factor(peaimortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(peaimodel1, pred = Treatment, modx = Neighbours, data = peaimortality)
cat_plot(peaimodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = peaimortality, main = "PEAI")

#PLDE
pldemodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = pldemortality)
summary(pldemodel1)

pldemodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = pldemortality)
summary(pldemodel2)

pldemodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = pldemortality)
summary(pldemodel3)

pldemodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = pldemortality)
summary(pldemodel4)

pldemodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = pldemortality)
summary(pldemodel5)

AIC(pldemodel1, pldemodel2, pldemodel3, pldemodel4, pldemodel5)
#Models 4 and 5 have similar fit. Model 5 most parsimonious. No significance.

#PLOT
pldemortality$Treatment <- factor(pldemortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(pldemodel1, pred = Treatment, modx = Neighbours, data = pldemortality)
cat_plot(pldemodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = pldemortality, main = "PLDE")

#POLE
polemodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = polemortality)
summary(polemodel1)

polemodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = polemortality)
summary(polemodel2)

polemodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = polemortality)
summary(polemodel3)

polemodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = polemortality)
summary(polemodel4)

polemodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = polemortality)
summary(polemodel5)

AIC(polemodel1, polemodel2, polemodel3, polemodel4, polemodel5)
#Model 5 has best fit but they aren't converging properly. No significance

#PLOT
polemortality$Treatment <- factor(polemortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(polemodel1, pred = Treatment, modx = Neighbours, data = polemortality)
cat_plot(polemodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = polemortality, main = "POLE")

#TRCY
trcymodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = trcymortality)
summary(trcymodel1)

trcymodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = trcymortality)
summary(trcymodel2)

trcymodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = trcymortality)
summary(trcymodel3)

trcymodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = trcymortality)
summary(trcymodel4)

trcymodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = trcymortality)
summary(trcymodel5)

AIC(trcymodel1, trcymodel2, trcymodel3, trcymodel4, trcymodel5)
#Model 3 best fit. Significant TreatmentWet:CoverSun and Neighbours1:CoverSun

#PLOT
trcymortality$Treatment <- factor(trcymortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(trcymodel1, pred = Treatment, modx = Neighbours, data = trcymortality)
cat_plot(trcymodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = trcymortality, main = "TRCY")

#TROR
trormodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = trormortality)
summary(trormodel1)

trormodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = trormortality)
summary(trormodel2)

trormodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = trormortality)
summary(trormodel3)

trormodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = trormortality)
summary(trormodel4)

trormodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = trormortality)
summary(trormodel5)

AIC(trormodel1, trormodel2, trormodel3, trormodel4, trormodel5)
#All models have similar fit but 1-3 didn't converge. Model 5 most parsimonious. No significance.

#PLOT
trormortality$Treatment <- factor(trormortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(trormodel1, pred = Treatment, modx = Neighbours, data = trormortality)
cat_plot(trormodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = trormortality, main = "TROR")

#VERO
veromodel1 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^3 + (1|Site/Plot),
                    family = binomial,
                    data = veromortality)
summary(veromodel1)

veromodel2 <- glmer(formula = Survival ~ (Neighbours + Treatment + Cover)^2 + (1|Site/Plot),
                    family = binomial,
                    data = veromortality)
summary(veromodel2)

veromodel3 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + Neighbours:Cover + (1|Site/Plot),
                    family = binomial,
                    data = veromortality)
summary(veromodel3)

veromodel4 <- glmer(formula = Survival ~ Neighbours + Treatment*Cover + (1|Site/Plot),
                    family = binomial,
                    data = veromortality)
summary(veromodel4)

veromodel5 <- glmer(formula = Survival ~ Neighbours + Treatment + Cover + (1|Site/Plot),
                    family = binomial,
                    data = veromortality)
summary(veromodel5)

AIC(veromodel1, veromodel2, veromodel3, veromodel4, veromodel5)
#Model 5 most parsimonious. Significant intercept and Neighbours. Models 1 - 4 didn't converge

#PLOT
veromortality$Treatment <- factor(veromortality$Treatment, level = c("Dry", "Control", "Wet"))
cat_plot(veromodel1, pred = Treatment, modx = Neighbours, data = veromortality)
cat_plot(veromodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = veromortality, main = "VERO")

#Plotting all of the above in a 3x3 grid
#### Note that this isn't working!
#Remove mod2 = Cover if I don't wan't them split between shade and sun

#par(mfrow=c(3,3))
cat_plot(arcamodel1, pred = Treatment, modx = Neighbours, data = arcamortality, main = "ARCA")
cat_plot(hyglmodel1, pred = Treatment, modx = Neighbours, data = hyglmortality, main = "HYGL")
cat_plot(laromodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = laromortality, main = "LARO")
cat_plot(peaimodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = peaimortality, main = "PEAI")
cat_plot(pldemodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = pldemortality, main = "PLDE")
cat_plot(polemodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = polemortality, main = "POLE")
cat_plot(trcymodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = trcymortality, main = "TRCY")
cat_plot(trormodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = trormortality, main = "TROR")
cat_plot(veromodel1, pred = Treatment, modx = Neighbours, mod2= Cover, data = veromortality, main = "VERO")

#Could do below but it looks terrible, repeating y and x axis labels and key
#library(gridExtra)
#grid.arrange(arcaplot, hyglplot, nrow = 1)

#### When did things die? ####
##Note that I need dates to be recognised as dates
library(lubridate)
mortalitydata$Date_mortality <- as.Date(mortalitydata$Date_mortality, "%d/%m/%y")
str(mortalitydata)

par(mfrow=c(2,1))
hist(mortalitydata$Date_mortality, breaks = "days")
hist(mortalitydata$Date_mortality, breaks = "weeks")

####### Coming back to this in March 2021 
#Want to see how survival varied across sites
survivalsite <- mortalitydata %>% group_by(Site) %>% filter(Survival == "1") %>% tally()
ggplot(survivalsite, aes(x = Site, y = n))+
  geom_point()+
  theme_classic()
survivalplot <- mortalitydata %>% group_by(Site, Plot) %>% filter(Survival == "1") %>% tally()
plot1 <- ggplot(survivalplot, aes(x = Site, y = n))+
  geom_point()+
  theme_classic()+
  ylab("Survival count")
#With data from germination_stocktake
plot2 <- ggplot(tallyplotgerm, aes(x = Site, y = n))+
  geom_point()+
  theme_classic()+
  ylab("Germination count")
gridExtra::grid.arrange(plot1, plot2)
#Can also look at mortality count by site

# Want to calculate survival rates for each species overall and split across 
# sun and shade
survivalcounts <- mortalitydata %>% group_by(Species, Cover, Survival) %>% tally()
survivalrates <- survivalcounts %>% filter(!is.na(Survival)) %>%
  group_by(Species, Cover) %>%
  mutate(number_survived = sum(n[Survival == "1"]),
            total_number = sum(n[Survival == "1"])+sum(n[Survival == "0"]),
            survival_rate = number_survived/total_number)
#Filter to one row
survivalcover <- survivalrates %>% group_by(Species, Cover) %>% filter(row_number() == 1)
#Overall
survivalspecies <- survivalcounts %>% filter(!is.na(Survival)) %>%
  group_by(Species) %>%
  mutate(number_survived = sum(n[Survival == "1"]),
         total_number = sum(n[Survival == "1"])+sum(n[Survival == "0"]),
         survival_rate = number_survived/total_number) %>%
  filter(row_number() == 1)
#Selecting only relevant columns
survivalsp <- survivalspecies %>% select(Species, survival_rate)

#Is survival rate related to SLA, LDMC or D13?
alltraits <- merge(speciestraits, summaryD13C)
survivaltraits <- merge(survivalsp, alltraits)

#Pulling out for survival for vero, looking by site and by plot
verosurvival <- mortalitydata %>% filter(Species == "VERO") %>% 
  group_by(Site, Plot, Survival) %>% tally()
veroplot <- verosurvival %>% filter(!is.na(Survival)) %>%
  group_by(Site, Plot) %>%
  mutate(number_survived = sum(n[Survival == "1"]),
         total_number = sum(n[Survival == "1"])+sum(n[Survival == "0"]),
         survival_rate = number_survived/total_number) %>%
  filter(row_number() == 1)

# Want to calculate survival rates by species across a gradient of canopy cover.
#Calculating suvivalrates by plot
#Merging mortality data with canopy cover continuous

canopylitterdata <- read_csv("Data/canopy_litter_2020.csv")
canopydata <- canopylitterdata  %>% rowwise() %>%
  mutate(meancc = mean(c(`Measurement 1`, `Measurement 2`, `Measurement 3`, `Measurement 4`, `Measurement 5`)),
         cc_percentage = meancc/24*100)
canopydatatrim <- canopydata %>% select(Site, Plot, cc_percentage)
mortalitydatacanopy <- merge(mortalitydata, canopydatatrim)

survivalcountsplot <- mortalitydatacanopy %>% group_by(Species, cc_percentage, Survival) %>% tally()
survivalratesplot <- survivalcountsplot %>% filter(!is.na(Survival)) %>%
  group_by(Species, cc_percentage) %>%
  mutate(number_survived = sum(n[Survival == "1"]),
         total_number = sum(n[Survival == "1"])+sum(n[Survival == "0"]),
         survival_rate = number_survived/total_number)
#Filter to one row
survivalcanopy <- survivalratesplot %>% group_by(Species, cc_percentage) %>% filter(row_number() == 1)

ggplot(survivalcanopy, aes(x = log(cc_percentage+1), y = survival_rate))+
  geom_point(alpha = 0.4)+
  theme_classic()+
  geom_smooth(method="lm")+
  facet_wrap(~Species)

#Can't have site within plot as a random effect without reworking this dataframe
for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      lm(survival_rate ~ log(cc_percentage+1),
           data = filter(survivalcanopy, Species == specieslist[i]))))
}

##### October 2021 - data from data_preparation file ####
source("data_preparation.R")
## Calculating rates of survival to produce at least one viable seed
survivalcount <- dataall %>% group_by(Species, ProducedSeeds) %>% tally()
survivalsp <- survivalcount %>%
  group_by(Species) %>%
  mutate(number_survived = sum(n[ProducedSeeds == "1"]),
         total_number = sum(n[ProducedSeeds == "1"])+sum(n[ProducedSeeds == "0"]),
         survival_rate = number_survived/total_number) %>%
  filter(row_number() == 1)
#Selecting only relevant columns
survivalsp <- survivalsp %>% select(Species, survival_rate)
#####And by shade and sun
survivalshadecount <- dataall %>% group_by(Species, Cover, ProducedSeeds) %>% tally()
survivalshadesp <- survivalshadecount %>%
  group_by(Species, Cover) %>%
  mutate(number_survived = sum(n[ProducedSeeds == "1"]),
         total_number = sum(n[ProducedSeeds == "1"])+sum(n[ProducedSeeds == "0"]),
         survival_rate = number_survived/total_number) %>%
  filter(row_number() == 1)
#Selecting only relevant columns
survivalshadesp <- survivalshadesp %>% select(Species, Cover, survival_rate)





## Unimportant - old analysis ####
#MY WORK/IDEAS BEFORE MEETING WITH JOHN

#Creating dataframe for counts of each group and by species
countall <- mortalitydata %>% drop_na %>% group_by(Species) %>% count(Survival, Neighbours)
head(countall)
#Making survival and Neighbours factors for plotting purposes
countall <- countall %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))

str(countall)

#ARCA
countarca <- arcamortality %>% drop_na %>% count(Survival, Neighbours)
countarca <- countarca %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))

#And for each watering treatment
#Want mortality as a percentage/proportion for each plot, then averaging that
countarcaplot <- arcamortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
# Number that died without neighbours / total for that plot with neighbours
# and number that died with neighbours/total for that plot with neighbours
# Then bring in the treatment info
countarcanb <- countarcaplot %>% filter(Neighbours == "1")
countarcanbmortality <- countarcanb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
arcamortnb <- full_join(countarcanbmortality, treatments)
countarcawonb <- countarcaplot %>% filter(Neighbours == "0")
countarcawonbmortality <- countarcawonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
arcamortwonb <- full_join(countarcawonbmortality, treatments)
arcamortperc <- full_join(arcamortnb, arcamortwonb) %>% drop_na()

#HYGL
counthygl <- hyglmortality %>% drop_na %>% count(Survival, Neighbours)
counthygl <- counthygl %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
counthyglplot <- hyglmortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
counthyglnb <- counthyglplot %>% filter(Neighbours == "1")
counthyglnbmortality <- counthyglnb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
hyglmortnb <- full_join(counthyglnbmortality, treatments)
counthyglwonb <- counthyglplot %>% filter(Neighbours == "0")
counthyglwonbmortality <- counthyglwonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
hyglmortwonb <- full_join(counthyglwonbmortality, treatments)
hyglmortperc <- full_join(hyglmortnb, hyglmortwonb) %>% drop_na()

#LARO
countlaro <- laromortality %>% drop_na %>% count(Survival, Neighbours)
countlaro <- countlaro %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
countlaroplot <- laromortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
countlaronb <- countlaroplot %>% filter(Neighbours == "1")
countlaronbmortality <- countlaronb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
laromortnb <- full_join(countlaronbmortality, treatments)
countlarowonb <- countlaroplot %>% filter(Neighbours == "0")
countlarowonbmortality <- countlarowonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
laromortwonb <- full_join(countlarowonbmortality, treatments)
laromortperc <- full_join(laromortnb, laromortwonb) %>% drop_na()

#PEAI
countpeai <- peaimortality %>% drop_na %>% count(Survival, Neighbours)
countpeai <- countpeai %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
countpeaiplot <- peaimortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
countpeainb <- countpeaiplot %>% filter(Neighbours == "1")
countpeainbmortality <- countpeainb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
peaimortnb <- full_join(countpeainbmortality, treatments)
countpeaiwonb <- countpeaiplot %>% filter(Neighbours == "0")
countpeaiwonbmortality <- countpeaiwonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
peaimortwonb <- full_join(countpeaiwonbmortality, treatments)
peaimortperc <- full_join(peaimortnb, peaimortwonb) %>% drop_na()

#PLDE
countplde <- pldemortality %>% drop_na %>% count(Survival, Neighbours)
countplde <- countplde %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
countpldeplot <- pldemortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
countpldenb <- countpldeplot %>% filter(Neighbours == "1")
countpldenbmortality <- countpldenb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
pldemortnb <- full_join(countpldenbmortality, treatments)
countpldewonb <- countpldeplot %>% filter(Neighbours == "0")
countpldewonbmortality <- countpldewonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
pldemortwonb <- full_join(countpldewonbmortality, treatments)
pldemortperc <- full_join(pldemortnb, pldemortwonb) %>% drop_na()

#POLE
countpole <- polemortality %>% drop_na %>% count(Survival, Neighbours)
countpole <- countpole %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
countpoleplot <- polemortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
countpolenb <- countpoleplot %>% filter(Neighbours == "1")
countpolenbmortality <- countpolenb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
polemortnb <- full_join(countpolenbmortality, treatments)
countpolewonb <- countpoleplot %>% filter(Neighbours == "0")
countpolewonbmortality <- countpolewonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
polemortwonb <- full_join(countpolewonbmortality, treatments)
polemortperc <- full_join(polemortnb, polemortwonb) %>% drop_na()

#TRCY
counttrcy <- trcymortality %>% drop_na %>% count(Survival, Neighbours)
counttrcy <- counttrcy %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
counttrcyplot <- trcymortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
counttrcynb <- counttrcyplot %>% filter(Neighbours == "1")
counttrcynbmortality <- counttrcynb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
trcymortnb <- full_join(counttrcynbmortality, treatments)
counttrcywonb <- counttrcyplot %>% filter(Neighbours == "0")
counttrcywonbmortality <- counttrcywonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
trcymortwonb <- full_join(counttrcywonbmortality, treatments)
trcymortperc <- full_join(trcymortnb, trcymortwonb) %>% drop_na()

#TROR
counttror <- trormortality %>% drop_na %>% count(Survival, Neighbours)
counttror <- counttror %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
counttrorplot <- trormortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
counttrornb <- counttrorplot %>% filter(Neighbours == "1")
counttrornbmortality <- counttrornb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
trormortnb <- full_join(counttrornbmortality, treatments)
counttrorwonb <- counttrorplot %>% filter(Neighbours == "0")
counttrorwonbmortality <- counttrorwonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
trormortwonb <- full_join(counttrorwonbmortality, treatments)
trormortperc <- full_join(trormortnb, trormortwonb) %>% drop_na()

#VERO
countvero <- veromortality %>% drop_na %>% count(Survival, Neighbours)
countvero <- countvero %>% mutate(Survival = as.factor(Survival)) %>% mutate(Neighbours = as.factor(Neighbours))
countveroplot <- veromortality %>% drop_na %>% group_by(Site, Plot) %>% count(Survival, Neighbours)
countveronb <- countveroplot %>% filter(Neighbours == "1")
countveronbmortality <- countveronb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
veromortnb <- full_join(countveronbmortality, treatments)
countverowonb <- countveroplot %>% filter(Neighbours == "0")
countverowonbmortality <- countverowonb %>% group_by(Site, Plot) %>% mutate(Mortality = n / sum(n)*100) %>% 
  mutate(n = sum(n)) %>% filter(Survival == "0") %>% select(-Survival)
veromortwonb <- full_join(countverowonbmortality, treatments)
veromortperc <- full_join(veromortnb, veromortwonb) %>% drop_na()

#### PLOTTING DATA

#Stacked barchart
ggplot(countall, aes(fill=Survival, y=n, x=Neighbours)) + 
  geom_bar(position="stack", stat="identity")

#Percentage stacked barchart
ggplot(countarca, aes(fill=Survival, y=n, x=Neighbours)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("ARCA")

ggplot(counthygl, aes(fill=Survival, y=n, x=Neighbours)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("HYGL")

#ARCA
# Boxplots for wet, control and dry plots
# Reording x-axis control, dry and wet
arcamortperc$Treatment <- factor(arcamortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(arcamortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw() +
  ggtitle("ARCA")

#HYGL
hyglmortperc$Treatment <- factor(hyglmortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(hyglmortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("HYGL")

#LARO
laromortperc$Treatment <- factor(laromortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(laromortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("LARO")

#PEAI
peaimortperc$Treatment <- factor(peaimortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(peaimortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("PEAI")

#PLDE
pldemortperc$Treatment <- factor(pldemortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(pldemortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("PLDE")

#POLE
polemortperc$Treatment <- factor(polemortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(polemortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("POLE")

#TRCY
trcymortperc$Treatment <- factor(trcymortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(trcymortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("TRCY")

#TROR
trormortperc$Treatment <- factor(trormortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(trormortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("TROR")

#VERO
veromortperc$Treatment <- factor(veromortperc$Treatment, level = c("Dry", "Control", "Wet"))
ggplot(veromortperc, aes(y = Mortality, x = Treatment))+
  geom_boxplot(aes(fill = Neighbours))+
  theme_bw()+
  ggtitle("VERO")
