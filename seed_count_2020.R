####### Seed prodution data
### May 2021
#### WA 2020 season data

library(tidyverse) 
library(lme4)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(gridExtra)
library(vegan)

my_theme <- theme(axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16),
                  axis.text = element_text(size = 16),
                  strip.text.x = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.title = element_text(size = 14))

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=seq(min(x, na.rm=T), max(x, na.rm=T), length.out=7))
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y, na.rm=T)
  rect(breaks[-nB], 0, breaks[-1], y,...)
}

###########Importing and cleaning data
seeddataraw <- read_csv("Data/seed_count_2020.csv")
head(seeddataraw)

#I also need to sum the number of seeds for different envelopes of same individual
#Example to test this is VERO 6C E1 - has 10 + 15 viable seeds, Site 1 C VERO C1, 10 + 17 viable, 0 +8 inviable; 
#2A VERO E1, 7 + 12 + 8 viable, 3 +0 +1 inviable; 
seeddatarawgrouped <- seeddataraw %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(No_viable_seeds_grouped = sum(No_viable_seeds), 
         No_inviable_seeds_grouped = sum(No_inviable_seeds)) %>% filter(row_number() == 1)
## Trying filter by top row thing -- this works!!
#test <- seeddatarawgrouped %>% filter(Site == "1", Plot == "C", Species == "VERO")

##Merging with information about treatments
treatments <- read_csv("Data/treatments_meta.csv")
seeddata <- merge(treatments, seeddatarawgrouped)
head(seeddata)
#Reordering in order of dry, control, wet.
level = c("Dry", "Control", "Wet")
seeddata$Treatment <- factor(seeddata$Treatment, level = c("Dry", "Control", "Wet"))
#Need Site, Treatment and Species to be factors
seeddata$Site <- factor(seeddata$Site)
seeddata$Treatment <- factor(seeddata$Treatment)
seeddata$Species <- factor(seeddata$Species)
str(seeddata)

# I think with requirements a continuous variable to be plotted
with(seeddata, pairs(No_viable_seeds_grouped ~ (No_inviable_seeds_grouped), diag.panel = panel.hist))
with(seeddata, pairs(log(No_viable_seeds_grouped+1) ~ log(No_inviable_seeds_grouped+1), diag.panel = panel.hist))

#Can't take a log of zero, so plus one then logging
#Note that log(1) = 0, so zeros on plots are true zeros
seeddata <- seeddata %>% mutate(log1_viable_seeds = log(No_viable_seeds_grouped+1))
glimpse(seeddata)

#Plots with neighbours - C. Plots without neighbours - E and T
seeddata$Neighbours <- ifelse(seeddata$C_E_or_T == "C", "Yes", "No")
str(seeddata)
seeddata$Neighbours <- factor(seeddata$Neighbours)

#Removing Site 1 B POLE E 2 since I wasn't sure it was my focal
seeddata <- seeddata %>% filter(!(Site == "1" & Plot == "B" & Species == "POLE" & C_E_or_T == "E" & Rep == "2"))


###Merging with community survey data
surveydata <- read_csv("Data/community_surveys_data.csv")
glimpse(surveydata)
surveydata$Site <- factor(surveydata$Site)
surveydata$Species <- factor(surveydata$Species)
seedsurveys <- full_join(seeddata, surveydata)
glimpse(seedsurveys)
#Replacing Neighbour_count NAs where Neighbour = Yes with 0. e.g. so Dodder now is counted as 0
# This is because I can't tally with NAs and treating Dodder as 0 because it is in my no-comp. subplots
seedsurveys <- within(seedsurveys, Neighbour_count[Neighbour_sp == 'Dodder'] <- '0')
seedsurveys$Neighbour_count <- as.numeric(seedsurveys$Neighbour_count)
#Renaming E and T Neighbours to Yes that actually have neighbours (neighbour_count >= 1)
seedsurveys$Neighbours <- ifelse(seedsurveys$Neighbour_count >= 1, "Yes", "No")
#Also want to plot the subplots with 0 neighbour_count, so need NAs to be 0s
#This filters to subplots without Neighbour_count (currently NAs, same with Neighbours)
testing <- seedsurveys %>% filter(!complete.cases(Neighbour_count)) 
testing$Neighbour_count <- "0"
testing$Neighbours <- "No"
testing2 <- seedsurveys %>% filter(complete.cases(Neighbour_count))
#Checking numbers are right - 1740 lines in seedsurveys, 633 !complete.cases, 1107 complete.cases. Phew!
testing$Neighbour_count <- as.numeric(testing$Neighbour_count)
seedsurveysfixed <- full_join(testing, testing2)

#Removing Site 1 B POLE E 2 since I wasn't sure it was my focal
seedsurveysfixed <- seedsurveysfixed %>% filter(!(Site == "1" & Plot == "B" & Species == "POLE" & C_E_or_T == "E" & Rep == "2"))

#I shouldn't have to do the below in two steps but I can't figure out one, and it works
neighbourabundance <- seedsurveysfixed %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  summarise(Total_abundance = sum(Neighbour_count))
seedabundance <- merge(neighbourabundance, seedsurveysfixed)

# #Number of plots with Dodder overall and that have no other neighbours
# # None! So I must have Dodder recorded separately...
# seedsurveys %>% tally(Neighbour_sp == "Dodder" & Neighbours == "Yes")
# seedsurveys %>% tally(Neighbour_sp == "Dodder")

#Some Species names to be fixed - PLDE' is actually a VERO.
seedabundance %>% filter(Species == "PLDE'")
#Changing the rep number to 9 (so it is unique) and species name to VERO.
seedabundance <- within(seedabundance, Rep[Species == "PLDE'"] <- '9')
seedabundance <- within(seedabundance, Species[Species == "PLDE'"] <- 'VERO')
# PLDE-POLE is actually POLE, changing rep number to 9 and species name to POLE
seedabundance %>% filter(Species == "PLDE-POLE")
seedabundance <- within(seedabundance, Rep[Species == "PLDE-POLE"] <- '9')
seedabundance <- within(seedabundance, Species[Species == "PLDE-POLE"] <- 'POLE')
#NAs, site 3, plot  - 756 - 759 are NAs. Manually changing them to Unks and removing for now
## Will check notes later
seedabundance %>% filter(Species == "UNK")
seedabundance <- seedabundance %>% filter(Species != "UNK")

# Separating abundances of intra- and inter-specific neighbours
#Isaac's code first two lines
seedabundance$Matching <- seedabundance$Species == seedabundance$Neighbour_sp
seedabundance$Matching <- ifelse(seedabundance$Matching == TRUE, 1, 0)

seedabundance <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  mutate(Intra_abundance = sum(Neighbour_count[Matching == "1"]))
seedabundance <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  mutate(Inter_abundance = sum(Neighbour_count[Matching == "0"]))

#Turning into a dataset with one row per subplot
allabundanceseeds <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep)%>%
                      filter(row_number() == 1)
### Creating a metric of diversity of neighbours per subplot
#Three rows exactly the same (89-90, 465-466 and 1436-1437) so removing one of each
#1apldee3, 2barcae2, 7atrore2, 

test <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep, log1_viable_seeds, 
                                   Neighbour_sp, Neighbour_count) %>%
                  filter(row_number() == 1)
#Note that there is a problem in above code, only expect it to remove 3 lines but removes 8
mf<- test %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  pivot_longer(names_to = "Neighbour_sp", values_to = "Neighbour_count")
###Didn't fix this stuff. Not working
#diversity(seedabundance$Neighbour_sp, index = "shannon")


######### Note that seeddata is not as cleaned as allabundanceseeds
#Filtering allabundanceseeds to things that produced seeds
seeddataclean <- allabundanceseeds %>% filter(!is.na(No_viable_seeds))

#Creating data sets by Species
seedarca <- seeddataclean %>% filter(Species == "ARCA")
seedhygl <- seeddataclean %>% filter(Species == "HYGL")
seedlaro <- seeddataclean %>% filter(Species == "LARO")
seedpeai <- seeddataclean %>% filter(Species == "PEAI")
seedplde <- seeddataclean %>% filter(Species == "PLDE")
seedpole <- seeddataclean %>% filter(Species == "POLE")
seedtrcy <- seeddataclean %>% filter(Species == "TRCY")
seedtror <- seeddataclean %>% filter(Species == "TROR")
seedvero <- seeddataclean %>% filter(Species == "VERO")

#Testing how to do this with a loop - works but outputs capitals - e.g. seedARCA
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
#This works but calls i 1-9 so seed1, seed 2... seed 9.
for (i in 1:length(specieslist)){
  nam <- paste0("seed", specieslist[i])
  assign(nam, filter(seeddata, Species == specieslist[i]))
}
#Can use either paste("seed", i, sep = "") or paste0("seeds", i) - same result

#Counts~
#Count of sample size by species
countsp <- seeddataclean %>% group_by(Species) %>% count()
countsp
#Count of number of plants with at least one seed per species
countseedssp <- seeddataclean %>% filter(No_viable_seeds_grouped > 0) %>% group_by(Species) %>% count()
countseedssp
differencecount <- countsp$n - countseedssp$n
differencecount
#Sample size by Site and by Species
countsite <- seeddata %>% group_by(Site) %>% count()
countsiteseeds <- seeddata %>% filter(No_viable_seeds_grouped > 0) %>% group_by(Site) %>% count()
countsitesp <- seeddata %>% filter(No_viable_seeds_grouped > 0) %>% group_by(Site, Species) %>% count()
#~Counts

#Plotting main relationships of interest
#Site~
ggplot(seeddata, aes(x = Site, y = log1_viable_seeds, group = Site))+
  geom_boxplot()+
  geom_jitter(alpha = 0.3)+
 # ylab("log(Number of viable seeds)")+
  theme_classic()

#Site model
modsite <- aov(log1_viable_seeds ~ Site, seeddata)
summary(modsite)
TukeyHSD(modsite)
#Suggests that Site 1, 2, 3, 7 and 8 are s.d. from Site 4
# Site 6 is s.d. from Site 1 and Site 8
# Site 8 s.d. from Site 5
# Site 4 has lowest viability, alongside 5/6
#~Site

#Cover~
ggplot(seeddata, aes(x = Cover, y = log1_viable_seeds))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()

t.test(log1_viable_seeds ~ Cover, data = seeddata)

ggplot(seeddata, aes(x = Cover, y = log1_viable_seeds))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species))

t.test(log1_viable_seeds ~ Cover, data = seedarca)
t.test(log1_viable_seeds ~ Cover, data = seedhygl)
t.test(log1_viable_seeds ~ Cover, data = seedlaro)
t.test(log1_viable_seeds ~ Cover, data = seedpeai)
t.test(log1_viable_seeds ~ Cover, data = seedplde)
t.test(log1_viable_seeds ~ Cover, data = seedpole)
t.test(log1_viable_seeds ~ Cover, data = seedtrcy)
t.test(log1_viable_seeds ~ Cover, data = seedtror)
t.test(log1_viable_seeds ~ Cover, data = seedvero)

#Testing how to do this with a loop -- it works!!
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
  
  for (i in 1:length(specieslist)){
    print(specieslist[i])
    print(
      t.test(log1_viable_seeds ~ Cover, data = filter(seeddata, Species == specieslist[i])))
  }
#~Cover

#Watering treatment~
ggplot(seeddata, aes(x = Treatment, y = log1_viable_seeds))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()

treatmentmod1 <- aov(log1_viable_seeds ~ Treatment, seeddata)
summary(treatmentmod1)
#bartlett.test(log1_viable_seeds ~ Treatment, seeddata)
TukeyHSD(treatmentmod1)

ggplot(seeddata, aes(x = Treatment, y = log1_viable_seeds))+
  geom_boxplot()+
  geom_jitter(alpha = 0.4)+
  theme_classic()+
  facet_wrap(vars(Species))

treatmentmodarca <- aov(log1_viable_seeds ~ Treatment, seedarca)
summary(treatmentmodarca)
treatmentmodhygl <- aov(log1_viable_seeds ~ Treatment, seedhygl)
summary(treatmentmodhygl)
treatmentmodlaro <- aov(log1_viable_seeds ~ Treatment, seedlaro)
summary(treatmentmodlaro)
treatmentmodpeai <- aov(log1_viable_seeds ~ Treatment, seedpeai)
summary(treatmentmodpeai)
treatmentmodplde <- aov(log1_viable_seeds ~ Treatment, seedplde)
summary(treatmentmodplde)
TukeyHSD(treatmentmodplde)
treatmentmodpole <- aov(log1_viable_seeds ~ Treatment, seedpole)
summary(treatmentmodpole)
treatmentmodtrcy <- aov(log1_viable_seeds ~ Treatment, seedtrcy)
summary(treatmentmodtrcy)
treatmentmodtror <- aov(log1_viable_seeds ~ Treatment, seedtror)
summary(treatmentmodtror)
treatmentmodvero <- aov(log1_viable_seeds ~ Treatment, seedvero)
summary(treatmentmodvero)
TukeyHSD(treatmentmodvero)

######### Loops!
#Both work - top gives Chr, unique gives Factor, but c gives specified order
specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
#specieslist<- unique(seeddata$Species)
model_summaries <- list()

for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(model_summaries <- summary(
    aov(log1_viable_seeds ~ Treatment, data = filter(seeddata, Species == specieslist[i]))))
}
#Below demonstrates that this works
# test <- aov(log1_viable_seeds ~ Treatment, seedarca)
# summary(test)

#~Watering treatment

#Models~

#### Isaac's code
####run intercept-only model to determine average seed production for this species in thinned subplots only
ARCA.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Plot), data=seedarca[seedarca$C_E_or_T=="E",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))

###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(ARCA.int.mod)),exp(fixef(ARCA.int.mod)+sqrt(vcov(ARCA.int.mod))),exp(fixef(ARCA.int.mod)-sqrt(vcov(ARCA.int.mod))))
################

#Having this in model or not doesn't seem to change model output (for ARCA at least):
# control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa")

#Can't get 1|Site/Plot to run

#ARCA
#Relative to ARCA, Dry, Shade
arcamodel1 <- glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Plot), data = seedarca)
summary(arcamodel1)

arcamodel2 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), seedarca)
summary(arcamodel2)

#HYGL
hyglmodel1 <- glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Plot), data = seedhygl)
summary(hyglmodel1)

hyglmodel2 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), seedhygl)
summary(hyglmodel2)

#LARO
laromodel2 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), seedlaro)
summary(laromodel2)
#PEAI
peaimodel1 <- glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Plot), data = seedpeai)
summary(peaimodel1)

peaimodel2 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), seedpeai)
summary(peaimodel2)

#VERO
veromodel2 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), seedvero)
summary(veromodel2)

#All species
spmodel1 <- glmer.nb(No_viable_seeds ~ Treatment + Cover + Species + (1|Site), data = seeddata)
summary(spmodel1)

spmodel2 <- glmer.nb(No_viable_seeds ~ Treatment:Cover + Treatment:Species + Cover:Species + (1|Site), data = seeddata)
summary(spmodel2)

spmodel2.1 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover + Species)^2 + (1|Site), data = seeddata)
summary(spmodel2.1)

#Now trying loops for negative binomial models
specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")

for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), data = filter(seeddata, Species == specieslist[i]))))
}

#print(specieslist[i]) at the top is immensely helpful because it indicates
# which output is for which species

#Can't really figure out how to export model output cleanly
#Thought it would be useful for lining up output with species more easily
# 
# for(i in 1:length(specieslist)){
#   out <- capture.output(
#     summary(
#       glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), data = filter(seeddata, Species == specieslist[i]))))
#   write.table(out, "Testing.csv", append = TRUE, sep =",")
# }

# This is a way to store it but then when I use print it only gives the last output
#nbmodel_summaries <- list()
# for(i in 1:length(specieslist)){
#   nbmodel_summaries <- summary(
#     glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), data = filter(seeddata, Species == specieslist[i])))
# }
#print(nbmodel_summaries)

#This is a more succinct way of printing it, but again seems to be only storing last (VERO) model output
# for(i in 1:length(specieslist)){
#   print(nbmodel_summaries <- summary(
#     glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), data = filter(seeddata, Species == specieslist[i]))))
#}

#These all give the same output
# testarca1 <-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), data = seedarca)
# summary(testarca1)
# testarca2 <-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Plot), data = filter(seeddata, Species == "ARCA"))
# summary(testarca2)
# test <- filter(seeddata, Species == "ARCA")


## Loop plotting ggplots
# Testing plotting log1_viable_seeds for LARO against Cover, Treatment and C_E_or_T
#These two work, cycle through graphs to see different var plots
testvars <- c("Treatment", "Cover", "C_E_or_T")
for(i in 1:length(testvars)){
  print(ggplot(seedlaro, aes_string(x = testvars[i], y = "log1_viable_seeds"))+
    geom_boxplot()+
      geom_jitter(alpha = 0.5)+
    theme_classic())
}

#Below only works for base R plots, must use gridExtra::grid.arrange for ggplot
par(mfrow = c(1,3))
#If I want multi-column for ggplot:
testvars <- c("Treatment", "Cover", "C_E_or_T")
plotlist <- list()
for(i in 1:length(testvars)){
  plotlist[[i]] <- ggplot(seedlaro, aes_string(x = testvars[i], y = "log1_viable_seeds"))+
          geom_boxplot()+
          geom_jitter(alpha = 0.5)+
          theme_classic()
}
grid.arrange(grobs=plotlist,ncol=2)

######### Answering some of my questions!
#How do neighbours influence survival and reproduction?

#Need to remove Unknowns
seeddatanounk <- seeddata %>% filter(C_E_or_T != "UNKNOWN")

ggplot(seeddatanounk, aes(x = Neighbours, y = log1_viable_seeds))+
         geom_boxplot()+
         geom_jitter(alpha = 0.5, colour = "forestgreen")+
         theme_classic()
t.test(log1_viable_seeds ~ Neighbours, data = seeddatanounk)
#Checking sample size
seeddatanounk %>% group_by(Neighbours) %>% count()

#By Species
ggplot(seeddatanounk, aes(x = Neighbours, y = log1_viable_seeds))+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, colour = "forestgreen")+
  theme_classic()+
  facet_wrap(vars(Species))

#T-tests for each species
#Swap data between seeddatanounk and seeddataproducedseed to compare at least one seed produced and not
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    t.test(log1_viable_seeds ~ Neighbours, data = filter(seeddataproducedseed, Species == specieslist[i])))
}

## Adding in watering treatment as a layer

ggplot(seeddatanounk, aes(x = Treatment, y = log1_viable_seeds, fill = Neighbours))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(), alpha = 0.3)+
  ylab("log(Fecundity+1)")+
  xlab("Watering treatment")+
  theme_classic()+
  scale_fill_manual(values=c("gold", "deepskyblue"))+
  my_theme+
  facet_wrap(vars(Species))
#I like forest green and orange

#Model watering treatment x neigbhours by species
specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ (Treatment + Neighbours)^2 + (1|Plot), data = filter(seeddatanounk, Species == specieslist[i]))))
}

#Comparing neighbours and cover
ggplot(seeddatanounk, aes(x = Cover, y = log1_viable_seeds, fill = Neighbours))+
  geom_boxplot()+
  geom_jitter(position = position_jitterdodge(), alpha = 0.3)+
  ylab("log(Fecundity+1)")+
  xlab("Cover")+
  theme_classic()+
  scale_fill_manual(values=c("gold", "deepskyblue"))+
  my_theme
  facet_wrap(vars(Species))

#Model cover x neigbhours by species
specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ (Cover + Neighbours)^2 + (1|Plot), data = filter(seeddatanounk, Species == specieslist[i]))))
}

#Model three way interaction of Treatment, Cover and Neighbours
# Won't run - not enough data?
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ (Treatment + Cover + Neighbours)^2 + (1|Plot), data = filter(seeddatanounk, Species == specieslist[i]))))
}
#This gives the exact same output
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ Treatment + Cover + Neighbours + Treatment:Cover + Treatment:Neighbours + Cover:Neighbours + (1|Plot), data = filter(seeddatanounk, Species == specieslist[i]))))
}

##Cath's code
install.packages("sJ.Plot")
library(sJ.Plot)


drydata <- seeddatanounk %>% filter(Treatment == "Dry")
controldata <- seeddatanounk %>% filter(Treatment == "Control")
wetdata <- seeddatanounk %>% filter(Treatment == "Wet")

t.test(log1_viable_seeds ~ Neighbours, data = drydata)
#Yes, confirming that the difference is between neighbours yes/no when dry, p-value = 0.03

shadedata <- seeddatanounk %>% filter(Cover == "Shade")
sundata <- seeddatanounk %>% filter(Cover == "Sun")

t.test(log1_viable_seeds ~ Neighbours, data = sundata)

# Model of how neighbours interacts with watering treatment and shade to affect seed production
fullmodel1 <- glmer.nb(No_viable_seeds ~ (Treatment + Cover + Neighbours + Species)^2 + (1|Site), data = seeddatanounk)
summary(fullmodel1)

# By species
specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ (Treatment + Cover + Neighbours)^2 + (1|Plot), data = filter(seeddatanounk, Species == specieslist[i]))))
}

# Plotting seed production as a function of neighbour abundance
ggplot(allabundanceseeds, aes(x = log(Total_abundance+1), y = log1_viable_seeds))+
  geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
  ylab("Fecundity")+
  xlab("Total neighbour abundance")+
  geom_smooth(method="lm")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16))+
  facet_wrap(vars(Species))

ggplot(allabundanceseeds, aes(x = log(Total_abundance+1), y = log1_viable_seeds, colour = Species))+
  geom_jitter(alpha = 0.4, width = 0.1, height = 0.1)+
  theme_classic()

#Number of plants that produced at least one seed - 657
tallying <- allabundanceseeds %>% group_by(Total_abundance >= "0" & log1_viable_seeds >= "0") %>% tally()
# Number of plants that didn't produce any seeds but I haven't brought in the survival datasheet yet
filteringnas <- allabundanceseeds %>% filter(!complete.cases(log1_viable_seeds))

#Testing plots by filtering to those will intraspecific interactions
#Only 27! Mostly vero.
filterintras <- allabundanceseeds %>% filter(Intra_abundance > 1)

#Logging + 1 so that all the zeros can be properly logged.
ggplot(allabundanceseeds, aes(x = log(Intra_abundance+1), y = log1_viable_seeds))+
  geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
  ylab("Fecundity")+
  xlab("Intraspecific neighbour abundance")+
  geom_smooth(method = "lm")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16))+
  facet_wrap(vars(Species))

ggplot(allabundanceseeds, aes(x = log(Inter_abundance+1), y = log1_viable_seeds))+
  geom_jitter(alpha = 0.4, width = 0.05, height = 0.05)+
  geom_smooth(method="lm")+
  ylab("Fecundity")+
  xlab("Interspecific neighbour abundance")+
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16))+
  facet_wrap(vars(Species))

# Plotting them on the same graph
ggplot(allabundanceseeds)+
  geom_jitter(aes(x = log(Inter_abundance+1), y = log1_viable_seeds, colour = "Heterospecific"), alpha = 0.4, width = 0.05, height = 0.05)+
  geom_jitter(aes(x = log(Intra_abundance+1), y = log1_viable_seeds, colour = "Conspecific"), alpha = 0.4, width = 0.05, height = 0.05)+
  geom_smooth(aes(x = log(Inter_abundance+1), y = log1_viable_seeds, colour = "Heterospecific"), method = "lm")+
  geom_smooth(aes(x = log(Intra_abundance+1), y = log1_viable_seeds, colour = "Conspecific"), method = "lm")+
  ylab("log(Fecundity+1)")+
  xlab("log(Neighbour abundance+1)")+
  theme_classic()+
  scale_colour_manual(values = c("Heterospecific"="forestgreen", "Conspecific"="orchid"), name = NULL)+
  my_theme+
  facet_wrap(vars(Species))

#Write loops for each species to see r-squared and p-values

specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(summary(lm(log1_viable_seeds ~ log(Total_abundance+1), data = filter(allabundanceseeds, Species == specieslist[i]))))
  print(specieslist[i])
  print(summary(lm(log1_viable_seeds ~ log(Intra_abundance+1), data = filter(allabundanceseeds, Species == specieslist[i]))))
  print(specieslist[i])
  print(summary(lm(log1_viable_seeds ~ log(Inter_abundance+1), data = filter(allabundanceseeds, Species == specieslist[i]))))
}

### Model for each species with total number of neighbours
specieslist<- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ Total_abundance +(1|Plot), data = filter(allabundanceseeds, Species == specieslist[i]))))
}

#Not sure if putting them together or separately is more appropriate! 
#Esp. when data is duplicated when together
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(summary(glmer.nb(No_viable_seeds ~ Total_abundance + (1|Plot), data = filter(allabundanceseeds, Species == specieslist[i]))))
        print(specieslist[i])
        print(summary(glmer.nb(No_viable_seeds ~ Intra_abundance + (1|Plot), data = filter(allabundanceseeds, Species == specieslist[i]))))
              print(specieslist[i])
              print(summary(glmer.nb(No_viable_seeds ~ Inter_abundance + (1|Plot), data = filter(allabundanceseeds, Species == specieslist[i]))))
}

for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ Total_abundance + Intra_abundance + (1|Plot), data = filter(allabundanceseeds, Species == specieslist[i]))))
}

###################### Using Isaac's 2020 paper code ####################################
###function to predict linear curves only
lm.predict<-function(mod, newdat){ 
  pred = predict(mod, newdata = newdat, interval = "confidence", level = 0.95)
  upper<-pred[,3]
  lower<-pred[,2] 
  return(data.frame(pred=pred[,1], upper, lower))
}

###function to predict linear mixed effects models
lme.predict<-function(mod, newdat, int.type=NA, se.mult){
  pred <- as.vector(predict(mod, newdat, level = 0))
  Designmat <- model.matrix(formula(mod)[-2], newdat)
  predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat)) 
  SE <- sqrt(predvar) 
  SE2 <- sqrt(predvar+mod$sigma^2)
  if(int.type=="CI"){
    upper<-pred + (se.mult*SE)
    lower<-pred - (se.mult*SE) 
  }
  if(int.type=="PI"){
    upper<-pred + (se.mult*SE2)
    lower<-pred - (se.mult*SE2)
  }
  
  return(data.frame(pred, upper, lower))
}	

####function to plot the curve and associated condfidence interval
plot.CI.func<- function(x.for.plot, pred, upper, lower, env.colour, env.trans=NA, line.colour, line.weight){
  colour.rgb<-col2rgb(col=env.colour)  
  polygon.coords<-data.frame(rbind(cbind(x.for.plot[1], lower[1]), 
                                   cbind(x.for.plot, upper), 
                                   cbind(x.for.plot, lower)[rev(order(x.for.plot)),]))
  names(polygon.coords)<-c("x", "y")							
  polygon(polygon.coords$x, polygon.coords$y, col=rgb(red=colour.rgb["red",],blue=colour.rgb["blue",], green=colour.rgb["green",] , alpha=env.trans, maxColorValue = 255), border=NA)
  lines(x.for.plot, pred, col=line.colour, lwd=line.weight)         
} 

###function to inverse logits
inv.logit<-function(x)(exp(x)/(1+exp(x)))

##creates a sequence of values frmo min to max by 100
seq.func<-function(x)(seq(min(x, na.rm=T), max(x, na.rm=T), length.out=100))

### plot with 95% CI using lmer  or lme objects
lmer.predict<-function(mod, newdat, se.mult, binom=NULL, poisson=NULL){
  pvar1 <- diag(as.matrix(newdat) %*% tcrossprod(vcov(mod),as.matrix(newdat)))
  newdat$y<- as.matrix(newdat) %*% fixef(mod)  
  newdat <- data.frame(newdat, plo = newdat$y-(se.mult*sqrt(pvar1)), phi = newdat$y+(se.mult*sqrt(pvar1)))
  
  ## if you have used binomial errors then this will back transform logits to the probability scale
  if(binom==T) {
    newdat$y<-plogis(newdat$y); newdat$plo<-plogis(newdat$plo); newdat$phi<-plogis(newdat$phi)
  } else 
    
    ## if you have used poisson errors or have log-transformed your response, then this will back transform to the original scale (e.g. abundance)
    ##  not used in this code
    if(poisson==T) {
      newdat$y<-exp(newdat$y); newdat$plo<-exp(newdat$plo); newdat$phi<-exp(newdat$phi)
    } 
  return(with(newdat, data.frame(y, phi, plo)))
}

plot.CI.func<- function(x.for.plot, pred, upper, lower, env.colour, env.trans=NA, line.colour, line.weight, line.type){
  colour.rgb<-col2rgb(col=env.colour)  
  polygon.coords<-data.frame(rbind(cbind(x.for.plot[1], lower[1]), 
                                   cbind(x.for.plot, upper), 
                                   cbind(x.for.plot, lower)[rev(order(x.for.plot)),]))
  names(polygon.coords)<-c("x", "y")							
  polygon(polygon.coords$x, polygon.coords$y, col=rgb(red=colour.rgb["red",],blue=colour.rgb["blue",], green=colour.rgb["green",] , alpha=env.trans, maxColorValue = 255), border=NA)
  lines(x.for.plot, pred, col=line.colour, lwd=line.weight, lty=line.type)         
} 


######Q1: Do species exhibit species-specific responses to the environment?
# Cover and water availability
#Plots that have no neighbours only, looking at fecundity
########## ARCA ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
ARCA.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedarca[seedarca$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))

###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(ARCA.int.mod)),exp(fixef(ARCA.int.mod)+sqrt(vcov(ARCA.int.mod))),exp(fixef(ARCA.int.mod)-sqrt(vcov(ARCA.int.mod))))

#####run full model with quadratic terms
ARCA.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedarca[seedarca$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))

#####remove non-significant quadratic terms
ARCA.mod.q<-glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Site), data=seedarca[seedarca$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(ARCA.mod.q)

plot(seedarca[seedarca$Neighbours=="No",]$No_viable_seeds~seedarca[seedarca$Neighbours=="No",]$Treatment, pch=19, xlab="Watering treatment", ylab="Fecundity (thinned)", cex.lab=2, cex.axis=2.00,tck=-0.01, main=substitute(paste('A. calendula')), cex.main=2.5)
####Note that this pred isn't going to work because I don't have a continuous predictor variable.
####Treatment and Cover are both factors, so not plotting regressions but boxplots
#x.for.plot<-seedarca[seedarca$Neighbours=="No",]$Treatment
#y.for.plot<-seedarca[seedarca$Neighbours=="No",]$No_viable_seeds
# pred.data<-data.frame(intercept=rep(1, 100), x.for.plot=seq.func(x.for.plot), Treatment=mean(seedarca[seedarca$Neighbours=="No",]$Treatment), 0.5)
# pred<-lmer.predict(seedarca, newdat=pred.data, se.mult=1.96, binom=F, poisson=T)
# plot.CI.func(x.for.plot=seq.func(x.for.plot), pred=pred$y, upper=pred$phi, lower=pred$plo, env.colour="black", env.trans=40, line.colour="black", line.type=1, line.weight=3)
# mtext("a)", side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)

#The way that I would do the above:
seedarca %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("arca")+
  theme_classic()+
  my_theme
seedarca %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("arca")+
  theme_classic()+
  my_theme

########## HYGL ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
hygl.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedhygl[seedhygl$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(hygl.int.mod)),exp(fixef(hygl.int.mod)+sqrt(vcov(hygl.int.mod))),exp(fixef(hygl.int.mod)-sqrt(vcov(hygl.int.mod))))
#####run full model with quadratic terms
hygl.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedhygl[seedhygl$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms
hygl.mod.q<-glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Site), data=seedhygl[seedhygl$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(hygl.mod.q)
seedhygl %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("hygl")+
  theme_classic()+
  my_theme
seedhygl %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("hygl")+
  theme_classic()+
  my_theme
########## laro ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
laro.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedlaro[seedlaro$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(laro.int.mod)),exp(fixef(laro.int.mod)+sqrt(vcov(laro.int.mod))),exp(fixef(laro.int.mod)-sqrt(vcov(laro.int.mod))))
#####run full model with quadratic terms
laro.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedlaro[seedlaro$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms
laro.mod.q<-glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Site), data=seedlaro[seedlaro$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(laro.mod.q)
seedlaro %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("laro")+
  theme_classic()+
  my_theme
seedlaro %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("laro")+
  theme_classic()+
  my_theme
########## peai ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
peai.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedpeai[seedpeai$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(peai.int.mod)),exp(fixef(peai.int.mod)+sqrt(vcov(peai.int.mod))),exp(fixef(peai.int.mod)-sqrt(vcov(peai.int.mod))))
#####run full model with quadratic terms
peai.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedpeai[seedpeai$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms (none)
peai.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedpeai[seedpeai$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(peai.mod.q)
seedpeai %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("peai")+
  theme_classic()+
  my_theme
seedpeai %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("peai")+
  theme_classic()+
  my_theme
########## plde ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
plde.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedplde[seedplde$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(plde.int.mod)),exp(fixef(plde.int.mod)+sqrt(vcov(plde.int.mod))),exp(fixef(plde.int.mod)-sqrt(vcov(plde.int.mod))))
#####run full model with quadratic terms
plde.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedplde[seedplde$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms
plde.mod.q<-glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Site), data=seedplde[seedplde$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(plde.mod.q)
seedplde %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("plde")+
  theme_classic()+
  my_theme
seedplde %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("plde")+
  theme_classic()+
  my_theme
########## pole ##########
#Is there just not enough data? Why would the simpler models have trouble running but (Treatment + Cover)^2 is okay?
####run intercept-only model to determine average seed production for this species in thinned subplots only
pole.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedpole[seedpole$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(pole.int.mod)),exp(fixef(pole.int.mod)+sqrt(vcov(pole.int.mod))),exp(fixef(pole.int.mod)-sqrt(vcov(pole.int.mod))))
#####run full model with quadratic terms
pole.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedpole[seedpole$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms
pole.mod.q<-glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Site), data=seedpole[seedpole$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(pole.mod.q)
seedpole %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("pole")+
  theme_classic()+
  my_theme
seedpole %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("pole")+
  theme_classic()+
  my_theme
########## trcy ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
trcy.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedtrcy[seedtrcy$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(trcy.int.mod)),exp(fixef(trcy.int.mod)+sqrt(vcov(trcy.int.mod))),exp(fixef(trcy.int.mod)-sqrt(vcov(trcy.int.mod))))
#####run full model with quadratic terms
trcy.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedtrcy[seedtrcy$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms
trcy.mod.q<-glmer.nb(No_viable_seeds ~ Treatment + Cover + (1|Site), data=seedtrcy[seedtrcy$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(trcy.mod.q)
seedtrcy %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("trcy")+
  theme_classic()+
  my_theme
seedtrcy %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("trcy")+
  theme_classic()+
  my_theme
########## tror ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
tror.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedtror[seedtror$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(tror.int.mod)),exp(fixef(tror.int.mod)+sqrt(vcov(tror.int.mod))),exp(fixef(tror.int.mod)-sqrt(vcov(tror.int.mod))))
#####run full model with quadratic terms
tror.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedtror[seedtror$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms (none)
tror.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedtror[seedtror$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(tror.mod.q)
seedtror %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("tror")+
  theme_classic()+
  my_theme
seedtror %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("tror")+
  theme_classic()+
  my_theme
########## vero ##########
####run intercept-only model to determine average seed production for this species in thinned subplots only
vero.int.mod<-glmer.nb(No_viable_seeds ~ 1 + (1|Site), data=seedvero[seedvero$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
###back transform intercept and standard deviation around intercept estimates
c(exp(fixef(vero.int.mod)),exp(fixef(vero.int.mod)+sqrt(vcov(vero.int.mod))),exp(fixef(vero.int.mod)-sqrt(vcov(vero.int.mod))))
#####run full model with quadratic terms
vero.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedvero[seedvero$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
#####remove non-significant quadratic terms (none)
vero.mod.q<-glmer.nb(No_viable_seeds ~ (Treatment + Cover)^2 + (1|Site), data=seedvero[seedvero$Neighbours=="No",],control= glmerControl(tolPwrss=1e-5, optimizer="bobyqa"))
summary(vero.mod.q)
seedvero %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Treatment, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("vero")+
  theme_classic()+
  my_theme
seedvero %>% filter(Neighbours == "No") %>%
  ggplot(aes(x = Cover, y = log(No_viable_seeds)))+
  geom_boxplot()+
  geom_jitter(colour = "dodgerblue", alpha = 0.5)+
  ggtitle("vero")+
  theme_classic()+
  my_theme

## Importing canopy cover data and seeing how it relates to seed set and survival
canopylitterdata <- read_csv("Data/canopy_litter_2020.csv")
#Calculating canopy cover percentage
canopydata <- canopylitterdata  %>% rowwise() %>%
  mutate(meancc = mean(c(`Measurement 1`, `Measurement 2`, `Measurement 3`, `Measurement 4`, `Measurement 5`)),
        cc_percentage = meancc/24*100)
#Merging with seed data
canopyseeds <- merge(allabundanceseeds, canopydata)
#Plotting canopy cover as a continuous variable by plot
ggplot(canopyseeds, aes(x = log(cc_percentage+1), y = log(No_viable_seeds+1)))+
  geom_point(alpha = 0.4)+
  geom_smooth(method = "lm")+
  theme_classic()+
  facet_wrap(vars(Species))

for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      lmer(log(No_viable_seeds+1) ~ log(cc_percentage+1) + (1|Site/Plot),
           data = filter(canopyseeds, Species == specieslist[i]))))
}

######Q: Can community composition explain fecundity in unthinned subplots?
#Filter to unthinned data
test <- seedarca %>% filter(Neighbours == "Yes")
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  nam <- paste0("seednon", specieslist[i])
  assign(nam, filter(seeddataclean, Species == specieslist[i], Neighbours == "Yes"))
}

####check distribution of total and intra abundance for each species and where necessary transform to normal distribution 
with(seednonLARO, pairs(No_viable_seeds ~ Total_abundance + Intra_abundance + Inter_abundance, diag.panel = panel.hist))

####analyse effect of total and intra abundance on each species - seperated into species because only some species require variable transformation
#####ARCA
ARCA.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonARCA))
summary(ARCA.tot.comp)

ARCA.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonARCA))
summary(ARCA.intra.abs.comp)

#####HYGL
HYGL.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonHYGL))
summary(HYGL.tot.comp)

HYGL.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonHYGL))
summary(HYGL.intra.abs.comp)
#####LARO
LARO.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonLARO))
summary(LARO.tot.comp)

LARO.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonLARO))
summary(LARO.intra.abs.comp)
#####PEAI
PEAI.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonPEAI))
summary(PEAI.tot.comp)

PEAI.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonPEAI))
summary(PEAI.intra.abs.comp)
#####PLDE
PLDE.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonPLDE))
summary(PLDE.tot.comp)

PLDE.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonPLDE))
summary(PLDE.intra.abs.comp)
#####POLE
POLE.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonPOLE))
summary(POLE.tot.comp)

POLE.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonPOLE))
summary(POLE.intra.abs.comp)
#####TRCY
TRCY.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonTRCY))
summary(TRCY.tot.comp)

TRCY.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonTRCY))
summary(TRCY.intra.abs.comp)
#####TROR
TROR.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonTROR))
summary(TROR.tot.comp)

TROR.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonTROR))
summary(TROR.intra.abs.comp)
#####VERO
VERO.tot.comp<-(glmer.nb(No_viable_seeds ~ log(Total_abundance+1)+(1|Site), seednonVERO))
summary(VERO.tot.comp)

VERO.intra.abs.comp<-(glmer.nb(No_viable_seeds~log(Intra_abundance+1)+(1|Site), seednonVERO))
summary(VERO.intra.abs.comp)

###create a list of the dataframes which only have C data
species.c.list<-list(seednonARCA, seednonHYGL, seednonLARO, seednonPEAI, seednonPLDE, seednonPOLE, seednonTRCY, seednonTROR, seednonVERO)

###create a list of names for figure headings
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

#dev.off()
#pdf("output/FIGURE_S5.pdf", width=21, height=21)
#par(mfrow=c(3,3))
#par(mar=c(2,6,2,1))
#par(pty="s")
for(i in 1:length(species.c.list)){
  plotted.data<-as.data.frame(species.c.list[i])
  plotted.data$Intra_abundance<-log(plotted.data$Intra_abundance+1)
  plot(plotted.data$No_viable_seeds~plotted.data$Intra_abundance, pch=19, col="grey60", ylab="Fecundity (unthinned)", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "log(Number of intra-specific neighbours + 1)", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer.nb(No_viable_seeds~Intra_abundance+(1|Site), plotted.data, na.action=na.omit)
  if(summary(model)$coefficients<0.05){
    x.for.plot<-plotted.data$Intra_abundance
    y.for.plot<-plotted.data$No_viable_seeds
    pred.data<-data.frame(intercept=rep(1, 100), x.for.plot=seq.func(x.for.plot))
    pred<-lmer.predict(model, newdat=pred.data, se.mult=1.96, binom=F, poisson=T)
    plot.CI.func(x.for.plot=seq.func(x.for.plot), pred=pred$y, upper=pred$phi, lower=pred$plo, env.colour="black", env.trans=40, line.colour="black", line.type=1, line.weight=3)
  }
}

#dev.off()
#The above code doesn't plot any lines, suggesting that there are no significant relationships
#Well, can only get it to plot with coefficients[1,2] or just coefficients (Isaac had coefficients[3,4])
# Checking that here:
#Need a dataset with all species and only neighbours
seeddatan <- seeddataclean %>% filter(Neighbours == "Yes")
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
glmer.nb(No_viable_seeds ~ log(Intra_abundance+1) + (1|Site), data = filter(seeddatan, Species == specieslist[i]))))
}
#Hm, nah, PEAI and POLE are significant... and half won't plot - not enough data
print(
  summary(
    glmer.nb(No_viable_seeds ~ log(Intra_abundance+1) + (1|Site), data = seednonPEAI)))
#So maybe just plot these myself! Would like to have the data jittered with an alpha value anyway
ggplot(seeddatan, aes(x = log(Intra_abundance+1), y = No_viable_seeds))+
  geom_jitter(alpha = 0.4)+
theme_classic()+
  my_theme+
  facet_wrap(vars(Species), scales = "free")

pdf("output/FIGURE_S6.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(6,6,2,2))
par(pty="s")

for(i in 1:length(species.c.list)){
  plotted.data<-as.data.frame(species.c.list[i])
  plotted.data$Total_abundance<-log(plotted.data$Total_abundance+1)
  plot(plotted.data$No_viable_seeds~plotted.data$Total_abundance, pch=19, col="grey60", ylab="Fecundity (unthinned)", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "log(Number of neighbours + 1)", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-glmer.nb(No_viable_seeds~Total_abundance+(1|Site), plotted.data, na.action=na.omit)
  if(summary(model)$coefficients[2,4]<0.05){
    x.for.plot<-plotted.data$Total_abundance
    y.for.plot<-plotted.data$No_viable_seeds
    pred.data<-data.frame(intercept=rep(1, 100), x.for.plot=seq.func(x.for.plot))
    pred<-lmer.predict(model, newdat=pred.data, se.mult=1.96, binom=F, poisson=T)
    plot.CI.func(x.for.plot=seq.func(x.for.plot), pred=pred$y, upper=pred$phi, lower=pred$plo, env.colour="black", env.trans=40, line.colour="black", line.type=1, line.weight=3)
  }
}
dev.off()
#Arca, pole, trcy and tror coming up with lines here
#And my way below agrees!
for(i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      glmer.nb(No_viable_seeds ~ log(Total_abundance+1) + (1|Site), data = filter(seeddatan, Species == specieslist[i]))))
}
ggplot(seeddatan, aes(x = log(Total_abundance+1), y = No_viable_seeds))+
  geom_jitter(alpha = 0.3)+
  theme_classic()+
  my_theme+
  facet_wrap(vars(Species), scales = "free")

