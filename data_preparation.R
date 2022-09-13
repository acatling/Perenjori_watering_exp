#### Alexandra Catling
#### Perenjori Watering Experiment
######### Putting all the data prep and notes here
#### Updated June 2022

#### Loading packages and functions####
library(tidyverse)
library(lmerTest)
library(glmmTMB)
source("R_functions/functions.R")

### Importing a list of each of these to be filtered out
to_remove_start <- read_csv("Data/subplots_to_remove_start.csv")
to_remove_after_germ <- read_csv("Data/to_remove_after_germination.csv")

################# Germination data #############
#Importing data
germinationdataraw <- read_csv("Data/final_germination_data.csv")
#Remove blank subplots
germinationdata <- germinationdataraw %>% filter(!Rep == "Blank")
######## To be removed for various reasons:
to_remove_start$Rep <- as.factor(to_remove_start$Rep)
germinationdata <- anti_join(germinationdata, to_remove_start)
### Fixing some duplication and other individual subplot issues
# 1 A PLDE T 2 is actually a vero. at 4, 9. It's in the mortalitydata set but not listsubplots. In listsubplots it's down as VERO T4, 4, 9.
germinationdata <- within(germinationdata, Rep[Site == '1' & Plot == 'A' & Species == 'PLDE' & C_E_or_T == 'T' & Rep == '2'] <- '4')
germinationdata <- within(germinationdata, Species[Site == '1' & Plot == 'A' & Species == 'PLDE' & C_E_or_T == 'T' & Rep == '4'] <- 'VERO')
#Accidentally swapped two things to become 1 C VERO C2. One of these didn't germinate at all, position 9, 6 (E 1), so renaming this one C 4
germinationdata <- within(germinationdata, Rep[Site == '1' & Plot == 'C' & Species == 'VERO' & C_E_or_T == 'C' & Rep == '2' & Number_germinated == '0'] <- '4')

#Merge with canopy data
canopylitterdata <- read_csv("Data/canopy_litter_2020.csv")
canopydata <- canopylitterdata  %>% rowwise() %>%
  mutate(meancc = mean(c(`Measurement 1`, `Measurement 2`, `Measurement 3`, `Measurement 4`, `Measurement 5`)),
         cc_percentage = meancc/24*100)
canopydatatrim <- canopydata %>% select(Site, Plot, cc_percentage, Litter_cover_percent)
germinationdata$Site <- as.factor(germinationdata$Site)

################# Survival/mortality data #############
#Importing mortality data 
mortalitydataraw <- read_csv("Data/mortalityofgerminated.csv")
#Need to rename 8 C VERO T 1 in mortality data to E1 since note says 'stick says E1'. Already changed in seed data
mortalitydata <- within(mortalitydataraw, C_E_or_T[Site == '8' & Plot == 'C' & Species == 'VERO' & C_E_or_T == 'T'] <- 'E')
# From survey notes PEAI T2 became C2 so I am changing this above in the mortality dataset (already changed for survey and seeds)
#Check this later for the germination data, will have to change it in when I merge all vital rates because has a diff number of seeds sown for T2 and C2*
mortalitydata <- within(mortalitydata, C_E_or_T[Site == '1' & Plot == 'A' & Species == 'PEAI' & C_E_or_T == 'T' & Rep == '2'] <- 'C')
### Fixing some duplication and other individual subplot issues
#1 A LARO C 2 is two separate subplots! 1, 3 and 10, 1. Accidentally renamed E2 and T1 both to C2. Renaming the T1 (10, 1) one to C 4 instead, but remember it was imported as C 2.
mortalitydata <- within(mortalitydata, Rep[Site == '1' & Plot == 'A' & Species == 'LARO' & Row == '10' & Column == '1'] <- '4')
#There are two VERO T 2s because the focal died before thinning but then there was another germinant which also died with no seeds
#pseudo replication if I include them both from one subplot? Might just remove the second line for now and check this/make a decision later**
dup_to_remove <- mortalitydata %>% filter(Notes == "Swap. One new germ after recorded death, plot never thinned. Next germ died as well with no seeds, from checklist. See duplication below")
mortalitydata <- anti_join(mortalitydata, dup_to_remove)
#8 C LARO T 1 mislabelled row, position 7,1 is meant to be E 2 (not second T 1)
mortalitydata <- within(mortalitydata, C_E_or_T[Notes == "Check what this is labelled as. T1? No seeds. Notes say 'swapped with VERO T1'. Also from 29/9 - or was I looking at wrong plot first time? +1g"] <- 'E')
mortalitydata <- within(mortalitydata, Rep[Notes == "Check what this is labelled as. T1? No seeds. Notes say 'swapped with VERO T1'. Also from 29/9 - or was I looking at wrong plot first time? +1g"] <- '2')
# 1 A PLDE T 2 is actually a vero, recorded as VERO T4, 4, 9.
mortalitydata <- within(mortalitydata, Rep[Site == '1' & Plot == 'A' & Species == 'PLDE' & C_E_or_T == 'T' & Rep == '2'] <- '4')
mortalitydata <- within(mortalitydata, Species[Site == '1' & Plot == 'A' & Species == 'PLDE' & C_E_or_T == 'T' & Rep == '4'] <- 'VERO')

mortalitydata$Site <- as.factor(mortalitydata$Site)
to_remove_start$Site <- as.factor(to_remove_start$Site)
to_remove_after_germ$Site <- as.factor(to_remove_after_germ$Site)
to_remove_after_germ$Rep <- as.factor(to_remove_after_germ$Rep)
mortalitydata <- anti_join(mortalitydata, to_remove_start)
mortalitydata <- anti_join(mortalitydata, to_remove_after_germ)
#Assigning Neighbour values to C, E and T NAs. Removing this for now, letting survey data determine this
#Test: 2 A 9, 1 VERO T1 has Neighbours (1) --> should stay one
#mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'C' & is.na(Neighbours)] <- '1')
#mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'T' & is.na(Neighbours)] <- '0')
#mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'E' & is.na(Neighbours)] <- '0')
#CHECK THIS: Should I assign Site 1 C TROR C 1 a 0 for survival? Died
#Where no_viable_seeds and no_viable_seeds > 0, SeedsProduced = 1, otherwise 0.
#For the survival dataset, will look at SeedsProduced 1/0.
mortalitydatatrim <- mortalitydata %>% select(Site, Plot, Species, C_E_or_T, Rep,
                                              No_germinated,Survival, Neighbours)
mortalitydatatrim$Site <- as.factor(mortalitydatatrim$Site)
mortalitydatatrim$Rep <- as.factor(mortalitydatatrim$Rep)

################# Seed count data ################
#Extra seeds that were collected (e.g. broke of by accident) and flower counts where collected early or seeds missed
flowers_seeds_to_add <- read_csv("Data/flowers_seeds_to_add.csv")
#### Only adding mean number of viable and inviable seeds for four LARO cases where
#I measured flower count and made a note that I missed some (all flower count = 1)
#PEAI only has one case with 9 flowers that doesn't have enough surrounding data to estimate
#1 LARO flower = 2 viable seeds, 3 inviable seeds. Rounding.
testlaro <- flowers_seeds_to_add %>% filter(Species == "LARO", !is.na(No_viable_seeds)) 
mean(testlaro$No_viable_seeds[testlaro$Number_flowers == 1])
mean(testlaro$No_inviable_seeds[testlaro$Number_flowers == 1])
###
seeddataraw <- read_csv("Data/seed_count_2020.csv")
#Summing the number of seeds for different envelopes of same individual
seeddatarawgrouped <- seeddataraw %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(No_viable_seeds_grouped = sum(No_viable_seeds), 
         No_inviable_seeds_grouped = sum(No_inviable_seeds)) %>% filter(row_number() == 1)
seeddata <- seeddatarawgrouped
#Need Site, Treatment and Species to be factors
seeddata$Site <- factor(seeddata$Site)
#seeddata$Treatment <- factor(seeddata$Treatment)
seeddata$Species <- factor(seeddata$Species)
#Removing subplots
to_remove_start$Rep <- factor(to_remove_start$Rep)
to_remove_start$Site <- factor(to_remove_start$Site)
seeddata <- anti_join(seeddata, to_remove_start)
to_remove_after_germ$Site <- factor(to_remove_after_germ$Site)
to_remove_after_germ$Rep <- factor(to_remove_after_germ$Rep)
seeddata <- anti_join(seeddata, to_remove_after_germ)
## Adding extra missed seeds from flowers_seeds_to_add
#Replacing values in dataframe manually because I couldn't figure out merge
seeddata <- within(seeddata, No_viable_seeds_grouped[Site == '2' & Plot == 'B' & Species == 'LARO' & C_E_or_T == "E" & Rep == '2'] <- '1')
seeddata <- within(seeddata, No_inviable_seeds_grouped[Site == '2' & Plot == 'B' & Species == 'LARO' & C_E_or_T == "E" & Rep == '2'] <- '3')
seeddata <- within(seeddata, No_viable_seeds_grouped[Site == '2' & Plot == 'C' & Species == 'LARO' & C_E_or_T == "C" & Rep == '2'] <- '1')
seeddata <- within(seeddata, No_inviable_seeds_grouped[Site == '2' & Plot == 'C' & Species == 'LARO' & C_E_or_T == "C" & Rep == '2'] <- '3')
seeddata <- within(seeddata, No_viable_seeds_grouped[Site == '3' & Plot == 'C' & Species == 'LARO' & C_E_or_T == "C" & Rep == '3'] <- '1')
seeddata <- within(seeddata, No_inviable_seeds_grouped[Site == '3' & Plot == 'C' & Species == 'LARO' & C_E_or_T == "C" & Rep == '3'] <- '3')
seeddata <- within(seeddata, No_viable_seeds_grouped[Site == '8' & Plot == 'C' & Species == 'LARO' & C_E_or_T == "T" & Rep == '1'] <- '1')
seeddata <- within(seeddata, No_inviable_seeds_grouped[Site == '8' & Plot == 'C' & Species == 'LARO' & C_E_or_T == "T" & Rep == '1'] <- '3')
#Replacing values for 2C VERO C2 (+4 seeds) and 6B VERO C2 (+4 seeds)
seeddata <- within(seeddata, No_viable_seeds_grouped[Site == '2' & Plot == 'C' & Species == 'VERO' & C_E_or_T == "C" & Rep == '2'] <- '45')
seeddata <- within(seeddata, No_viable_seeds_grouped[Site == '6' & Plot == 'B' & Species == 'VERO' & C_E_or_T == "C" & Rep == '2'] <- '43')
#Plots with neighbours - C. Plots without neighbours - E and T
#DO THIS - for mortality dataset. removing this, basing neighbour presence off surveys.
#seeddata$Neighbours <- ifelse(seeddata$C_E_or_T == "C", "Yes", "No")
#seeddata$Neighbours <- factor(seeddata$Neighbours)
#Making a new column for logged viable seeds +1
seeddata$No_viable_seeds_grouped <- as.numeric(seeddata$No_viable_seeds_grouped)
seeddata$No_inviable_seeds_grouped <- as.numeric(seeddata$No_inviable_seeds_grouped)
#Check- do I need this row below? not kept in seeddatatrim
seeddata <- seeddata %>% mutate(logp1_viable_seeds = log(No_viable_seeds_grouped+1))
seeddatatrim <- seeddata %>% select(Site, Plot, Species, C_E_or_T, Rep,
                                    No_viable_seeds_grouped, No_inviable_seeds_grouped)
seeddatatrim$Rep <- as.factor(seeddatatrim$Rep)
seeddatatrim$Site <- as.factor(seeddatatrim$Site)

################# Community survey data ################
surveydataraw <- read_csv("Data/community_surveys_data.csv")
surveydataraw$Site <- factor(surveydataraw$Site)
surveydataraw$Species <- factor(surveydataraw$Species)
#Removing subplots
surveydata <- anti_join(surveydataraw, to_remove_start)
surveydata <- anti_join(surveydata, to_remove_after_germ)
#test <- surveydataraw %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(row_number() == 1)
#1065 subplots surveyed
#Fixing where subplot names don't match other datasets
#1 A LARO T1 - this is where I accidentally named two things C 2. Changed mort and germ data such that T1 (10, 1) became C 4 and 1, 3 was C 2.
#Changing names in surveys to match. Note that names and positions don't match up so going off position in notes.
surveydata <- within(surveydata, Rep[Site == '1' & Plot == 'A' & Species == 'LARO' & C_E_or_T == "C" & Rep == '2'] <- '4')
surveydata <- within(surveydata, Rep[Site == '1' & Plot == 'A' & Species == 'LARO' & C_E_or_T == "T" & Rep == '1'] <- '123')
surveydata <- within(surveydata, C_E_or_T[Site == '1' & Plot == 'A' & Species == 'LARO' & C_E_or_T == "T" & Rep == '123'] <- 'C')
surveydata <- within(surveydata, Rep[Site == '1' & Plot == 'A' & Species == 'LARO' & C_E_or_T == "C" & Rep == '123'] <- '2')
# 2 B HYGL T? 1? is indeed T 1
surveydata <- within(surveydata, Rep[Site == '2' & Plot == 'B' & Species == 'HYGL' & C_E_or_T == "T?" & Rep == '1?'] <- '1')
surveydata <- within(surveydata, C_E_or_T[Site == '2' & Plot == 'B' & Species == 'HYGL' & C_E_or_T == "T?" & Rep == '1'] <- 'T')
#2 C PLDE-POLE T 2 was renamed to POLE C4
surveydata <- within(surveydata, C_E_or_T[Site == '2' & Plot == 'C' & Species == 'PLDE-POLE' & C_E_or_T == "T" & Rep == '2'] <- 'C')
surveydata <- within(surveydata, Rep[Site == '2' & Plot == 'C' & Species == 'PLDE-POLE' & C_E_or_T == "C" & Rep == '2'] <- '4')
surveydata <- within(surveydata, Species[Site == '2' & Plot == 'C' & Species == 'PLDE-POLE' & C_E_or_T == "C" & Rep == '4'] <- 'POLE')
#Check this later - where 3 A NA NA NA data belongs 'not sure which plot this is', does not merge into vitaldata
#Not sure what TROR E 4 is, check this later, does not merge into vitaldata
# Dodder is listed under Neighbour_sp but has a Neighbour_count of NA
#Replacing Neighbour_count NAs where Neighbour_sp = Dodder with 0. e.g. so Dodder now is counted as 0
# This is because I can't tally with NAs and treating Dodder as 0 because it is in my no-comp. subplots
surveydata <- within(surveydata, Neighbour_count[Neighbour_sp == 'Dodder'] <- '0')
surveydata$Neighbour_count <- as.numeric(surveydata$Neighbour_count)
#Also want to plot the subplots with 0 neighbour_count, so need NAs to be 0s
#This filters to subplots without Neighbour_count (currently NAs, same with Neighbours)
## Check and fix 3 C PEAI C1 which has Cheilanthes and UNK 8 as NA neighbour count
## Also check and fix 7 C PEAI C2 which has UNK 1 as NA neighbour count
surveytest <- surveydata %>% filter(!complete.cases(Neighbour_count)) 
surveytest$Neighbour_count <- "0"
surveytest2 <- surveydata %>% filter(complete.cases(Neighbour_count))
## Separating abundances of intra- and inter-specific neighbours
#Can't calculate this with NAs from Neighbour_sp, which are coming from subplots with no neighbours
#Isaac's code first two lines
surveytest2$Matching <- surveytest2$Species == surveytest2$Neighbour_sp
surveytest2$Matching <- ifelse(surveytest2$Matching == TRUE, 1, 0)
surveytest2 <- surveytest2 %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  mutate(Intra_abundance = sum(Neighbour_count[Matching == "1"]),
         Inter_abundance = sum(Neighbour_count[Matching == "0"]))
#Checking numbers are right - 2157 lines in surveytest, 470 !complete.cases, 1687 complete.cases. Phew
surveytest$Neighbour_count <- as.numeric(surveytest$Neighbour_count)
surveydatafixed <- full_join(surveytest, surveytest2)
neighbourabundance <- surveydatafixed %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(Total_abundance = sum(Neighbour_count))
#Making a column for whether or not Dodder was present in subplot (1 = yes, 0 = no)
#Then filtering to one row per subplot
#Creating the below dataframe surveyabundancerows to use later
surveyabundancerows <- neighbourabundance
neighbourabundance <- neighbourabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
                  mutate(Dodder01 = case_when(any(Neighbour_sp == "Dodder") ~ "1",
                         TRUE ~ "0")) %>% filter(row_number() == 1)
#Renaming E and T Neighbours to Yes that actually have neighbours (neighbour_count >= 1)
neighbourabundance$Neighbours <- ifelse(neighbourabundance$Total_abundance > 0, "Yes", "No")
#Trimming down to just what is required to merge
surveytrim <- neighbourabundance %>% select(Site, Plot, Species, C_E_or_T, Rep, Neighbours,
                                           Total_abundance, Intra_abundance, Inter_abundance, Dodder01)
#To deal with NA problems in scaling abund values, I will change NAs to 0
#Have to remove one Species plot that is NA
surveytrim <- surveytrim %>% filter(!is.na(Species))
surveytrim <- surveytrim %>% replace(is.na(.),0)
#Creating dataframe with neighbour info only
surveytrimn <- surveytrim %>% filter(Neighbours == "Yes")
#1100 subplots surveyed, 534 of which have neighbours

#### Calculating species richness and diversity ####
#Not including dodder because we don't have an abundance count for it so can't calculate diversity inc dodder
#Species richness
diversitydata <- surveyabundancerows %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  add_tally(Neighbour_count > 0, name = "sp_richness")
#SDI = Simpsons diversity index
diversitydata <- diversitydata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(SDI_step = Neighbour_count/Total_abundance,
         SDI_step2 = SDI_step^2,
         SDI = 1-(sum(SDI_step2)))
diversitydataonerow <- diversitydata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(row_number() == 1)
SDItomerge <- diversitydataonerow %>% select(Site, Plot, Species, C_E_or_T, Rep, sp_richness, SDI)

#Shannons diversity index = shannon
#Shannon index is not being calculated correctly because of rows where neighbour_count is 0 
# including Dodder has 0 count -> NAs, and plots without neighbours
#SDI is not a problem
#remove 0 neighbour count rows, calculate shannon, merge shannon, after merge make NA shannons 0
diversitydata <- diversitydata %>% filter(Neighbour_count > 0)
diversitydata <- diversitydata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(shannon = -sum(SDI_step*(log(SDI_step))))
diversitydata <- diversitydata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(row_number() == 1)
shannontomerge <- diversitydata %>% select(Site, Plot, Species, C_E_or_T, Rep, shannon)

############### Soil/abiotic data #####
soildata <- read_csv("Data/soil_analysis_2020.csv")
soildata$Site <- as.factor(soildata$Site)
soildata <- soildata %>% mutate(log_NH4N = log(`NH4-N`), log_NO3N = log(`NO3-N`), log_Al = log(Al), log_Ca = log(Ca), log_Cu = log(Cu),
                                log_Fe = log(Fe), log_K = log(K), log_Mn = log(Mn), log_Na = log(Na), log_P = log(P), N = `NH4-N` + `NO3-N`, log_N = log(N))
soilpcadata <- soildata %>% select(Site, Plot, pH, log_NH4N, log_NO3N, log_P, log_K)
abioticpcadata <- merge(soilpcadata, canopydatatrim)
abioticpcadata <- abioticpcadata %>% unite("plotid", Site:Plot, remove = "false")
abioticpcatrim <- abioticpcadata %>% select(-c(Site, Plot))
rownames(abioticpcatrim) <- abioticpcatrim$plotid
abioticpcatrim <- as.data.frame(abioticpcatrim[,-1])
soil_pca <- princomp(abioticpcatrim, cor = TRUE)
#summary(soil_pca)
#loadings(soil_pca)
abioticpcadata$PC1<-soil_pca$scores[,1]
abioticpcadata$PC2<-soil_pca$scores[,2]
abioticpcadata$PC3<-soil_pca$scores[,3]

###### Trait data ####
##Isotope data
isotopedataraw <- read_csv("Data/leaf_isotope_data.csv")
IDkey <- read_csv("Data/leaf_isotope_ID_key.csv")
isotopedata <- cbind(isotopedataraw, IDkey)
#Transformation that John used to calculate delta from raw values (x):
# (-8-x)/(1+x/1000). This should probably be in functions
isotopedata <- isotopedata %>% mutate(delta = (-8-D13C)/(1+D13C/1000))
##Leaf area, SLA, LDMC -- leaf traits
leafdataraw <- read_csv("Data/wa_2020_leaf_traits.csv")
#Calculating LDMC (dry mass/fresh mass)
leafdata <- leafdataraw %>% filter(!is.na(Dry_mass_mg)) %>% mutate(LDMC = Dry_mass_mg/Fresh_mass_mg)
#Calculating SLA
leafdata <- leafdata %>% filter(!is.na(Leaf_area_cm2)) %>% mutate(SLA = Leaf_area_cm2/Dry_mass_mg)
traitkey <- read_csv("Data/Trait_key.csv")
coverisotope <- merge(isotopedata, traitkey)
#Count of sample size by species
countall <- isotopedata %>% group_by(Species) %>% count()
#Merging with means and averages
meanD13C <- isotopedata %>% group_by(Species) %>%
  summarise(mean_D13C = mean(delta, na.rm= TRUE),
            sd_D13C = sd(delta, na.rm=TRUE))
summaryD13C <- merge(meanD13C, countall)
#### SLA and LDMC
#Need to summarise data, raw data is three leaves per individual (not independent)
# Need a mean per individual with standard error/sd
#Then need to look at distribution of data (probably needs to be logged)
leafdatamean <- leafdata %>% group_by(Species, Collection_site, Collection_plot, Individual) %>%
  mutate(mean_SLA = mean(SLA, na.rm= TRUE),
         sd_SLA = sd(SLA, na.rm=TRUE))
leafdatamean <- leafdatamean %>% group_by(Species, Collection_site, Collection_plot, Individual) %>%
  mutate(mean_LDMC = mean(LDMC, na.rm= TRUE),
         sd_LDMC = sd(LDMC, na.rm=TRUE))
#Reducing to one row
leafsimple <- leafdatamean %>% group_by(Species, Collection_site, Collection_plot, 
                                        Individual, mean_SLA, sd_SLA, mean_LDMC, sd_LDMC) %>%
  filter(row_number() == 1)
#Calculating one value for each species, based on average of individual averages
speciestraits <- leafsimple %>% group_by(Species) %>% summarise(SLA = mean(mean_SLA),
                                                                sd_SLA = sd(mean_SLA),
                                                                LDMC = mean(mean_LDMC),
                                                                sd_LDMC = mean(mean_LDMC))
#Merging with isotope species summaries
alltraits <- merge(speciestraits, summaryD13C)
#This gives a summary value for sun and shade traits per species
leafsimple2 <- merge(leafsimple, traitkey)
covertraits <- leafsimple2 %>% group_by(Species, Cover) %>% summarise(SLA = mean(mean_SLA),
                                                                      sd_SLA = sd(mean_SLA),
                                                                      LDMC = mean(mean_LDMC),
                                                                      sd_LDMC = mean(mean_LDMC))
### This gives info per individual on sun vs shade
leafsimplecover <- leafsimple2 %>% group_by(Species, Collection_site, Collection_plot, Individual) %>% 
  summarise(SLA = mean(mean_SLA),
            sd_SLA = sd(mean_SLA),
            LDMC = mean(mean_LDMC),
            sd_LDMC = mean(mean_LDMC))
#Merge back in Cover info
leafsimplecover <- merge(leafsimplecover, traitkey)
#####Vero - can look within plots 
verotraits <- leafsimple %>% filter(Species == "VERO")
veroplot <- verotraits %>% group_by(Collection_site, Collection_plot) %>% 
                          summarise(SLA = mean(mean_SLA),
                          sd_SLA = sd(mean_SLA),
                          LDMC = mean(mean_LDMC),
                          sd_LDMC = mean(mean_LDMC),
                          area = mean(Leaf_area_cm2),
                          sd_area = sd(Leaf_area_cm2))
#Merge back with info on treatments
verotreatments <- merge(veroplot, traitkey)
#Plotting traits against canopy cover as continuous variable
verositeplot <- verotreatments %>% mutate (Site = Collection_site, Plot = Collection_plot)
verocover <- merge(verositeplot, canopydata)
#And isotope data
veroisotope <- isotopedata %>% filter(Species == "VERO") %>% mutate(Site = Collection_site,
                                                                    Plot = Collection_plot)
verocoverisotope <- merge(veroisotope, canopydata)

#### All data together ####
#### Trying to merge germ, surv and mort data onto complete list of subplots 
#1614 rows with germination data
#1139 rows with survival data
#653 rows with seed production data
listsubplots <- germinationdata %>% select(Site, Plot, Species, C_E_or_T, Rep, Number_seeds_sown, February_germination, Number_germinated, Notes_germination)
#1614 subplots with seeds sown, after excluding blanks and things that were to be removed at start
### Merge in survival data, mortalitydatatrim which has 1139 rows
mortalitydatatrim$Site <- as.factor(mortalitydatatrim$Site)
listsubplots$Rep <- as.factor(listsubplots$Rep)
survtomerge <- mortalitydatatrim %>% select(Site, Plot, Species, C_E_or_T, Rep, Survival)
germsurv <- left_join(listsubplots, survtomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"))

## Looking into why there are discrepancies more closely:
#So basically I want to use left join to keep all the rows from listsubplots,
# but then I also need to check which rows were missed from the dataframe I am adding in, and fix that
#library(arsenal)
#test <- comparedf(listsubplots, survtomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"), factor.as.char = TRUE)
#summary(test)
# The 'observations not shared' table shows rows in data frame x (listsubplots) that aren't in y (mortalitydatatrim) and vice versa. y is at the bottom

### Adding in seed production data, seeddatatrim which has 653 rows
germsurv$Rep <- as.factor(germsurv$Rep)
vitaldata <- left_join(germsurv, seeddatatrim)
# Note that it is not merging in the six unknown plots in seeddatatrim, fine for now
### Adding in info on abiotic environment and traits
vitaldata <- left_join(vitaldata, abioticpcadata)
traitdata <- alltraits %>% select(Species, SLA, LDMC, mean_D13C)
vitaldata <- left_join(vitaldata, traitdata)

### Adding in survey data
#1057 rows with neighbourhood surveys done
vitaldata <- left_join(vitaldata, surveytrim, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"))
# Subplots in surveytrim that aren't in vitaldata: 1 A LARO T1 (fixed, check), 2 B HYGL T? 1? (fixed), 2 C PLDE-POLE T 2 (fixed), 3 A NA NA NA (check later), 4 A TROR E 4 (check later)
### Adding in diversity data from surveys
vitaldata <- left_join(vitaldata, SDItomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"))
vitaldata <- left_join(vitaldata, shannontomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"))
#Make SDI and shannon NA rows 0.
#NAs are where there are no neighbours (or just dodder for SDI)
vitaldata <- within(vitaldata, SDI[is.na(SDI)] <-  '0')
vitaldata <- within(vitaldata, shannon[is.na(shannon)] <-  '0')
#Dodder is currently NA for plots that weren't surveyed, making it 0 for these (I recorded where it was present, so for everywhere I didn't record it it's a 0)
vitaldata <- within(vitaldata, Dodder01[is.na(Dodder01)] <-  '0')

## Want to calculate the relevant things, e.g. Total_abundance and ProducedSeeds
# But don't want to override NAs for subplots where it's not relevant - i.e. those that never germinated
# So going to calculate them separately and then merge it back in
### Germination relevant calculations ##
#For subplots that germinated, replace all the NAs for seed counts with 0s - below code replaces everything in dataframe
# that is an NA with 0 - that's okay here, where only No_viable and inviable seeds have NAs
germinatedsubplots <- vitaldata %>% filter(Number_germinated > 0 | February_germination > 0)
#1171 subplots germinated
germcalcs <- germinatedsubplots %>% select(Site, Plot, Species, C_E_or_T, Rep, Number_seeds_sown, February_germination, Number_germinated, No_viable_seeds_grouped, No_inviable_seeds_grouped)
germcalcs <- germcalcs %>% replace(is.na(.),0)
#Make a new column for survival to produce seeds (1 or 0)
germcalcs <- germcalcs %>% mutate(surv_to_produce_seeds = case_when(No_viable_seeds_grouped > "0" | No_inviable_seeds_grouped > "0" ~ "1",
                                                                    No_viable_seeds_grouped == "0" & No_inviable_seeds_grouped == "0" ~ "0"))
#Merge back into big dataframe
germcalcstomerge <- germcalcs %>% select(Site, Plot, Species, C_E_or_T, Rep, surv_to_produce_seeds, No_viable_seeds_grouped)
vitaltomerge <- vitaldata %>% select(-c(Notes_germination, Survival, No_viable_seeds_grouped, No_inviable_seeds_grouped))
vitaldata <- left_join(vitaltomerge, germcalcstomerge, by =c("Site", "Plot", "Species", "C_E_or_T", "Rep"))

#604 produced seeds (at least one inviable or viable), and 551 didn't produce any
#Want to use cbind for number of successes and failures, using binomial glm.
#Create a row for total number germinated, total number didn't
vitaldata <- vitaldata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(total_germ = sum(Number_germinated+February_germination),
         total_no_germ = Number_seeds_sown-total_germ,
         percent_germ = total_germ/Number_seeds_sown)
###Scaling the number of seeds produced across species to a percentage of maximum produced
#Need these values to be integers for negative binomial models
#note that below isn't working
#vitaldata <- vitaldata %>% group_by(Species) %>% mutate(seeds_percent = round(No_viable_seeds_grouped/max(No_viable_seeds_grouped)*100))

#Merging with information on treatments (watering and cover (shade/sun) per plot)
treatments <- read_csv("Data/treatments_meta.csv")
treatments$Site <- as.factor(treatments$Site)
vitaldata <- left_join(vitaldata, treatments)

#Assigning 1 to subplots that had neighbours, based on surveys where total_abundance > 0
vitaldata$Neighbours01 <- ifelse(vitaldata$Total_abundance > 0, "1", "0")
#test <- vitaldata %>% filter(!(is.na(surv_to_produce_seeds)) & is.na(Neighbours01))
#Why are there 86 plots where focal germinated but neighbour count is not recorded?
#Mix of plants that survived and didn't, mix of Cs, Es and Ts
# Don't have survey data for them, so assuming no neighbours
vitaldata <- within(vitaldata, Total_abundance[!(is.na(surv_to_produce_seeds)) & is.na(Neighbours01)] <- '0')
vitaldata <- within(vitaldata, Intra_abundance[!(is.na(surv_to_produce_seeds)) & is.na(Neighbours01)] <- '0')
vitaldata <- within(vitaldata, Inter_abundance[!(is.na(surv_to_produce_seeds)) & is.na(Neighbours01)] <- '0')
#Also updating Neighbours01 and dropping Neighbours column to avoid confusion
vitaldata <- within(vitaldata, Neighbours01[!(is.na(surv_to_produce_seeds)) & is.na(Neighbours01)] <- '0')
vitaldata <- vitaldata %>% select(-(Neighbours))

#Need total abund to be numeric
vitaldata$Total_abundance <- as.numeric(vitaldata$Total_abundance)
vitaldata$Intra_abundance <- as.numeric(vitaldata$Intra_abundance)
vitaldata$Inter_abundance <- as.numeric(vitaldata$Inter_abundance)

### Checking distributions and standardising data ##
#Distributions of traits and PC1, have to use datanotscaled to log
#SLA better logged?
#LDMC still strange logged
#WUE/D13C better logged but still left-skewed
#PC1 - right-skewed, logging doesn't help
#with(vitaldata, pairs(~SLA + LDMC + mean_D13C + PC1, diag.panel = panel.hist))
#with(vitaldata, pairs(~log(SLA+2) + LDMC + log(mean_D13C+2) + PC1, diag.panel = panel.hist))
#Logging and standardising variables
vitaldata$log_SLA <-log(vitaldata$SLA)
vitaldata$log_D13C <-log(vitaldata$mean_D13C)
vitaldata$logp1_totalabund <- log(vitaldata$Total_abundance+1)
vitaldata$logp1_interabund <- log(vitaldata$Inter_abundance+1)
vitaldata$logp1_intraabund <-log(vitaldata$Intra_abundance+1)

#Standardising continuous explanatory variables to a mean of 0 and SD of 1
vitaldata$std_cc <- scale(vitaldata$cc_percentage, center = TRUE, scale = TRUE)[,1]
vitaldata$std_PC1 <- scale(vitaldata$PC1, center = TRUE, scale = TRUE)[,1]
vitaldata$std_PC2 <- scale(vitaldata$PC2, center = TRUE, scale = TRUE)[,1]
vitaldata$std_PC3 <- scale(vitaldata$PC3, center = TRUE, scale = TRUE)[,1]
vitaldata$std_log_SLA <- scale(vitaldata$log_SLA, center = TRUE, scale = TRUE)[,1]
vitaldata$std_LDMC <- scale(vitaldata$LDMC, center = TRUE, scale = TRUE)[,1]
vitaldata$std_log_D13C <- scale(vitaldata$log_D13C, center = TRUE, scale = TRUE)[,1]

#Standardising continuous explanatory variables to a mean of 0 and SD of 1
vitaldata$std_logp1_totalabund <- scale(vitaldata$logp1_totalabund, center = TRUE, scale = TRUE)[,1]
vitaldata$std_logp1_interabund <- scale(vitaldata$logp1_interabund, center = TRUE, scale = TRUE)[,1]
vitaldata$std_logp1_intraabund <- scale(vitaldata$logp1_intraabund, center = TRUE, scale = TRUE)[,1]

#Need surv_to_produce_seeds to be numeric, not a factor to plot curves
vitaldata$surv_to_produce_seeds <- as.numeric(vitaldata$surv_to_produce_seeds)

#Renaming Control watering treatment to Ambient
vitaldata <- within(vitaldata, Treatment[Treatment == 'Control'] <- 'Ambient')

## Need SDI and shannon to be numeric
vitaldata$SDI <- as.numeric(vitaldata$SDI)
vitaldata$shannon <- as.numeric(vitaldata$shannon)

#Adding unique row ID
vitaldata <- rownames_to_column(vitaldata, var = "rowID") %>% as_tibble()

#Splitting data by species
arcadata <- vitaldata %>% filter(Species == "ARCA")
hygldata <- vitaldata %>% filter(Species == "HYGL")
larodata <- vitaldata %>% filter(Species == "LARO")
peaidata <- vitaldata %>% filter(Species == "PEAI")
pldedata <- vitaldata %>% filter(Species == "PLDE")
poledata <- vitaldata %>% filter(Species == "POLE")
trcydata <- vitaldata %>% filter(Species == "TRCY")
trordata <- vitaldata %>% filter(Species == "TROR")
verodata <- vitaldata %>% filter(Species == "VERO")

#Seed production model data, which is filtered to things that produced
#at least one seed (viable or inviable)
seedmodeldata <- vitaldata %>% filter(surv_to_produce_seeds == "1")
seedarca <- seedmodeldata %>% filter(Species == "ARCA")
seedhygl <- seedmodeldata %>% filter(Species == "HYGL")
seedlaro <- seedmodeldata %>% filter(Species == "LARO")
seedpeai <- seedmodeldata %>% filter(Species == "PEAI")
seedplde <- seedmodeldata %>% filter(Species == "PLDE")
seedpole <- seedmodeldata %>% filter(Species == "POLE")
seedtrcy <- seedmodeldata %>% filter(Species == "TRCY")
seedtror <- seedmodeldata %>% filter(Species == "TROR")
seedvero <- seedmodeldata %>% filter(Species == "VERO")

### Creating dataset with neighbours only
#Check why some plots weren't surveyed - I thought I surveyed everything that germinated
#test <- anti_join(mortalitydatatrim, surveytrim) #107 rows
datanonly <- vitaldata %>% filter(Total_abundance > 0)
#Reordering watering treatments to  Dry, Ambient, Wet for plotting
datanonly$Treatment <- factor(datanonly$Treatment, level = c("Dry", "Ambient", "Wet"))
vitaldata$Treatment <- factor(vitaldata$Treatment, level = c("Dry", "Ambient", "Wet"))
datanonly$SDI <- as.numeric(datanonly$SDI)
datanonly$shannon <- as.numeric(datanonly$shannon)

#Creating species lists for loops
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
species.list.s<-list(arcadata, hygldata, larodata, peaidata, pldedata, poledata, trcydata, trordata, verodata)
species.list.f<-list(seedarca, seedhygl, seedlaro, seedpeai, seedplde, seedpole, seedtrcy, seedtror, seedvero)
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

### Dropping POLE
specieslist.nop <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "TRCY", "TROR", "VERO")
species.list.s.nop <-list(arcadata, hygldata, larodata, peaidata, pldedata, trcydata, trordata, verodata)
species.list.f.nop <-list(seedarca, seedhygl, seedlaro, seedpeai, seedplde, seedtrcy, seedtror, seedvero)
species.name.list.nop <-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

#### Calculating population growth rate ####
#popdata is the dataframe to use in models

### Extracting seed survival values from Isaac Tower's Oecologia 2022 paper
#Isaac's code with slight modifications, including calculating average:
#####bring in seed fill information
seed_fill<-read_csv("Data_and_code_from_others/seed_fill_towers_oecologia_2022.csv") %>%
  select(Species = 'species', plot, num_seeds = 'no.seeds', num_seeds_visible = 'no.seeds.visible', num_seeds_unfilled = 'no.seeds.unfilled') %>%
  mutate(num_filled = num_seeds_visible - num_seeds_unfilled, percentage_filled = num_filled/num_seeds_visible)
###find the species-level median fill rate just to check (i.e. species-level fill) and the upper and lower 95% quantile as a measure of range, we will use plot-specific fill rates to adjust seed numbers
species_level_seed_fill<-seed_fill%>%
  drop_na(percentage_filled)%>%
  group_by(Species)%>%
  summarise(median.species.fill=median(percentage_filled),l95.species.fill=quantile(percentage_filled, 0.05),h95median.species.fill=quantile(percentage_filled, 0.95),
            mean.species.fill = mean((percentage_filled))) %>%
  mutate_if(is.numeric, round, digits =2)
#Selecting my species only and calculating an average from those values (for remaining 2 species, PEAI and POLE)
#Creating a table with species names
speciestable <- c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides", "Plantago debilis","Podolepis lessonii", "Trachymene cyanopetala","Trachymene ornata","Velleia rosea")
speciestable <- data.frame(speciestable)
speciestable <- speciestable %>% select(Species = 'speciestable')
#Merging in seed fill data
speciestable <- left_join(speciestable, species_level_seed_fill) %>% select(Species, mean.species.fill)
#Calculating an average for PEAI and POLE
speciestable %>% summarise(mean(mean.species.fill, na.rm = TRUE))
#Mean is 0.96, assigning this value to PEAI and POLE
speciestable <- within(speciestable, mean.species.fill[Species == "Pentameris airoides"] <- 0.96)
speciestable <- within(speciestable, mean.species.fill[Species == "Podolepis lessonii"] <- 0.96)
#Replacing names with my abbreviations
speciestable <- within(speciestable, Species[Species == "Arctotheca calendula"] <- 'ARCA')
speciestable <- within(speciestable, Species[Species == "Hyalosperma glutinosum"] <- 'HYGL')
speciestable <- within(speciestable, Species[Species == "Lawrencella rosea"] <- 'LARO')
speciestable <- within(speciestable, Species[Species == "Pentameris airoides"] <- 'PEAI')
speciestable <- within(speciestable, Species[Species == "Plantago debilis"] <- 'PLDE')
speciestable <- within(speciestable, Species[Species == "Podolepis lessonii"] <- 'POLE')
speciestable <- within(speciestable, Species[Species == "Trachymene cyanopetala"] <- 'TRCY')
speciestable <- within(speciestable, Species[Species == "Trachymene ornata"] <- 'TROR')
speciestable <- within(speciestable, Species[Species == "Velleia rosea"] <- 'VERO')
speciestable <- speciestable %>% select(Species, seed_survival = 'mean.species.fill')

# Per capita growth rate of a given population  i = seed survival*(1-germination)+number of viable seeds produced per germinant*germination
#Calculate this at the plot level

## Survival -- TRCY not converging, fixed
#survival, glmer, binom
#plogis instead of exp

#Create dataframe
Species <- rep(c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO"), each = 3, times = 8)
Site <- rep(c("1", "2", "3", "4", "5", "6", "7", "8"), each = 9, times = 3)
Plot <- rep(c("A", "B", "C"), times = 72)
survpopdataframe <- cbind(Species, Site, Plot)
survpopdataframe <- data.frame(survpopdataframe)
survpopdataframe <- survpopdataframe %>% unite("idforjoining", Plot:Site, sep = ":", remove = "false")

#ARCA 
arcapopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, arcadata)
arcasurvplot_means <- coef(arcapopsurvmodel)$Plot
arcasurvpop <- cbind(idforjoining = rownames(arcasurvplot_means), arcasurvplot_means)
arcasurvpop$means_no_neighbours <- plogis(arcasurvpop[,2])
arcasurvpop$means_neighbours <- plogis(arcasurvpop[,2] + arcasurvpop[,3])
arcasurvpop$Species <- 'ARCA'
arcasurvpoptomerge <- arcasurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframe, arcasurvpoptomerge)

#HYGL
hyglpopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, hygldata)
hyglsurvplot_means <- coef(hyglpopsurvmodel)$Plot
hyglsurvpop <- cbind(idforjoining = rownames(hyglsurvplot_means), hyglsurvplot_means)
hyglsurvpop$means_no_neighbours <- plogis(hyglsurvpop[,2])
hyglsurvpop$means_neighbours <- plogis(hyglsurvpop[,2] + hyglsurvpop[,3])
hyglsurvpop$Species <- 'HYGL'
hyglsurvpoptomerge <- hyglsurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, hyglsurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#LARO
laropopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, larodata)
larosurvplot_means <- coef(laropopsurvmodel)$Plot
larosurvpop <- cbind(idforjoining = rownames(larosurvplot_means), larosurvplot_means)
larosurvpop$means_no_neighbours <- plogis(larosurvpop[,2])
larosurvpop$means_neighbours <- plogis(larosurvpop[,2] + larosurvpop[,3])
larosurvpop$Species <- 'LARO'
larosurvpoptomerge <- larosurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, larosurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#PEAI
peaipopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, peaidata)
peaisurvplot_means <- coef(peaipopsurvmodel)$Plot
peaisurvpop <- cbind(idforjoining = rownames(peaisurvplot_means), peaisurvplot_means)
peaisurvpop$means_no_neighbours <- plogis(peaisurvpop[,2])
peaisurvpop$means_neighbours <- plogis(peaisurvpop[,2] + peaisurvpop[,3])
peaisurvpop$Species <- 'PEAI'
peaisurvpoptomerge <- peaisurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, peaisurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#PLDE
pldepopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, pldedata)
pldesurvplot_means <- coef(pldepopsurvmodel)$Plot
pldesurvpop <- cbind(idforjoining = rownames(pldesurvplot_means), pldesurvplot_means)
pldesurvpop$means_no_neighbours <- plogis(pldesurvpop[,2])
pldesurvpop$means_neighbours <- plogis(pldesurvpop[,2] + pldesurvpop[,3])
pldesurvpop$Species <- 'PLDE'
pldesurvpoptomerge <- pldesurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, pldesurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#POLE
polepopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, poledata)
polesurvplot_means <- coef(polepopsurvmodel)$Plot
polesurvpop <- cbind(idforjoining = rownames(polesurvplot_means), polesurvplot_means)
polesurvpop$means_no_neighbours <- plogis(polesurvpop[,2])
polesurvpop$means_neighbours <- plogis(polesurvpop[,2] + polesurvpop[,3])
polesurvpop$Species <- 'POLE'
polesurvpoptomerge <- polesurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, polesurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#TRCY
trcypopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, trcydata)
trcysurvplot_means <- coef(trcypopsurvmodel)$Plot
trcysurvpop <- cbind(idforjoining = rownames(trcysurvplot_means), trcysurvplot_means)
trcysurvpop$means_no_neighbours <- plogis(trcysurvpop[,2])
trcysurvpop$means_neighbours <- plogis(trcysurvpop[,2] + trcysurvpop[,3])
trcysurvpop$Species <- 'TRCY'
trcysurvpoptomerge <- trcysurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, trcysurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#TROR
trorpopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, trordata)
trorsurvplot_means <- coef(trorpopsurvmodel)$Plot
trorsurvpop <- cbind(idforjoining = rownames(trorsurvplot_means), trorsurvplot_means)
trorsurvpop$means_no_neighbours <- plogis(trorsurvpop[,2])
trorsurvpop$means_neighbours <- plogis(trorsurvpop[,2] + trorsurvpop[,3])
trorsurvpop$Species <- 'TROR'
trorsurvpoptomerge <- trorsurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, trorsurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#VERO
veropopsurvmodel <- glmer(surv_to_produce_seeds ~ Neighbours01 + (Neighbours01|Site/Plot),family = binomial, verodata)
verosurvplot_means <- coef(veropopsurvmodel)$Plot
verosurvpop <- cbind(idforjoining = rownames(verosurvplot_means), verosurvplot_means)
verosurvpop$means_no_neighbours <- plogis(verosurvpop[,2])
verosurvpop$means_neighbours <- plogis(verosurvpop[,2] + verosurvpop[,3])
verosurvpop$Species <- 'VERO'
verosurvpoptomerge <- verosurvpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
survpopdataframemerged <- left_join(survpopdataframemerged, verosurvpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

### Seed production - trcy not converging
#Create dataframe
Species <- rep(c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO"), each = 3, times = 8)
Site <- rep(c("1", "2", "3", "4", "5", "6", "7", "8"), each = 9, times = 3)
Plot <- rep(c("A", "B", "C"), times = 72)
merge <- cbind(Species, Site, Plot)
seedpopdataframe <- data.frame(merge)
seedpopdataframe <- seedpopdataframe %>% unite("idforjoining", Plot:Site, sep = ":", remove = "false")

#ARCA
arcapopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|Site/Plot), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedarca)
arcaseedplot_means<-coef(arcapopseedmodel)$Plot
arcaseedpop <- cbind(idforjoining = rownames(arcaseedplot_means), arcaseedplot_means)
arcaseedpop$means_no_neighbours<-exp(arcaseedpop[,2])
arcaseedpop$means_neighbours<-exp(arcaseedpop[,2] + arcaseedpop[,3])
arcaseedpop$Species <- 'ARCA'
arcaseedpoptomerge <- arcaseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframe, arcaseedpoptomerge)

#HYGL
hyglpopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|Site/Plot), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedhygl)
hyglseedplot_means<-coef(hyglpopseedmodel)$Plot
hyglseedpop <- cbind(idforjoining = rownames(hyglseedplot_means), hyglseedplot_means)
hyglseedpop$means_no_neighbours<-exp(hyglseedpop[,2])
hyglseedpop$means_neighbours<-exp(hyglseedpop[,2] + hyglseedpop[,3])
hyglseedpop$Species <- 'HYGL'
hyglseedpoptomerge <- hyglseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, hyglseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#LARO - won't converge, separated Site/Plot
laropopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (1|Site) + (Neighbours01|plotid), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedlaro)
laroseedplot_means<-coef(laropopseedmodel)$plotid
rownames(laroseedplot_means) <- c("A:1", "B:1", "C:1", "A:2", "B:2", "C:2", "A:3", "B:3", "C:3", "A:4", "B:4", "C:4", "A:5", "B:5", "A:6", "B:6", "C:6", "A:7", "B:7", "C:7", "A:8", "B:8", "C:8")
laroseedpop <- cbind(idforjoining = rownames(laroseedplot_means), laroseedplot_means)
laroseedpop$means_no_neighbours<-exp(laroseedpop[,2])
laroseedpop$means_neighbours<-exp(laroseedpop[,2] + laroseedpop[,3])
laroseedpop$Species <- 'LARO'
laroseedpoptomerge <- laroseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, laroseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#PEAI - won't converge, even with separated site and plot, even with site removed
#will run with lmer.........
peaipopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|plotid), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedpeai)
#laroseedplot_means<-coef(peaipopseedmodel)$plotid
### here 03/08/2022*
# ggplot(seedpeai, aes(x = Neighbours01, y = log(No_viable_seeds_grouped+1)))+
#   geom_boxplot()+
#   geom_jitter()+
#   theme_classic()
peaipopseedmodel2 <- lmer(log(No_viable_seeds_grouped+1) ~ Neighbours01 + (Neighbours01|Site/Plot), seedpeai)
peaiseedplot_means<-coef(peaipopseedmodel2)$Plot
peaiseedpop <- cbind(idforjoining = rownames(peaiseedplot_means), peaiseedplot_means)
peaiseedpop$means_no_neighbours<-exp(peaiseedpop[,2])
peaiseedpop$means_neighbours<-exp(peaiseedpop[,2] + peaiseedpop[,3])
peaiseedpop$Species <- 'PEAI'
peaiseedpoptomerge <- peaiseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, peaiseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#PLDE
pldepopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|Site/Plot), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedplde)
pldeseedplot_means<-coef(pldepopseedmodel)$Plot
pldeseedpop <- cbind(idforjoining = rownames(pldeseedplot_means), pldeseedplot_means)
pldeseedpop$means_no_neighbours<-exp(pldeseedpop[,2])
pldeseedpop$means_neighbours<-exp(pldeseedpop[,2] + pldeseedpop[,3])
pldeseedpop$Species <- 'PLDE'
pldeseedpoptomerge <- pldeseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, pldeseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#POLE
polepopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|Site/Plot), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedpole)
poleseedplot_means<-coef(polepopseedmodel)$Plot
poleseedpop <- cbind(idforjoining = rownames(poleseedplot_means), poleseedplot_means)
poleseedpop$means_no_neighbours<-exp(poleseedpop[,2])
poleseedpop$means_neighbours<-exp(poleseedpop[,2] + poleseedpop[,3])
poleseedpop$Species <- 'POLE'
poleseedpoptomerge <- poleseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, poleseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#TRCY - won't converge, even with separated site and plot, will converge if I omit Site
trcypopseedmodel3 <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|plotid), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedtrcy)
trcyseedplot_means<-coef(trcypopseedmodel3)$plotid
#Rename rownames to fit this different model type, problem with this method is I need to check exact names/numbers of plots
rownames(trcyseedplot_means) <- c("A:1", "B:1", "C:1", "A:2", "B:2", "C:2", "A:3", "B:3", "C:3", "A:4", "B:4", "C:4", "A:5", "B:5", "C:5", "A:6", "B:6", "C:6", "A:7", "B:7", "C:7", "A:8", "B:8", "C:8")
trcyseedpop <- cbind(idforjoining = rownames(trcyseedplot_means), trcyseedplot_means)
trcyseedpop$means_no_neighbours<-exp(trcyseedpop[,2])
trcyseedpop$means_neighbours<-exp(trcyseedpop[,2] + trcyseedpop[,3])
trcyseedpop$Species <- 'TRCY'
trcyseedpoptomerge <- trcyseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, trcyseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#TROR
trorpopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|Site/Plot), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedtror)
trorseedplot_means<-coef(trorpopseedmodel)$Plot
trorseedpop <- cbind(idforjoining = rownames(trorseedplot_means), trorseedplot_means)
trorseedpop$means_no_neighbours<-exp(trorseedpop[,2])
trorseedpop$means_neighbours<-exp(trorseedpop[,2] + trorseedpop[,3])
trorseedpop$Species <- 'TROR'
trorseedpoptomerge <- trorseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, trorseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#VERO, have to simplify to just Neighbours01|plotid to run
veropopseedmodel <- glmer.nb(No_viable_seeds_grouped ~ Neighbours01 + (Neighbours01|plotid), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), seedvero)
veroseedplot_means<-coef(veropopseedmodel)$plotid
rownames(veroseedplot_means) <- c("A:1", "B:1", "C:1", "A:2", "B:2", "C:2", "A:3", "B:3", "C:3", "B:4", "C:4", "A:5", "B:5", "A:6", "B:6", "C:6", "A:7", "B:7", "C:7", "A:8", "C:8")
veroseedpop <- cbind(idforjoining = rownames(veroseedplot_means), veroseedplot_means)
veroseedpop$means_no_neighbours<-exp(veroseedpop[,2])
veroseedpop$means_neighbours<-exp(veroseedpop[,2] + veroseedpop[,3])
veroseedpop$Species <- 'VERO'
veroseedpoptomerge <- veroseedpop %>% select(Species, idforjoining, means_no_neighbours, means_neighbours)
seedpopdataframemerged <- left_join(seedpopdataframemerged, veroseedpoptomerge, by= c("Species", "idforjoining")) %>%
  mutate(means_no_neighbours = coalesce(means_no_neighbours.x, means_no_neighbours.y),
         means_neighbours = coalesce(means_neighbours.x, means_neighbours.y)) %>%
  select(Species, idforjoining, Site, Plot, means_no_neighbours, means_neighbours)

#Renaming column headings to differentiate survival and seed production in big dataframe
seedpop <- seedpopdataframemerged %>% rename(seed_means_no_neighbours = means_no_neighbours,
                                             seed_means_neighbours = means_neighbours)
survpop <- survpopdataframemerged %>% rename(surv_means_no_neighbours = means_no_neighbours,
                                             surv_means_neighbours = means_neighbours)
popgrowthratedata <- merge(seedpop, survpop)
#Adding in germination rates per plot
plotgermrates <- vitaldata %>% group_by(Species, Site, Plot) %>% summarise(plot_germ = mean(percent_germ, na.rm = TRUE))
popgrowthratedata <- full_join(plotgermrates, popgrowthratedata)

#### Comparing survival probability and proportion ####
# Calculating survival values as proportion survival at plot level
# to see if this concurs with calculated probability of survivals
plotsurvtest <- vitaldata %>% group_by(Species, Site, Plot, Neighbours01, surv_to_produce_seeds) %>% count()
#Why so many Neighbours01 NAs?
#These are for the ones that didn't germinate, which is why they are survival NAs too

## Want to create a dataframe from scratch (9 species, 3 plots, 8 sites = 216)
#9 species, 8 sites, 3 plots, 2 levels neighbours, 2 levels survival = 864
test1000 <- survpop %>% select(Species, Site, Plot)
test1000a <- test1000 %>% mutate(Neighbours01=1, surv_to_produce_seeds=1)
test1000b <- test1000 %>% mutate(Neighbours01=0, surv_to_produce_seeds=0)
test1000c <-test1000 %>% mutate(Neighbours01=0, surv_to_produce_seeds=1)
test1000d <-test1000 %>% mutate(Neighbours01=1, surv_to_produce_seeds=0)
#Also need where Neighbours01=0, surv_to_produce_seeds=1 and
#where Neighbours01=1, surv_to_produce_seeds=0
test100 <- rbind(test1000a,test1000b, test1000c, test1000d)
#Left_join counts
counts <- plotsurvtest %>% select(Species, Site, Plot, Neighbours01, surv_to_produce_seeds, n)
#Making Neighbours01 and survival characters
test100$Neighbours01 <- as.factor(test100$Neighbours01)
test100$surv_to_produce_seeds <- as.factor(test100$surv_to_produce_seeds)
counts$Neighbours01 <- as.factor(counts$Neighbours01)
counts$surv_to_produce_seeds <- as.factor(counts$surv_to_produce_seeds)
test5 <- left_join(test100, counts)
#Making NAs 0 (this is necessary to correctly calculate proportions, otherwise returns NAs)
test5 <- within(test5, n[is.na(n)] <- '0')
#Calculate proportions
test5$n <- as.numeric(test5$n)
test6 <- test5 %>% group_by(Species, Site, Plot, Neighbours01)  %>% 
  summarise(prop_survival = n[surv_to_produce_seeds==1]/(n[surv_to_produce_seeds==1] + n[surv_to_produce_seeds==0]))

#Create identifier for Site|Plot|Neighbours01
#plotsurvtest2 <- plotsurvtest2 %>% unite("id", Species:Site:Plot:Neighbours01, remove = "false")

test60 <- test6 %>% group_by(Species, Site, Plot, Neighbours01) %>% 
  pivot_wider(names_from='Neighbours01', values_from='prop_survival')
test60 <- test60 %>% select(Species, Site, Plot, surv_prop_no_neighbours='1', surv_prop_neighbours='0')

testing <- left_join(survpop, test60)
## They are very different - why??

par(mfrow=c(1,2))
ggplot(testing, aes(x= surv_prop_no_neighbours, y = surv_means_no_neighbours))+
  geom_point(alpha=0.4)+
  theme_classic()+
  facet_wrap(~Species)
ggplot(testing, aes(x= surv_prop_neighbours, y = surv_means_neighbours))+
  geom_point(alpha=0.4)+
  theme_classic()+
  facet_wrap(~Species)

### Calculating population growth values as lambda_no_nbh and lambda_nbh
#Need seed survival to be species-specific, adding in seed survival values
#speciestable is where seed survival data is
popgrowthratedata <- left_join(popgrowthratedata, speciestable)

#Going to split the dataframe and then merge so I can pivot seeds and survival separately with the same names_to column
survtopivot <- popgrowthratedata %>% select(Species, Site, Plot, plot_germ, seed_survival, surv_means_no_neighbours, surv_means_neighbours)
seedstopivot <- popgrowthratedata %>% select(Species, Site, Plot, seed_means_no_neighbours, seed_means_neighbours)

#pivot_longer survival
survtopivot <- survtopivot %>% pivot_longer(cols = c(surv_means_no_neighbours, surv_means_neighbours), names_to = "neighbours01", values_to = "plot_survival")
#Renaming neighbours01 values to neighbours0 and neighbours1
survtopivot <- within(survtopivot, neighbours01[neighbours01 == "surv_means_no_neighbours"] <- "neighbours0")
survtopivot <- within(survtopivot, neighbours01[neighbours01 == "surv_means_neighbours"] <- "neighbours1")

#Pivoting seeds
seedstopivot <- seedstopivot %>% pivot_longer(cols = c(seed_means_no_neighbours, seed_means_neighbours), names_to = "neighbours01", values_to = "plot_fecundity")
#Renaming
seedstopivot <- within(seedstopivot, neighbours01[neighbours01 == "seed_means_no_neighbours"] <- "neighbours0")
seedstopivot <- within(seedstopivot, neighbours01[neighbours01 == "seed_means_neighbours"] <- "neighbours1")

#Merge them back together
poplongdata <- left_join(survtopivot, seedstopivot)

#Calculate population growth rates
# Per capita growth rate of a given population  i = seed survival*(1-germination)+number of viable seeds produced per germinant*germination
# Isaac's (1-germ)*seed survival + rate of germ*prob of survival to reproductive maturity*seed production of survivors
# prob of germination / germination fraction. # germinated / total germination, currently I have this as percent_germ which is a proportion (not percentage, despite the name)
poplongdata <- poplongdata %>% group_by(neighbours01) %>% mutate(lambda = seed_survival*(1-plot_germ)+plot_fecundity*plot_survival*plot_germ)
#Making a column for plotid in same format as other dataframes
poplongdata <- poplongdata %>% unite("plotid", Site:Plot, remove = "false")

#To analyse this, need to bring in other plot-level data - PC1, PC2 values.
#Creating dataframe with PCA values, trait values and Treatment info to merge
pcaplot <- vitaldata %>% select(Site, Plot, plotid, PC1, PC2, std_PC1, std_PC2, Treatment, std_log_SLA, std_LDMC, std_log_D13C) %>% 
  group_by(plotid) %>% filter(row_number() == 1) %>% select(plotid, PC1, PC2, std_PC1, std_PC2, Treatment, std_log_SLA, std_LDMC, std_log_D13C)
popdata <- left_join(poplongdata, pcaplot, by = "plotid")

#What is the distribution of lambda data?
#left-skewed, log it
#hist(popdata$lambda)
#hist(log(popdata$lambda)+1)
#Create column for log(lambda)+1
popdata <- popdata %>% mutate(log_lambda_p1 = log(lambda)+1)

# For modelling, want Ambient to the be reference - coming across from vitaldata with dry first
popdata$Treatment <- factor(popdata$Treatment, level = c("Ambient", "Dry", "Wet"))

#Make a dataframe per species
lambdaarca <- popdata %>% filter(Species == "ARCA")
lambdahygl <- popdata %>% filter(Species == "HYGL")
lambdalaro <- popdata %>% filter(Species == "LARO")
lambdapeai <- popdata %>% filter(Species == "PEAI")
lambdaplde <- popdata %>% filter(Species == "PLDE")
lambdapole <- popdata %>% filter(Species == "POLE")
lambdatrcy <- popdata %>% filter(Species == "TRCY")
lambdatror <- popdata %>% filter(Species == "TROR")
lambdavero <- popdata %>% filter(Species == "VERO")

species.list.l<-list(lambdaarca, lambdahygl, lambdalaro, lambdapeai, lambdaplde, lambdapole, lambdatrcy, lambdatror, lambdavero)

#Reordering watering treatments to  Dry, Ambient, Wet for plotting
#popdata$Treatment <- factor(popdata$Treatment, level = c("Dry", "Ambient", "Wet"))

###03/08/2022
### Individual level random effects can help with overdispersion
# Try adding row ID for each subplot
test <- vitaldata %>% rownames_to_column()

