######### Putting all the data prep and notes here

library(tidyverse)
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
canopydatatrim <- canopydata %>% select(Site, Plot, cc_percentage)
germinationdata$Site <- as.factor(germinationdata$Site)

################# Survival/mortality data #############
#Importing mortality data 
mortalitydataraw <- read_csv("Data/mortalityofgerminated.csv")
#Need to rename 8 C VERO T 1 in mortality data to E1 since note says 'stick says E1'. Already changed in seed data
mortalitydata <- within(mortalitydataraw, C_E_or_T[Site == '8' & Plot == 'C' & Species == 'VERO' & C_E_or_T == 'T'] <- 'E')
# From survey notes PEAI T2 became C2 so I am changing this above in the mortality dataset (already changed for survey and seeds)
#Check this later for the germination data, will have to change it in when I merge all vital rates because has a diff number of seeds sown for T2 and C2
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
#Assigning Neighbour values to C, E and T NAs
#Test: 2 A 9, 1 VERO T1 has Neighbours (1) --> should stay one
mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'C' & is.na(Neighbours)] <- '1')
mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'T' & is.na(Neighbours)] <- '0')
mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'E' & is.na(Neighbours)] <- '0')
#CHECK THIS: Should I assign Site 1 C TROR C 1 a 0 for survival? Died
#Where no_viable_seeds and no_viable_seeds > 0, SeedsProduced = 1, otherwise 0.
#For the survival dataset, will look at SeedsProduced 1/0.
mortalitydatatrim <- mortalitydata %>% select(Site, Plot, Species, C_E_or_T, Rep,
                                              No_germinated,Survival, Neighbours)
mortalitydatatrim$Site <- as.factor(mortalitydatatrim$Site)
mortalitydatatrim$Rep <- as.factor(mortalitydatatrim$Rep)

################# Seed count data ################
##DO THIS - 
#Extra seeds that were collected (e.g. broke of by accident) and
#flower counts where collected early or seeds missed
#flowers_seeds_to_add <- read_csv("Data/flowers_seeds_to_add.csv")
##
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
#Plots with neighbours - C. Plots without neighbours - E and T
#DO THIS - for mortality dataset
seeddata$Neighbours <- ifelse(seeddata$C_E_or_T == "C", "Yes", "No")
seeddata$Neighbours <- factor(seeddata$Neighbours)
#Making a new column for logged viable seeds +1
seeddata <- seeddata %>% mutate(log1_viable_seeds = log(No_viable_seeds_grouped+1))
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
#I shouldn't have to do the below in two steps but I can't figure out one, and it works
neighbourabundance <- surveydatafixed %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  summarise(Total_abundance = sum(Neighbour_count))
surveyabundance <- left_join(surveydatafixed, neighbourabundance)
#Making a column for whether or not Dodder was present in subplot (1 = yes, 0 = no)
#Then filtering to one row per subplot
#Creating the below dataframe seedabundancerows to use later
surveyabundancerows <- surveyabundance
surveyabundance <- surveyabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
                  mutate(Dodder01 = case_when(any(Neighbour_sp == "Dodder") ~ "1",
                         TRUE ~ "0")) %>% filter(row_number() == 1)
#Renaming E and T Neighbours to Yes that actually have neighbours (neighbour_count >= 1)
surveyabundance$Neighbours <- ifelse(surveyabundance$Total_abundance > 0, "Yes", "No")
#Trimming down to just what is required to merge
surveytrim <- surveyabundance %>% select(Site, Plot, Species, C_E_or_T, Rep, Neighbours,
                                           Total_abundance, Intra_abundance, Inter_abundance, Dodder01)
#To deal with NA problems in scaling abund values, I will change NAs to 0
#Have to remove one Species plot that is NA
surveytrim <- surveytrim %>% filter(!is.na(Species))
surveytrim <- surveytrim %>% replace(is.na(.),0)
#Creating dataframe with neighbour info only
surveytrimn <- surveytrim %>% filter(Neighbours == "Yes")
#1057 subplots, 517 of which have neighbours

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
survtomerge <- mortalitydatatrim %>% select(Site, Plot, Species, C_E_or_T, Rep, Survival, Neighbours)
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
#604 produced seeds (at least one inviable or viable), and 551 didn't produce any
#Want to use cbind for number of successes and failures, using binomial glm.
#Create a row for total number germinated, total number didn't
germcalcs <- germcalcs %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(total_germ = sum(Number_germinated+February_germination),
         total_no_germ = Number_seeds_sown-total_germ,
         percent_germ = total_germ/Number_seeds_sown)
###Scaling the number of seeds produced across species to a percentage of maximum produced
#Need these values to be integers for negative binomial models
germcalcs <- germcalcs %>% group_by(Species) %>% mutate(seeds_percent = round(No_viable_seeds_grouped/max(No_viable_seeds_grouped)*100))
germcalcstomerge <- germcalcs %>% select(Site, Plot, Species, C_E_or_T, Rep, surv_to_produce_seeds, total_germ, total_no_germ, percent_germ)
#Merge back into big dataframe
vitaltomerge <- vitaldata %>% select(-c(Notes_germination, Survival, No_viable_seeds_grouped, No_inviable_seeds_grouped))
vitaldata <- left_join(vitaltomerge, germcalcstomerge, by =c("Site", "Plot", "Species", "C_E_or_T", "Rep"))

#Merging with information on treatments (watering and cover (shade/sun) per plot)
treatments <- read_csv("Data/treatments_meta.csv")
treatments$Site <- as.factor(treatments$Site)
vitaldata <- left_join(vitaldata, treatments)

### Checking distributions and standardising data ##
#Distributions of traits and PC1, have to use datanotscaled to log
#SLA better logged?
#LDMC still strange logged
#WUE/D13C better logged but still left-skewed
#PC1 - right-skewed, logging doesn't help
with(vitaldata, pairs(~SLA + LDMC + mean_D13C + PC1, diag.panel = panel.hist))
with(vitaldata, pairs(~log(SLA+2) + LDMC + log(mean_D13C+2) + PC1, diag.panel = panel.hist))
#Logging and standardising variables
vitaldata$log_SLA <-log(vitaldata$SLA)
vitaldata$log_D13C <-log(vitaldata$mean_D13C)

#Standardising continuous explanatory variables to a mean of 0 and SD of 1
vitaldata$std_cc <- scale(vitaldata$cc_percentage, center = TRUE, scale = TRUE)
vitaldata$std_PC1 <- scale(vitaldata$PC1, center = TRUE, scale = TRUE)
vitaldata$std_PC2 <- scale(vitaldata$PC2, center = TRUE, scale = TRUE)
vitaldata$std_PC3 <- scale(vitaldata$PC3, center = TRUE, scale = TRUE)
vitaldata$std_log_SLA <- scale(vitaldata$log_SLA, center = TRUE, scale = TRUE)
vitaldata$std_LDMC <- scale(vitaldata$LDMC, center = TRUE, scale = TRUE)
vitaldata$std_log_D13C <- scale(vitaldata$log_D13C, center = TRUE, scale = TRUE)

# Transforming explanatory variables
vitaldata$logp1_totalabund <- log(vitaldata$Total_abundance+1)
vitaldata$logp1_interabund <- log(vitaldata$Inter_abundance+1)
vitaldata$logp1_intraabund <-log(vitaldata$Intra_abundance+1)
vitaldata$log_SLA <-log(vitaldata$SLA)
vitaldata$log_D13C <-log(vitaldata$mean_D13C)
#Standardising continuous explanatory variables to a mean of 0 and SD of 1
vitaldata$std_cc <- scale(vitaldata$cc_percentage, center = TRUE, scale = TRUE)
vitaldata$std_logp1_totalabund <- scale(vitaldata$logp1_totalabund, center = TRUE, scale = TRUE)
vitaldata$std_logp1_interabund <- scale(vitaldata$logp1_interabund, center = TRUE, scale = TRUE)
vitaldata$std_logp1_intraabund <- scale(vitaldata$logp1_intraabund, center = TRUE, scale = TRUE)

#Need surv_to_produce_seeds to be numeric, not a factor to plot curves
vitaldata$surv_to_produce_seeds <- as.numeric(vitaldata$surv_to_produce_seeds)

#Renaming Control watering treatment to Ambient
vitaldata <- vitaldata %>% mutate(Treatment = recode(Treatment, Control = 'Ambient'))

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
