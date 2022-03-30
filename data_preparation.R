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
germinationdata$Rep <- as.numeric(germinationdata$Rep)
######## To be removed for various reasons:
germinationdata <- anti_join(germinationdata, to_remove_start)
#Merge with canopy data
canopylitterdata <- read_csv("Data/canopy_litter_2020.csv")
canopydata <- canopylitterdata  %>% rowwise() %>%
  mutate(meancc = mean(c(`Measurement 1`, `Measurement 2`, `Measurement 3`, `Measurement 4`, `Measurement 5`)),
         cc_percentage = meancc/24*100)
canopydatatrim <- canopydata %>% select(Site, Plot, cc_percentage)
germinationdata <- merge(germinationdata, canopydatatrim)
#Create a column for plotid
germinationdata <- germinationdata %>% unite("plotid", Site:Plot, remove = "false")
germinationdata$Site <- as.factor(germinationdata$Site)
#Calculating germination rate.
## Accounting for February germination. For now, total number germinated / total number seeds
germinationdata <- germinationdata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(Germination_rate = (February_germination + Number_germinated)/Number_seeds_sown)

#First checking distribution of germination rate
#Should be logged
ggplot(germinationdata, aes(x = log(Germination_rate)+1))+
  geom_histogram()+
  theme_classic()+
  facet_wrap(~Species)

################# Survival/mortality data #############
#Importing mortality data 
mortalitydataraw <- read_csv("Data/mortalityofgerminated.csv")
#Removing subplots
mortalitydata <- anti_join(mortalitydataraw, to_remove_start)
mortalitydata <- anti_join(mortalitydataraw, to_remove_after_germ)
#Merging with treatments
treatments <- read_csv("Data/treatments_meta.csv")
mortalitydata <- full_join(mortalitydata, treatments)
#Assigning Neighbour values to C, E and T NAs
#Test: 2 A 9, 1 VERO T1 has Neighbours (1) --> should stay one
mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'C' & is.na(Neighbours)] <- '1')
mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'T' & is.na(Neighbours)] <- '0')
mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'E' & is.na(Neighbours)] <- '0')
#CHECK THIS: Should I assign Site 1 C TROR C 1 a 0 for survival? Died
#Merge with canopy data cleaned above
mortalitydata <- merge(mortalitydata, canopydatatrim)
#Create a column for plotid
mortalitydata <- mortalitydata %>% unite("plotid", Site:Plot, remove = "false")

######## Plots that were never thinned, checked:
#Site 4 A PEAI E1 (yes neighbours, all good); Site 7 B VERO T 2 (no neighbours, all good); Site 2 C POLE T 9 (not thinned, actually POLE; can't find this, must have renamed);
#Site 5 C ARCA T1 (died before thinning, all good); Site 2 A VERO T 2; Site 2 B PEAI T 2 (yes neighbours, all good); Site 4 A PEAI E 1 (yes neighbours, all good)
#Site 4 B ARCA T1 (no neighbours, all good)

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
#Merging with treatments
seeddata <- merge(treatments, seeddatarawgrouped)
#Reordering in order of dry, control, wet.
level = c("Dry", "Control", "Wet")
seeddata$Treatment <- factor(seeddata$Treatment, level = c("Dry", "Control", "Wet"))
#Need Site, Treatment and Species to be factors
seeddata$Site <- factor(seeddata$Site)
seeddata$Treatment <- factor(seeddata$Treatment)
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

##### Merging mortality data and seed count data ####
#Where no_viable_seeds and no_viable_seeds > 0, SeedsProduced = 1, otherwise 0.
#For the survival dataset, will look at SeedsProduced 1/0.
mortalitydatatrim <- mortalitydata %>% select(plotid, Site, Plot, Species, C_E_or_T, Rep, Treatment,
                                              Cover, No_germinated,Survival, Neighbours, cc_percentage)
names(mortalitydatatrim)[names(mortalitydatatrim) == 'Neighbours'] <- 'Neighbours01'
# I don't think Neighbours01 is the most meaningful or up to date. Neighbours Yes or No should be used
mortalitydatatrim$Site <- as.factor(mortalitydatatrim$Site)
seeddatatrim$Site <- as.factor(seeddatatrim$Site)
mortalitydatatrim$Rep <- as.factor(mortalitydatatrim$Rep)
seeddatatrim$Rep <- as.factor(seeddatatrim$Rep)
seedmortality <- left_join(mortalitydatatrim, seeddatatrim)
#Replace all the NAs for seed counts with 0s - below code replaces everything in dataframe
# that is an NA with 0 - that's okay here, where only No_viable and inviable seeds have NAs
seedmortality <- seedmortality %>% replace(is.na(.),0)
#Make a new column for ProducedSeeds (1 or 0)
seedmortality <- seedmortality %>% mutate(ProducedSeeds = case_when(No_viable_seeds_grouped > "0" | No_inviable_seeds_grouped > "0" ~ "1",
                                                                   No_viable_seeds_grouped == "0" & No_inviable_seeds_grouped == "0" ~ "0"))
#It works! 609 produced seeds (at least one inviable or viable), and 534 didn't produce any

################# Community survey data ################
surveydataraw <- read_csv("Data/community_surveys_data.csv")
surveydataraw$Site <- factor(surveydataraw$Site)
surveydataraw$Species <- factor(surveydataraw$Species)
#Removing subplots
surveydata <- anti_join(surveydataraw, to_remove_start)
surveydata <- anti_join(surveydataraw, to_remove_after_germ)
#test <- surveydataraw %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(row_number() == 1)
#1065 subplots surveyed
#Merging with seed count data
# Probably don't need to merge yet, but keeping it like this
seedsurveys <- full_join(seeddata, surveydata)
# Dodder is listed under Neighbour_sp but has a Neighbour_count of NA
#Replacing Neighbour_count NAs where Neighbour_sp = Dodder with 0. e.g. so Dodder now is counted as 0
# This is because I can't tally with NAs and treating Dodder as 0 because it is in my no-comp. subplots
seedsurveys <- within(seedsurveys, Neighbour_count[Neighbour_sp == 'Dodder'] <- '0')
seedsurveys$Neighbour_count <- as.numeric(seedsurveys$Neighbour_count)
#Renaming E and T Neighbours to Yes that actually have neighbours (neighbour_count >= 1)
seedsurveys$Neighbours <- ifelse(seedsurveys$Neighbour_count >= 1, "Yes", "No")
#Also want to plot the subplots with 0 neighbour_count, so need NAs to be 0s
#This filters to subplots without Neighbour_count (currently NAs, same with Neighbours)
## Check and fix 3 C PEAI C1 which has Cheilanthes and UNK 8 as NA neighbour count
## Also check and fix 7 C PEAI C2 which has UNK 1 as NA neighbour count
surveytest <- seedsurveys %>% filter(!complete.cases(Neighbour_count)) 
surveytest$Neighbour_count <- "0"
surveytest$Neighbours <- "No"
surveytest2 <- seedsurveys %>% filter(complete.cases(Neighbour_count))
#Checking numbers are right - 2393 lines in seedsurveys, 706 !complete.cases, 1687 complete.cases.
surveytest$Neighbour_count <- as.numeric(surveytest$Neighbour_count)
seedsurveysfixed <- full_join(surveytest, surveytest2)
#I shouldn't have to do the below in two steps but I can't figure out one, and it works
neighbourabundance <- seedsurveysfixed %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  summarise(Total_abundance = sum(Neighbour_count))
seedabundance <- merge(neighbourabundance, seedsurveysfixed)
#NAs, site 3, plot  - 756 - 759 are NAs. Manually changing them to Unks and removing for now
## Will check notes later
seedabundance %>% filter(Species == "UNK")
seedabundance <- seedabundance %>% filter(Species != "UNK")
# Separating abundances of intra- and inter-specific neighbours
#Isaac's code first two lines
seedabundance$Matching <- seedabundance$Species == seedabundance$Neighbour_sp
seedabundance$Matching <- ifelse(seedabundance$Matching == TRUE, 1, 0)
seedabundance <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
  mutate(Intra_abundance = sum(Neighbour_count[Matching == "1"]),
         Inter_abundance = sum(Neighbour_count[Matching == "0"]))
#Making a column for whether or not Dodder was present in subplot (1 = yes, 0 = no)
#Then filtering to one row per subplot
#Creating the below dataframe seedabundancerows to use later
seedabundancerows <- seedabundance
seedabundance <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>%
                  mutate(Dodder01 = case_when(any(Neighbour_sp == "Dodder") ~ "1",
                         TRUE ~ "0")) %>% filter(row_number() == 1)
#Trimming down to just what is required to merge with seedmortality
surveytrim <- seedabundance %>% select(Site, Plot, Species, C_E_or_T, Rep, Neighbours,
                                           Total_abundance, Intra_abundance, Inter_abundance, Dodder01)
#There are 87 rows that say No to neighbours but have a neighbour total abundance >0
# Check this, why??
#test <- surveytrim %>% filter(Neighbours == "No" & Total_abundance > 0)
surveytrim <- within(surveytrim, Neighbours[Total_abundance > "0"] <- "Yes")
surveytrimn <- surveytrim %>% filter(Neighbours == "Yes")
#1130 subplots, 516 of which have neighbours

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
#test <- anti_join(seedmortality, surveytrim)
##53 rows that seedmortality has but surveytrim doesn't. None of these plants produced seeds
# but all did germinate
#test2 <- surveytrim %>% filter(Site == "1", Plot == "A", Species == "PEAI", C_E_or_T == "T", Rep == "2")
## Not sure why using surveytrimn before
dataall <- left_join(seedmortality, surveytrim)
dataall <- within(dataall, Neighbours[is.na(Total_abundance)] <- "No")
#Removing this NA thing atm because there are some subplots that weren't surveyed
# and I don't want to incorrectly give them 0 neighbour abundance
#dataall <- dataall %>% replace(is.na(.), 0)
#Adding in soil data
dataall <- left_join(dataall, abioticpcadata)
#Adding in trait data
traitdata <- alltraits %>% select(Species, SLA, LDMC, mean_D13C)
dataall <- left_join(dataall, traitdata)
#Distributions of traits and PC1, have to use datanotscaled to log
#SLA better logged?
#LDMC still strange logged
#WUE/D13C better logged but still left-skewed
#PC1 - right-skewed, logging doesn't help
with(dataall, pairs(~SLA + LDMC + mean_D13C + PC1, diag.panel = panel.hist))
with(dataall, pairs(~log(SLA+2) + LDMC + log(mean_D13C+2) + PC1, diag.panel = panel.hist))
# Transforming explanatory variables
dataall$logp1_totalabund <- log(dataall$Total_abundance+1)
dataall$logp1_interabund <- log(dataall$Inter_abundance+1)
dataall$logp1_intraabund <-log(dataall$Intra_abundance+1)
dataall$log_SLA <-log(dataall$SLA)
dataall$log_D13C <-log(dataall$mean_D13C)
#Standardising continuous explanatory variables to a mean of 0 and SD of 1
dataall$std_cc <- scale(dataall$cc_percentage, center = TRUE, scale = TRUE)
dataall$std_logp1_totalabund <- scale(dataall$logp1_totalabund, center = TRUE, scale = TRUE)
dataall$std_logp1_interabund <- scale(dataall$logp1_interabund, center = TRUE, scale = TRUE)
dataall$std_logp1_intraabund <- scale(dataall$logp1_intraabund, center = TRUE, scale = TRUE)
dataall$std_PC1 <- scale(dataall$PC1, center = TRUE, scale = TRUE)
dataall$std_PC2 <- scale(dataall$PC2, center = TRUE, scale = TRUE)
dataall$std_PC3 <- scale(dataall$PC3, center = TRUE, scale = TRUE)
dataall$std_log_SLA <- scale(dataall$log_SLA, center = TRUE, scale = TRUE)
dataall$std_LDMC <- scale(dataall$LDMC, center = TRUE, scale = TRUE)
dataall$std_log_D13C <- scale(dataall$log_D13C, center = TRUE, scale = TRUE)

###Scaling the number of seeds produced across species to a percentage of maximum produced
#Need these values to be integers for negative binomial models
dataall <- dataall %>% group_by(Species) %>% mutate(seeds_percent = round(No_viable_seeds_grouped/max(No_viable_seeds_grouped)*100))

#Need ProducedSeeds to be numeric, not a factor to plot curves
dataall$ProducedSeeds <- as.numeric(dataall$ProducedSeeds)
#Splitting data by species
arcadata <- dataall %>% filter(Species == "ARCA")
hygldata <- dataall %>% filter(Species == "HYGL")
larodata <- dataall %>% filter(Species == "LARO")
peaidata <- dataall %>% filter(Species == "PEAI")
pldedata <- dataall %>% filter(Species == "PLDE")
poledata <- dataall %>% filter(Species == "POLE")
trcydata <- dataall %>% filter(Species == "TRCY")
trordata <- dataall %>% filter(Species == "TROR")
verodata <- dataall %>% filter(Species == "VERO")

#Seed production model data, which is filtered to things that produced
#at least one seed (viable or inviable)
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

#### Adding trait and abiotic data to germination data ####
#Adding in info on abiotic environment and traits
germdata <- left_join(germinationdata, abioticpcadata)
germdata <- left_join(germdata, traitdata)
#Logging and standardising
germdata$log_SLA <-log(germdata$SLA)
germdata$log_D13C <-log(germdata$mean_D13C)
germdata$log_germ <- log(germdata$Germination_rate+1)
#Standardising continuous explanatory variables to a mean of 0 and SD of 1
germdata$std_cc <- scale(germdata$cc_percentage, center = TRUE, scale = TRUE)
germdata$std_PC1 <- scale(germdata$PC1, center = TRUE, scale = TRUE)
germdata$std_PC2 <- scale(germdata$PC2, center = TRUE, scale = TRUE)
germdata$std_PC3 <- scale(germdata$PC3, center = TRUE, scale = TRUE)
germdata$std_log_SLA <- scale(germdata$log_SLA, center = TRUE, scale = TRUE)
germdata$std_LDMC <- scale(germdata$LDMC, center = TRUE, scale = TRUE)
germdata$std_log_D13C <- scale(germdata$log_D13C, center = TRUE, scale = TRUE)

#Want to use cbind for number of successes and failures, using binomial glm.
#Create a row for total number germinated, total number didn't
germdata <- germdata %>% mutate(total_germ = sum(Number_germinated+February_germination))
germdata <- germdata %>% mutate(total_no_germ = Number_seeds_sown-total_germ)
germdata <- germdata %>% mutate(percent_germ = total_germ/Number_seeds_sown)

