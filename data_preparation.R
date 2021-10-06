######### Putting all the data prep and notes here

library(tidyverse)

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

######## CHECK THIS: Plots that were never thinned:
#Site 4 A PEAI E1; Site 7 B VERO T 2; Site 2 C POLE T 9 (not thinned, actually POLE);
#Site 5 C ARCA T1 (died before thinning); Site 2 A VERO T 2; Site 2 B PEAI T 2; Site 4 A PEAI E 1
#Site 4 B ARCA T1

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
mortalitydatatrim <- mortalitydatatrim %>% rename(Neighbours01 = "Neighbours")
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
seedmortality$ProducedSeeds <- as.factor(seedmortality$ProducedSeeds)
#It works! 609 produced seeds (at least one inviable or viable), and 534 didn't produce any

################# Community survey data ################
surveydataraw <- read_csv("Data/community_surveys_data.csv")
surveydataraw$Site <- factor(surveydataraw$Site)
surveydataraw$Species <- factor(surveydataraw$Species)
#Removing subplots
surveydata <- anti_join(surveydataraw, to_remove_start)
surveydata <- anti_join(surveydataraw, to_remove_after_germ)
#Merging with seed count data
# Probably don't need to merge yet, but keeping it like this
seedsurveys <- full_join(seeddata, surveydata)
#Replacing Neighbour_count NAs where Neighbour = Yes with 0. e.g. so Dodder now is counted as 0
# This is because I can't tally with NAs and treating Dodder as 0 because it is in my no-comp. subplots
seedsurveys <- within(seedsurveys, Neighbour_count[Neighbour_sp == 'Dodder'] <- '0')
seedsurveys$Neighbour_count <- as.numeric(seedsurveys$Neighbour_count)
#Renaming E and T Neighbours to Yes that actually have neighbours (neighbour_count >= 1)
seedsurveys$Neighbours <- ifelse(seedsurveys$Neighbour_count >= 1, "Yes", "No")
#Also want to plot the subplots with 0 neighbour_count, so need NAs to be 0s
#This filters to subplots without Neighbour_count (currently NAs, same with Neighbours)
surveytest <- seedsurveys %>% filter(!complete.cases(Neighbour_count)) 
surveytest$Neighbour_count <- "0"
surveytest$Neighbours <- "No"
surveytest2 <- seedsurveys %>% filter(complete.cases(Neighbour_count))
#Checking numbers are right - 1740 lines in seedsurveys, 633 !complete.cases, 1107 complete.cases. Phew!
surveytest$Neighbour_count <- as.numeric(surveytest$Neighbour_count)
seedsurveysfixed <- full_join(surveytest, surveytest2)
#I shouldn't have to do the below in two steps but I can't figure out one, and it works
neighbourabundance <- seedsurveysfixed %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  summarise(Total_abundance = sum(Neighbour_count))
seedabundance <- merge(neighbourabundance, seedsurveysfixed)
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
allabundanceseeds <- seedabundance %>% group_by(Site, Plot, Species, C_E_or_T, Rep, 
                                                Total_abundance,Intra_abundance, Inter_abundance)%>%
  filter(row_number() == 1)
########### TO DO? - Change this below bit in light of ProducedSeeds
#Trimming down to just what is required to merge with seedmortality
surveytrim <- allabundanceseeds %>% select(Site, Plot, Species, C_E_or_T, Rep, Neighbours,
                                           Total_abundance, Intra_abundance, Inter_abundance)
#There are 50 rows that say No to neighbours but have a neighbour total abundance >0
surveytrim <- within(surveytrim, Neighbours[Total_abundance > "0"] <- "Yes")
surveytrimn <- surveytrim %>% filter(Neighbours == "Yes")
#### All data together for Full_model_WA sheet ####
dataall <- left_join(seedmortality, surveytrimn)
dataall <- within(dataall, Neighbours[is.na(Total_abundance)] <- "No")
dataall <- dataall %>% replace(is.na(.), 0)
dataall$Survival <- as.factor(dataall$Survival)
datanotscaled <- dataall
#Standardising continuous explanatory variables to a mean of 0 and SD of 1
dataall$cc_percentage <- scale(dataall$cc_percentage, center = TRUE, scale = TRUE)
dataall$Total_abundance <- scale(dataall$Total_abundance, center = TRUE, scale = TRUE)
dataall$Inter_abundance <- scale(dataall$Inter_abundance, center = TRUE, scale = TRUE)
dataall$Intra_abundance <- scale(dataall$Intra_abundance, center = TRUE, scale = TRUE)
#Splitting data by species
arcamortality <- dataall %>% filter(Species == "ARCA")
hyglmortality <- dataall %>% filter(Species == "HYGL")
laromortality <- dataall %>% filter(Species == "LARO")
peaimortality <- dataall %>% filter(Species == "PEAI")
pldemortality <- dataall %>% filter(Species == "PLDE")
polemortality <- dataall %>% filter(Species == "POLE")
trcymortality <- dataall %>% filter(Species == "TRCY")
trormortality <- dataall %>% filter(Species == "TROR")
veromortality <- dataall %>% filter(Species == "VERO")
