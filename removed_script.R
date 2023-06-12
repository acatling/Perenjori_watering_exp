### Perenjori Watering Exp removed code

#### Code removed from the data_preparation doc
##### not needed #####
#Assigning Neighbour values to C, E and T NAs. Removing this for now, letting survey data determine this
#Test: 2 A 9, 1 VERO T1 has Neighbours (1) --> should stay one
#mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'C' & is.na(Neighbours)] <- '1')
#mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'T' & is.na(Neighbours)] <- '0')
#mortalitydata <- within(mortalitydata, Neighbours[C_E_or_T == 'E' & is.na(Neighbours)] <- '0')

#Plots with neighbours - C. Plots without neighbours - E and T
#DO THIS - for mortality dataset. removing this, basing neighbour presence off surveys.
#seeddata$Neighbours <- ifelse(seeddata$C_E_or_T == "C", "Yes", "No")
#seeddata$Neighbours <- factor(seeddata$Neighbours)

#Check- do I need this row below? not kept in seeddatatrim
seeddata <- seeddata %>% mutate(logp1_viable_seeds = log(No_viable_seeds_grouped+1))

#test <- surveydataraw %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(row_number() == 1)
#1065 subplots surveyed

## Looking into why there are discrepancies more closely:
#So basically I want to use left join to keep all the rows from listsubplots,
# but then I also need to check which rows were missed from the dataframe I am adding in, and fix that
#library(arsenal)
#test <- comparedf(listsubplots, survtomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"), factor.as.char = TRUE)
#summary(test)
# The 'observations not shared' table shows rows in data frame x (listsubplots) that aren't in y (mortalitydatatrim) and vice versa. y is at the bottom

###Scaling the number of seeds produced across species to a percentage of maximum produced
#Need these values to be integers for negative binomial models
#note that below isn't working
#vitaldata <- vitaldata %>% group_by(Species) %>% mutate(seeds_percent = round(No_viable_seeds_grouped/max(No_viable_seeds_grouped)*100))

## See previous versions for old script calculating survival and seed production from models

#this was resolved:
## Problem here: Getting 0 from top code, a 1 from bottom code if the mean is 0.
#Because exp(0) = 1. Fixed, except:
#Top code, ARCA 1 C neighbours is an NA. Bottom code, it's a 0.
#Fine - bottom code is what is needed for lambda anyway
#Top code, ARCA 1 C no neighbours is 0. Bottom code, 0. But it should be exp(mean(c(log(1),0,0))) = 1. Or flat mean 0.333
#ARCA 1 C no neighbours: mean(c(1,0,0)) = 0.333
# ARCA 1 C neighbours: NA. None.

# Isaac's (1-germ)*seed survival + rate of germ*prob of survival to reproductive maturity*seed production of survivors


#exp(0) = 1 fixes problem where it was 1 seed, log(1)=0, then exp(0)=1.
##### not used ####
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

### when merging datasets together:
### Adding in diversity data from surveys
vitaldata <- left_join(vitaldata, SDItomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"))
vitaldata <- left_join(vitaldata, shannontomerge, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep"))
#Make SDI and shannon NA rows 0.
#NAs are where there are no neighbours (or just dodder for SDI)
vitaldata <- within(vitaldata, SDI[is.na(SDI)] <-  '0')
vitaldata <- within(vitaldata, shannon[is.na(shannon)] <-  '0')
#Dodder is currently NA for plots that weren't surveyed, making it 0 for these (I recorded where it was present, so for everywhere I didn't record it it's a 0)
vitaldata <- within(vitaldata, Dodder01[is.na(Dodder01)] <-  '0')

## Need SDI and shannon to be numeric
vitaldata$SDI <- as.numeric(vitaldata$SDI)
vitaldata$shannon <- as.numeric(vitaldata$shannon)

datanonly$SDI <- as.numeric(datanonly$SDI)
datanonly$shannon <- as.numeric(datanonly$shannon)

#Creating dataframe with trait values at plot-level to merge
plottraits <- vitaldata %>% select(Site, Plot, plotid, PC1, PC2, std_PC1, std_PC2, Treatment, std_log_SLA, std_LDMC, std_log_D13C) %>% 
  group_by(plotid) %>% filter(row_number() == 1) %>% select(plotid, PC1, PC2, std_PC1, std_PC2, Treatment, std_log_SLA, std_LDMC, std_log_D13C)
popdata <- left_join(poplongdata, pcaplot, by = "plotid")


##### Still need to check these - code still in data prep script #####
# From survey notes PEAI T2 became C2 so I am changing this above in the mortality dataset (already changed for survey and seeds)
#Check this later for the germination data, will have to change it in when I merge all vital rates because has a diff number of seeds sown for T2 and C2*
mortalitydata <- within(mortalitydata, C_E_or_T[Site == '1' & Plot == 'A' & Species == 'PEAI' & C_E_or_T == 'T' & Rep == '2'] <- 'C')

#Check this later - where 3 A NA NA NA data belongs 'not sure which plot this is', does not merge into vitaldata
#Not sure what TROR E 4 is, check this later, does not merge into vitaldata

## Check and fix 3 C PEAI C1 which has Cheilanthes and UNK 8 as NA neighbour count
## Also check and fix 7 C PEAI C2 which has UNK 1 as NA neighbour count

# Note that it is not merging in the six unknown plots in seeddatatrim, fine for now

# Subplots in surveytrim that aren't in vitaldata: 1 A LARO T1 (fixed, check), 2 B HYGL T? 1? (fixed), 2 C PLDE-POLE T 2 (fixed), 3 A NA NA NA (check later), 4 A TROR E 4 (check later)

#removed this comment:
#Check why some plots weren't surveyed - I thought I surveyed everything that germinated
##mortalitydatatrim$Neighbours <- as.factor(mortalitydatatrim$Neighbours)
#test2 <- anti_join(mortalitydatatrim, surveytrim, by = c("Site", "Plot", "Species", "C_E_or_T", "Rep")) #107 rows

#### Code removed from paper_figures ####
