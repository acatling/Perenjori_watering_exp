### Stocktake of germination July 2020
# C (competition), E (environment) and T (trait) subplots
#This is telling me how many had zero germination

library(tidyverse)
library(plyr)
library(readr)

site1A <-read_csv("Data/Site_1_A_germination.csv")
tallysite1A <- site1A %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site1B <-read_csv("Data/Site_1_B_germination.csv")
tallysite1B <- site1B %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site1C <-read_csv("Data/Site_1_C_germination.csv")
tallysite1C <- site1C %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site2A <-read_csv("Data/Site_2_A_germination.csv")
head(site2A)
tallysite2A <- site2A %>% group_by(Species, C_E_or_T, No_germinants_210720) %>% tally() %>% filter(No_germinants_210720 == "0")

site2B <-read_csv("Data/Site_2_B_germination.csv")
head(site2B)
tallysite2B <- site2B %>% group_by(Species, C_E_or_T, No_germinants_210720) %>% tally() %>% filter(No_germinants_210720 == "0")

site2C <-read_csv("Data/Site_2_C_germination.csv")
head(site2C)
tallysite2C <- site2C %>% group_by(Species, C_E_or_T, No_germinants_210720) %>% tally() %>% filter(No_germinants_210720 == "0")

site3A <-read_csv("Data/Site_3_A_germination.csv")
head(site3A)
tallysite3A <- site3A %>% group_by(Species, C_E_or_T, No_germinants_190720) %>% tally() %>% filter(No_germinants_190720 == "0")

site3B <-read_csv("Data/Site_3_B_germination.csv")
head(site3B)
tallysite3B <- site3B %>% group_by(Species, C_E_or_T, No_germinants_190720) %>% tally() %>% filter(No_germinants_190720 == "0")

site3C <-read_csv("Data/Site_3_C_germination.csv")
head(site3C)
tallysite3C <- site3C %>% group_by(Species, C_E_or_T, No_germinants_190720) %>% tally() %>% filter(No_germinants_190720 == "0")

site4A <-read_csv("Data/Site_4_A_germination.csv")
tallysite4A <- site4A %>% group_by(Species, C_E_or_T, No_germinants_180720) %>% tally() %>% filter(No_germinants_180720 == "0")

site4B <-read_csv("Data/Site_4_B_germination.csv")
tallysite4B <- site4B %>% group_by(Species, C_E_or_T, No_germinants_180720) %>% tally() %>% filter(No_germinants_180720 == "0")

site4C <-read_csv("Data/Site_4_C_germination.csv")
tallysite4C <- site4C %>% group_by(Species, C_E_or_T, No_germinants_180720) %>% tally() %>% filter(No_germinants_180720 == "0")

site5A <-read_csv("Data/Site_5_A_germination.csv")
tallysite5A <- site5A %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site5B <-read_csv("Data/Site_5_B_germination.csv")
tallysite5B <- site5B %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site5C <-read_csv("Data/Site_5_C_germination.csv")
tallysite5C <- site5C %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site6A <-read_csv("Data/Site_6_A_germination.csv")
head(site6A)
tallysite6A <- site6A %>% group_by(Species, C_E_or_T, No_germinants_200720_210720) %>% tally() %>% filter(No_germinants_200720_210720 == "0")

site6B <-read_csv("Data/Site_6_B_germination.csv")
head(site6B)
tallysite6B <- site6B %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site6C <-read_csv("Data/Site_6_C_germination.csv")
head(site6C)
tallysite6C <- site6C %>% group_by(Species, C_E_or_T, No_germinants_200720) %>% tally() %>% filter(No_germinants_200720 == "0")

site7A <-read_csv("Data/Site_7_A_germination.csv")
head(site7A)
tallysite7A <- site7A %>% group_by(Species, C_E_or_T, No_germinants_180720) %>% tally() %>% filter(No_germinants_180720 == "0")

site7B <-read_csv("Data/Site_7_B_germination.csv")
head(site7B)
tallysite7B <- site7B %>% group_by(Species, C_E_or_T, No_germinants_180720) %>% tally() %>% filter(No_germinants_180720 == "0")

site7C <-read_csv("Data/Site_7_C_germination.csv")
head(site7C)
tallysite7C <- site7C %>% group_by(Species, C_E_or_T, No_germinants_180720) %>% tally() %>% filter(No_germinants_180720 == "0")

site8A <-read_csv("Data/Site_8_A_germination.csv")
head(site8A)
tallysite8A <- site8A %>% group_by(Species, C_E_or_T, No_germinants_190720) %>% tally() %>% filter(No_germinants_190720 == "0")

site8B <-read_csv("Data/Site_8_B_germination.csv")
head(site8B)
tallysite8B <- site8B %>% group_by(Species, C_E_or_T, No_germinants_190720) %>% tally() %>% filter(No_germinants_190720 == "0")

site8C <-read_csv("Data/Site_8_C_germination.csv")
head(site8C)
tallysite8C <- site8C %>% group_by(Species, C_E_or_T, No_germinants_190720) %>% tally() %>% filter(No_germinants_190720 == "0")

#### Filtering species to look at how many subplots

sp <- site1A %>% filter(Species == "HYGL")
#view(sp)

# Plotting the number of subplots with neighbours, alone and for traits for each species

subplotdata <-read_csv("Data/Stocktake_subplots_August_2020.csv")
head(subplotdata)

ARCAsub <- subplotdata %>% filter(Species == "ARCA")
view(ARCAsub)
HYGLsub <- subplotdata %>% filter(Species == "HYGL")
view(HYGLsub)
LAROsub <- subplotdata %>% filter(Species == "LARO")
view(LAROsub)
PEAIsub <- subplotdata %>% filter(Species == "PEAI")
view(PEAIsub)
PLDEsub <- subplotdata %>% filter(Species == "PLDE")
view(PLDEsub)
POLEsub <- subplotdata %>% filter(Species == "POLE")
view(POLEsub)
TRCYsub <- subplotdata %>% filter(Species == "TRCY")
view(TRCYsub)
TRORsub <- subplotdata %>% filter(Species == "TROR")
view(TRORsub)
VEROsub <- subplotdata %>% filter(Species == "VERO")
view(VEROsub)

#Below gives me a count for all. I need a grouped bar chart where site/plot are different Xs. Also need to add bars for
# other points - No_alone etc.
#ggplot(data = ARCAsub) +
 # geom_bar(mapping = aes(x = No_with_neighbours))


## Merging all files to export as a csv

library(plyr)

allgerminationdata <- rbind.fill(site1A, site1B, site1C, site2A, site2B, site2C, site3A, site3B, site3C,  site4A, site4B, site4C, 
                       site5A, site5B, site5C,  site6A, site6B, site6C,  site7A, site7B, site7C,  site8A, site8B, site8C)
#Making a list of all subplots
listsubplots <- allgerminationdata %>% select(1:5)
write_csv(listsubplots, "Output/listsubplots.csv")

#Below is equivalent to rbind.fill above
#Need to format all the Notes as factors before I can merge them
# str(site1A)
# str(site1B)
# str(site1C)
# site1C$No_germinants_200720 <- as.factor(site1C$No_germinants_200720)
# site1C$No_germinants_290720 <- as.factor(site1C$No_germinants_290720)
# str(site2A)
# site2A$No_germinants_130820 <- as.factor(site2A$No_germinants_130820)
# site2A$No_germinants_210720 <- as.factor(site2A$No_germinants_210720)
# str(site2B)
# str(site2C)
# site2C$No_germinants_210720 <- as.factor(site2C$No_germinants_210720)
# site2C$No_germinants_130820 <- as.factor(site2C$No_germinants_130820)
# str(site3A)
# site3A$No_germinants_130820 <- as.factor(site3A$No_germinants_130820)
# str(site3B)
# site3B$No_germinants_130820 <- as.factor(site3B$No_germinants_130820)
# site3B$No_germinants_190720 <- as.factor(site3B$No_germinants_190720)
# str(site3C)
# str(site4A)
# site4A$No_germinants_180720 <- as.factor(site4A$No_germinants_180720)
# site4A$No_germinants_130820 <- as.factor(site4A$No_germinants_130820)
# str(site4B)
# site4B$No_germinants_180720 <- as.factor(site4B$No_germinants_180720)
# site4B$No_germinants_130820 <- as.factor(site4B$No_germinants_130820)
# str(site4C)
# site4C$No_germinants_180720 <- as.factor(site4C$No_germinants_180720)
# site4C$No_germinants_130820 <- as.factor(site4C$No_germinants_130820)
# str(site5A)
# str(site5B)
# site5B$No_germinants_200720 <- as.factor(site5B$No_germinants_200720)
# site5B$No_germinants_130820 <- as.factor(site5B$No_germinants_130820)
# str(site5C)
# str(site6A)
# site6A$No_germinants_130820 <- as.factor(site6A$No_germinants_130820)
# str(site6B)
# site6B$No_germinants_200720 <- as.factor(site6B$No_germinants_200720)
# site6B$No_germinants_130820 <- as.factor(site6B$No_germinants_130820)
# str(site6C)
# site6C$No_germinants_200720 <- as.factor(site6C$No_germinants_200720)
# str(site7A)
# site7A$No_germinants_130820 <- as.factor(site7A$No_germinants_130820)
# site7A$No_germinants_180720 <- as.factor(site7A$No_germinants_180720)
# str(site7B)
# site7B$No_germinants_130820 <- as.factor(site7B$No_germinants_130820)
# site7B$No_germinants_180720 <- as.factor(site7B$No_germinants_180720)
# str(site7C)
# site7C$No_germinants_130820 <- as.factor(site7C$No_germinants_130820)
# site7C$No_germinants_180720 <- as.factor(site7C$No_germinants_180720)
# str(site8A)
# str(site8B)
# site8B$No_germinants_190720 <- as.factor(site8B$No_germinants_190720)
# site8B$No_germinants_130820 <- as.factor(site8B$No_germinants_130820)
# str(site8C)
# site8C$No_germinants_130820 <- as.factor(site8C$No_germinants_130820)
# 
# allgerminationdata <- bind_rows(site1A, site1B, site1C, site2A, site2B, site2C, site3A, site3B, site3C,  site4A, site4B, site4C, 
#                                 site5A, site5B, site5C,  site6A, site6B, site6C,  site7A, site7B, site7C,  site8A, site8B, site8C)

## Filtering checklist to just things that germinated (removing the 0s)
checklist <- read_csv("Data/Manually_edited/Checklist_germinated_plants.csv")
onlygerminatedsubplots <- checklist %>% filter(Germinate != "0")
write_csv(onlygerminatedsubplots, "Output/onlygerminatedsubplots.csv")

#Note that the edited onlygerminated has some more observations that the non-edited one because I hadn't updated all the
# germination data and manually imported some rows
remaining <-read_csv("Data/Manually_edited/onlygerminatedsubplotsedited.csv")
#NAs are a bit tricky to filter but this works:
trueremaining <- remaining %>% filter(is.na(Done))
write_csv(trueremaining, "Output/remainingchecklist.csv")

