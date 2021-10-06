### August 2021

testchecklist <- read_csv("Data/Manually_edited/Checklist_germinated_plants.csv")
testonlygerminated <- read_csv("Data/Manually_edited/onlygerminatedsubplotsedited.csv")

testing <- full_join(testchecklist, testonlygerminated, by = c("Site", "Plot", "Row", "Column"))
write.csv(testing, file = "Output/combined_checklist_onlygerminatededited.csv")

#Don't need any of the below anymore, after fixing data into new germination_final_data file
#There are two PEAI E1s (one of which was a POLE where PEAI was accidentally planted): manually remaining this one to PEAI E4 (position 3,3).
#Changed in onlygerminatedsubplotsedited, mortalityofgerminated, checklist
#They both seemingly survived? by seed_count_2020 only has PEAI E 1.
#Two HYGL C1s with different germination values, they swapped. Renaming only C1 to T2 (what it was swapped to)
#Two blanks in Checklist (coming up as duplicate rows because no C_E_or_T or rep assignment)

#Most up to date is onlygerminatedsubplotsedited, however note that this filtered out
# plants with 0 germination. Checklist minus things that didn't germinate, THEN manually updated
#Note that the edited onlygerminated has some more observations that the non-edited one because I hadn't updated all the
# germination data and manually imported some rows

#Need to merge my most up to date file, onlygerminatedsubplotsedited, back with checklist data
# to get all subplots (things with 0 germination as well)
checklist <- read_csv("Data/Manually_edited/Checklist_germinated_plants.csv")
#Need to drop the blanks before I merge with onlygerminated (where NAs become 0s)
#This is also removing a lot of PLDEs (which I never planted anyway!)
#Have to sort out two  weird ones popping up - 1 B POLE E3 and 4 A PLDE T@
#test <- checklist %>% filter(is.na(Germinate))
checklist <- checklist %>% drop_na(Germinate)
onlygerminated <- read_csv("Data/Manually_edited/onlygerminatedsubplotsedited.csv")
#Two duplicate rows with the same germination values, so reducing to one row each
onlygerminated <- onlygerminated %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(row_number() == 1)
onlygerminatedselect <- onlygerminated %>% select(Site, Plot, Species, C_E_or_T, Rep, Germinate, Notes)
subplotschecklist <- checklist %>% select(Site, Plot, Species, C_E_or_T, Rep)
germinationdata <- left_join(subplotschecklist, onlygerminatedselect)
germinationdata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% filter(!row_number() == 1)