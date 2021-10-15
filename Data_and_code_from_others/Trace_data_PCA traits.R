#PCA
# Ambient: Species, SLA, WUE, fecundity, height, density - dry
# Wet: Species, SLA (amb-wet), fecundity/height/density (amb-wet)
#Need long data
#Will also add Trace's instrinsic fecundity as a trait, for wet PCA
#Following this guide - great! 
#https://www.datacamp.com/community/tutorials/pca-analysis-r

library(tidyverse)
library(devtools)
library(ggbiplot)

#Note that I can't display labels on ggbiplot with NAs (first line below)
pcadatawithNAs <- read_csv("Data/Perenjori_PCA.csv")
pcadata <- read_csv("Data/Perenjori_PCA_dry_NArm.csv")

#Just ambient (dry) data
pcadatadry <- pcadata %>% select(Code, SLA_dry, WUE, Seeds_dry, 
                                 Height_dry, Density_dry)
  
#Dry PCA
pcadry <- prcomp(na.omit(pcadatadry[,c(2:6)]), center = TRUE,scale. = TRUE)
summary(pcadry)
str(pcadry)

install_github("vqv/ggbiplot")
#Note that I have to say none (3) in response to above line

ggbiplot(pcadry)
ggbiplot(pcadry, labels=row.names(pcadatadry))+
  theme_classic()

#Base R
biplot(pcadry, choices = 1:2, scale = 1, pc.biplot=FALSE)

#I think it's just restricting to the 7 species without NAs, but labels
# Aren't working so I'm just going to try that myself
#Nevermind, it was working! Just labelled with numbers for some reason.

pcadatadry
#Match numbers to rows in table (very first column of 1-7)
#So 6 and 7, GOBE and POCA are most similar? 









