######### Trace's 2017 data
###### Intrinsic plots (competitors thinned) watered and unwatered

library(tidyverse)
library(ggplot2)
library(ggpubr)

intrinsic_seedno<-read_csv("Data/2017.Seed.Count.Intrinsic.Data.7.6.2018.TM.csv")

#All species viewing plotted data points (not averaged)
#Note that ylim removes some points
ggplot(intrinsic_seedno, aes(x = Species, y = Num.Seed, group = Wet.Dry)) + 
  geom_jitter(aes(colour = Wet.Dry))+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+
  ylim(0,500)

#All species box plot
ggplot(data=intrinsic_seedno, aes(x=Species, y=Num.Seed)) + 
  geom_boxplot(aes(fill = Wet.Dry))+
  ylim(0,500) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Counting the number of observations (plants) for seed amount per species
#sort = TRUE sorts highest to lowest, FALSE is useful to compare number 
# wet and dry for each species.
intrinsic_seedno %>% 
  count(Species, Wet.Dry, sort = TRUE) %>% arrange(Species) 

#Filter to only these species that have at least 3 measurements for 
# each of wet and dry
intrinsic_seedno_sub <- intrinsic_seedno %>% filter(Species %in% c("VERO", 
                                "ARCA", "POCA", "PLDE", "TRCY", "WAAC", "MEDI",
                                "EROD", "GITE","GOBE", "GOPU", "HYPO", "DAGL",
                                "HYGL", "POLE"))
intrinsic_seedno_sub <- intrinsic_seedno_sub %>% mutate(logoneseed = log(Num.Seed+1))

#Plotting subsetted data
ggplot(data=intrinsic_seedno_sub, aes(x=Species, y=logoneseed)) + 
  geom_boxplot(aes(fill = Wet.Dry))+
  geom_jitter(aes(colour = Wet.Dry), position=position_jitter(0.2))+
  #coord_cartesian(ylim=c(0,1500)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(method = "t.test", (aes(group = Wet.Dry, label = ..p.format..)),
                     size = 2.5, label.y = -1)+
  labs(y='log(number of seeds +1', title = 'Intrinsic fecundity of species in watered and ambient plots')

