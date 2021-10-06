########## Abigail's 2017  data
## Some plots watered by Trace, camp showers
#SLA outliers removed
# Fecundity data
#l1 is, I think, the length from the crown to the tip of the longest leaf?
# So, longest length and then l2 is the length perpendicular to that

library(tidyverse)
library(ggplot2)
library(ggpubr)

water2017data <-read_csv("MasterAllData16.7.19outliersRM.csv")

#Need to get an average per species for wet and dry
# Note that in below code we get NA for meanseeds if there is only NA seed data
# for that species x environment and we get NA for SD in that situation and where
# there are less than 3 data points

meanseeds <- water2017data %>% group_by(focal.sp, environment) %>% 
  summarise(meanseeds = mean(seeds, na.rm=TRUE), sdseeds = sd(seeds, na.rm=TRUE))

#Also want a count of how many data points we have for each species x env in df
#Not properly filtering NAs atm, this below count includes them

water2017data %>% count(focal.sp, environment)


#Can't figure it out, so just going to remove NAs from the tibble
#Starting from the top! There are NAs and 0s but 0s are important!
#There is one random row with NA for focal species so dropping that

meanseeds <- water2017data %>% select(focal.sp, environment, seeds) %>%   
  #mutate(seeds= na_if(seeds, 0)) %>% 
  drop_na(seeds, focal.sp) %>% 
  group_by(focal.sp, environment) %>% 
  summarise(meanseeds = mean(seeds, na.rm=TRUE), sdseeds = sd(seeds, na.rm=TRUE))


#Count without NAs

water2017dataNArm <- water2017data %>% 
  #mutate(seeds = na_if(seeds, 0)) %>%
                    drop_na(seeds) %>% drop_na(focal.sp) 

water2017dataNArmcount <- water2017dataNArm %>% count(focal.sp, environment)

#Adding count column to the other dataset

meanseeds <- meanseeds %>% left_join(water2017dataNArmcount)

#Means continued
#SLA
meansla <- water2017data %>% select(focal.sp, environment, sla) %>%   
  mutate(sla= na_if(sla, 0)) %>% 
  drop_na(sla, focal.sp) %>% 
  group_by(focal.sp, environment) %>% 
  summarise(meansla = mean(sla, na.rm=TRUE), sdsla = sd(sla, na.rm=TRUE))
water2017dataNArm2 <- water2017data %>% 
  mutate(sla = na_if(sla, 0)) %>%
  drop_na(sla) %>% drop_na(focal.sp) 
water2017dataNArmcount2 <- water2017dataNArm2 %>% count(focal.sp, environment)
meansla <- meansla %>% left_join(water2017dataNArmcount2)

#Height
#Forcing height to be numerical

meanheight <- water2017data %>% select(focal.sp, environment, height) %>% 
  mutate(height = na_if(height, 0)) %>%
  mutate(height = na_if(height, '#N/A')) %>%
  drop_na(height, focal.sp) %>% 
  mutate(height = as.numeric(height)) %>%
  group_by(focal.sp, environment) %>% 
  summarise(meanheight = mean(height, na.rm=TRUE), sdheight = sd(height, na.rm=TRUE))
water2017dataNArm3 <- water2017data %>% 
  mutate(height = na_if(height, 0)) %>%
  drop_na(height) %>% drop_na(focal.sp) 
water2017dataNArmcount3 <- water2017dataNArm2 %>% count(focal.sp, environment)
meanheight <- meanheight %>% left_join(water2017dataNArmcount3)

#density
meandensity <- water2017data %>% select(focal.sp, environment, density) %>%   
  mutate(density = na_if(density, 0)) %>% 
  drop_na(density, focal.sp) %>% 
  group_by(focal.sp, environment) %>% 
  summarise(meandensity = mean(density, na.rm=TRUE), sddensity = sd(density, na.rm=TRUE))
water2017dataNArm4 <- water2017data %>% 
  mutate(density = na_if(density, 0)) %>%
  drop_na(density) %>% drop_na(focal.sp) 
water2017dataNArmcount4 <- water2017dataNArm2 %>% count(focal.sp, environment)
meandensity <- meandensity %>% left_join(water2017dataNArmcount4)


#All species box plot
ggplot(water2017dataNArm, aes(x=focal.sp, y=seeds)) + 
  geom_boxplot(aes(fill = environment))+
  #ylim(0,50) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Filter to only these species that have at least 3 fecundity measurements for
# each of wet and dry
waterdatasub <- water2017dataNArm %>% 
  filter(focal.sp %in% c("VERO","ARCA", "POCA", "PLDE","MEDI", "AUEL", "MOMO",
                        "GOBE","HYGL", "POLE", "PEAI", "PTGA"))

waterdatasub <- waterdatasub %>% mutate(logoneseeds = log(seeds+1))

#Plotting subsetted data
ggplot(waterdatasub, aes(x=focal.sp, y=logoneseeds)) + 
  geom_boxplot(aes(fill = environment))+
  geom_jitter() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Zooming in, note that cartesian coordinate does not clip my data
#Note for stat_compare_means: change label to ..p.signif.. to get * and nas
#..p.format.. shows p-values

ggplot(waterdatasub, aes(x=focal.sp, y=logoneseeds)) + 
  geom_boxplot(aes(fill = environment))+
  geom_jitter(aes(colour = environment), position=position_jitter(0.05)) +
  theme_bw() +
 coord_cartesian(ylim = c(0, 7.5))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
stat_compare_means(method = "t.test", (aes(group = environment, 
                                           label = ..p.format..)),
                                    label.y = 0.5, size = 2.5)+
  ylab('log(number of seeds+1)')+
  labs(title='Species fecundity in watered or ambient unweeded neighbourhoods') +
  facet_wrap(~ focal.sp, scales="free")

ggplot(waterdatasub, aes(x=focal.sp, y=seeds)) + 
  geom_boxplot(aes(fill = environment))+
  theme_bw() +
  coord_cartesian(ylim = c(0, 150))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(method = "t.test", (aes(group = environment, 
                                             label = ..p.format..)),
                     label.y = -3, size = 3)

#Trace said that boxplots show you the medians and interquartile ranges of
# the data, so displaying it with a t-test (testing means and variance)
# could be misleading. So, trying mean and SD instead
# for SEEDS data

ggplot(waterdatasub, aes(x=focal.sp, y=log(seeds+1), fill = environment)) + 
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.01), stat="summary", 
             fun.y="mean",size=1.5, show.legend=FALSE) +
  geom_errorbar(position=position_jitterdodge(0.01), stat="summary", 
                fun.data="mean_sd",width=0,size=0.7) +
  facet_wrap(~ focal.sp, scales="free")+
  stat_compare_means(method = "t.test", (aes(group = environment, 
                                             label = ..p.format..)),
                     label.x.npc = 'middle', label.y.npc = 'top', size = 3) +
  labs(y = "log(number of seeds+1)", x = "", 
       title = "Seed number in watered or ambient unweeded neighbourhoods",
       subtitle = "Abigail's 2017 Data - Subset of species")

#Stats time! 

for(x in unique(waterdatasub$focal.sp)) {
stats <- t.test(waterdatasub$seeds[which(waterdatasub$focal.sp== x & 
                                  waterdatasub$environment == "Water")], 
       waterdatasub$seeds[which(waterdatasub$focal.sp== x & 
                                      waterdatasub$environment == "Dry")])
print(x)
print(stats)
}

####### Now looking at SLA

water2017dataNArmsla <- water2017data %>% drop_na(sla) %>% drop_na(focal.sp) 

#All species box plot
ggplot(water2017dataNArmsla, aes(x=focal.sp, y=sla)) + 
  geom_boxplot(aes(fill = environment))+
  #ylim(0,50) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Note that below is only filtering out wet or dry for some species, rather than
# entire species (ideal)
#Not sure exactly what this code does anymore
water2017dataNArmslathree <- water2017dataNArmsla %>% group_by(focal.sp, environment) %>% 
  filter(n() >= 3)

#Filter to my species of interest, not including GOBE which 
#doesn't have enough wet data

water2017datasubNArmsla <- water2017dataNArmsla %>% 
  filter(focal.sp %in% c("VERO","ARCA", "POCA", "PLDE", "HYGL",
                         "POLE", "PEAI"))

#Plot for my species, with stats

ggplot(water2017datasubNArmsla, aes(x=focal.sp, y=log(sla), fill = environment)) + 
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.01), stat="summary", 
             fun.y="mean",size=1.5, show.legend=FALSE) +
  geom_errorbar(position=position_jitterdodge(0.01), stat="summary", 
                fun.data="mean_sd",width=0,size=0.7) +
  facet_wrap(~ focal.sp, scales="free")+
  stat_compare_means(method = "t.test", (aes(group = environment, 
                                             label = ..p.format..)),
                     label.x.npc = 'middle', label.y.npc = 'top', size = 3) +
  labs(y = "log(SLA)", x = "", 
       title = "SLA in watered or ambient unweeded neighbourhoods",
       subtitle = "Abigail's 2017 Data - Subset of species")

#Same, for SLA

ggplot(waterdatasub, aes(x=focal.sp, y=log(sla), fill = environment)) + 
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.01), stat="summary", 
             fun.y="mean",size=1.5, show.legend=FALSE) +
  geom_errorbar(position=position_jitterdodge(0.01), stat="summary", 
                fun.data="mean_sd",width=0,size=0.7) +
  facet_wrap(~ focal.sp, scales="free")+
  stat_compare_means(method = "t.test", (aes(group = environment, 
                                             label = ..p.format..)),
                     label.x.npc = 'middle', label.y.npc = 'top', size = 3) +
  labs(y = "log(SLA)", x = "", 
       title = "SLA in watered or ambient unweeded neighbourhoods",
       subtitle = "Abigail's 2017 Data - Subset of species")

#Now for height

water2017NArmheight <- water2017data %>% mutate(height = na_if(height, '#N/A')) %>%
  drop_na(height) %>% drop_na(focal.sp)
  #Note that below gives warning message, NAs introduced by coercion
  # Not sure where or if it's a problem
water2017NArmheight$height = as.numeric(water2017NArmheight$height)

#All species box plot
ggplot(water2017NArmheight, aes(x=focal.sp, y=height)) + 
  geom_boxplot(aes(fill = environment))+
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Confused about my filtering method way above, for sla. Not specifically n is sla
# anywhere

#Note that TRCY and GOBE don't have enough wet height data

water2017NArmheightsub <- water2017NArmheight %>% 
  filter(focal.sp %in% c("VERO","ARCA", "POCA", "PLDE", "HYGL",
                         "POLE", "PEAI"))

ggplot(water2017NArmheightsub, aes(x=focal.sp, y=height, fill = environment)) + 
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.01), stat="summary", 
             fun.y="mean",size=1.5, show.legend=FALSE) +
  geom_errorbar(position=position_jitterdodge(0.01), stat="summary", 
                fun.data="mean_sd",width=0,size=0.7) +
  facet_wrap(~ focal.sp, scales="free")+
  stat_compare_means(method = "t.test", (aes(group = environment, 
                                             label = ..p.format..)),
                     label.x.npc = 'middle', label.y.npc = 'top', size = 3) +
  labs(y = "Height (cm)", x = "", 
       title = "Height in watered or ambient unweeded neighbourhoods",
       subtitle = "Abigail's 2017 Data - Subset of species")
#I think HYGL and VERO definitely have outliers. Nonetheless, significant.
#Is t-test assuming same variance?

#Now for density


water2017NArmdensity <- water2017data %>% drop_na(density) %>% drop_na(focal.sp)

#All species box plot
ggplot(water2017NArmdensity, aes(x=focal.sp, y=density)) + 
  geom_boxplot(aes(fill = environment))+
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#TRCY and GOBE don't have enough data
water2017NArmdensitysub <- water2017NArmdensity %>% 
  filter(focal.sp %in% c("VERO","ARCA", "POCA", "PLDE", "HYGL",
                         "POLE", "PEAI"))

ggplot(water2017NArmdensitysub, aes(x=focal.sp, y=log(density+1), fill = environment)) + 
  geom_jitter(aes(colour = environment), size=1.6, alpha = 0.3, 
              position=position_jitterdodge(0.1)) +
  theme_bw() +
  geom_point(position=position_jitterdodge(0.01), stat="summary", 
             fun.y="mean",size=1.5, show.legend=FALSE) +
  geom_errorbar(position=position_jitterdodge(0.01), stat="summary", 
                fun.data="mean_sd",width=0,size=0.7) +
  facet_wrap(~ focal.sp, scales="free")+
  stat_compare_means(method = "t.test", (aes(group = environment, 
                                             label = ..p.format..)),
                     label.x.npc = 'middle', label.y.npc = 'top', size = 3) +
  labs(y = "log(Density (?)+1)", x = "", 
       title = "Density in watered or ambient unweeded neighbourhoods",
       subtitle = "Abigail's 2017 Data - Subset of species")


