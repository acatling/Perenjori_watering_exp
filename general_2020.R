#General info Perenjori Watering Experiment 2020
## Alexandra Catling March 2021

library(tidyverse)

watering <- read_csv("Data/Watering_amounts_2020.csv")
#Note that where watered wasn't added, I've put '0' under water_can_mm (to avoid NA issues)

#Adding a column to calculate the volume of water added, currently just the height (mm)
# in the watering can.
#Formula for volume in a cylinder is volume = pi * radius^2 * height
#Diameter of watering can from Bunnings website is 225 mm, so radius is 112.5
#This will give cubic millimeters with my units. 1mm^3 = 1e-6 L.
#So need to divide it all by 1e+6 to get L.

wateringdata <- watering %>% mutate(volume = (pi*112.5^2*Water_can_mm/2)/1000000)

wateringbysite <- wateringdata %>% group_by(Site) %>% summarise(total_volume = sum(volume))

str(wateringbysite)
wateringbysite$Site <- as.factor(wateringbysite$Site)

ggplot(wateringbysite, aes(x = Site, y = total_volume))+
  geom_point()+
  theme_classic()

## Want to see if watering treatment affected survival
#Wet plots: 1B, 2C, 3B, 4C, 5C, 6B, 7B, 8B
#Dry plots: 1C, 2A, 3A, 4B, 5B, 6A, 7A, 8A

treatments <- read_csv("Data/treatments_meta.csv")
survivaltreatment <- merge(treatments, survivalplot)

#Reordering in order of dry, control, wet.

level = c("Dry", "Control", "Wet")
survivaltreatment$Treatment <- factor(survivaltreatment$Treatment, level = c("Dry", "Control", "Wet"))

ggplot(survivaltreatment, aes(x = Treatment, y = n))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  ylab("Number of focals that survived to set seed")

ggplot(survivaltreatment, aes (x = Cover, y = n))+
  geom_boxplot()+
  geom_point()+
  theme_classic()+
  ylab("Number of focals that survived to set seed")

ggplot(survivaltreatment, aes (x = Treatment, y = n, colour = Cover))+
  geom_boxplot()+
  geom_jitter(position = position_jitter(width = 0.1, height = NULL))+
  theme_classic()+
  ylab("Number of focals that survived to set seed")

model1 <- lm(n ~ Treatment, survivaltreatment)
summary(model1)

model2 <- lm(n ~ Treatment + Cover, survivaltreatment)
summary(model2)

model3 <- lm(n ~ Treatment + Cover + Treatment*Cover, survivaltreatment)
summary(model3)

AIC(model1, model2, model3)

#BUT need to factor in that some of these died before I imposed the watering treatment



