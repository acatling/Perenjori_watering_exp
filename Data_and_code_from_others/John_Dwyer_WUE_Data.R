##### John Dwyer's Water Use Efficiency C13 isotope discrimination data
### Western Australia annuals

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(mosaic)
 
wuedata <-read_csv("Data/YGW_carbon_isotopes_March_2019_JD.csv")

notes_to_remove <-c('duplicate', 'not c3')
#Consider also removing ones ith 'bulked' in notes

wuedatanodups <- wuedata %>% filter(!notes%in%notes_to_remove)

mean_wue <- wuedatanodups %>%
  group_by(Code) %>%
  summarise(mean_wue = mean(delta, na.rm= TRUE),
            sd_wue = sd(delta, na.rm=TRUE))

  ggplot(wuedatanodups, aes(x=Code, y = delta)) + 
  geom_jitter(position=position_jitter(0.1), colour = "dodgerblue") + 
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

#Filter out species with less than three data points
wuedatanodupsthree <- wuedatanodups %>% group_by(Code) %>% filter(n() >= 3)
  
#Plot those with >=3 data points
ggplot(wuedatanodupsthree, aes(x=Code, y = delta)) + 
  geom_jitter(position=position_jitter(0.1), colour = "dodgerblue") + 
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab('Carbon isotope discrimination (delta, ppm)')
  
  #Only looking at the subset that I'm interested in
wuedatanodupssubthree <- wuedatanodupsthree %>% 
  filter(Code %in% c("VERO","ARCA", "POCA", "PLDE","DAGL", "HYPO", "TRCY",
                         "GOBE","HYGL", "POLE", "PEAI", "TROR", "LARO", "SCCA"))

#Plot my subset which have >3 data points
ggplot(wuedatanodupssubthree, aes(x=Code, y = delta)) + 
  geom_jitter(position=position_jitter(0.1), colour = "dodgerblue") + 
  geom_boxplot(fill = NA, outlier.shape = NA)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab('Carbon isotope discrimination (delta, ppm)')+
  labs(title = 'Water use efficiency data for subset of species', 
       subtitle = 'One-way ANOVA with Tukey HSD post-hoc suggests POLE is significantly different')

# Stats time
#One way ANOVA
#Only looking at data with more than three points

model1 <- aov(delta ~ Code, data = wuedatanodupssubthree)
summary(model1)
aov(model1)
confint(model1)
#Data appears to be more or less? normally distributed
hist(residuals(model1))
modelwue = data.frame(Fitted = fitted(model1), Residuals = resid(model1), 
                      Treatment = wuedatanodupssubthree$Code)
#Looks like there are unequal variances between species.
ggplot(modelwue, aes(Fitted, Residuals, colour = Treatment)) + geom_point()

#Bartlett’s test and Levene’s test can be used to check the homoscedasticity 
#of groups from a one-way anova.  A significant result for these tests 
#(p < 0.05) suggests that groups are heteroscedastic
#https://rcompanion.org/rcompanion/d_05.html
#Heteroscedastic means that have different variances.
bartlett.test(delta ~ Code, data = wuedatanodupssubthree)
#p-value is 0.08, so I think we are okay treating them as equal variances.

#Performing Tukey Honest Significant Differences post-hoc test
TukeyHSD(model1)
#POLE is significantly different to ARCA, GOBE, LARO and POCA
              

              