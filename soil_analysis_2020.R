###################### Analysis of Soil from 2020 Perenjori Experiment
########## Alexandra Catling
### March 2021

library(tidyverse) 
library(lme4)
library(ggplot2)
library(emmeans)

#Importing data
soildata <- read_csv("Data/soil_analysis_2020.csv")
head(soildata)
str(soildata)
soildata$Site <- as.factor(soildata$Site)
str(soildata$Site)

###################Visualising and transforming distributions of data

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=seq(min(x, na.rm=T), max(x, na.rm=T), length.out=7))
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y, na.rm=T)
  rect(breaks[-nB], 0, breaks[-1], y,...)
}

with(soildata, pairs(Site ~ pH + `NH4-N` + `NO3-N` + Al + B + Ca + Cu, diag.panel = panel.hist))
#pH normally distributed
#NH4-N has one huge outlier in Sample 1, logging helps
#NO3-N also has a big outlier (not as big) in Sample 1, logging helps
#Al left-skewed, logging maaaybe helps
#B has an outlier, logging doesn't really help
#Ca left-skewed, logging helps?
#Cu maybe left-skewed, logging helps
with(soildata, pairs(Site ~ pH + log(`NH4-N`) + log(`NO3-N`) + log(Al) + log(B) + log(Ca) + log(Cu), diag.panel = panel.hist))

with(soildata, pairs(Site ~ Fe + K + Mg + Mn + Na + P + S + Zn, diag.panel = panel.hist))
#Fe left-skewed, logging  kinda helps?
#K left-skewed, logging helps
#Mg more or less consistent
#Mn left-skewed with an outlier, logging helps
#Na left-skewed, logging kinda helps?
#P bit left-skewed, logging helps
#S normally distributed
#Zn left-skewed, logging doesn't really help
with(soildata, pairs(Site ~ log(Fe) + log(K) + Mg + log(Mn) + log(Na) + log(P) + S + log(Zn), diag.panel = panel.hist))

#Generating new columns for transformed variables
#Combining ammonium and nitrate to approximate plant-available nitrogen (N)

soildata <- soildata %>% mutate(log_NH4N = log(`NH4-N`), log_NO3N = log(`NO3-N`), log_Al = log(Al), log_Ca = log(Ca), log_Cu = log(Cu),
                                log_Fe = log(Fe), log_K = log(K), log_Mn = log(Mn), log_Na = log(Na), log_P = log(P), N = `NH4-N` + `NO3-N`, log_N = log(N))
glimpse(soildata)

#Making a list of variables (transformed) that I want to analyse (two different approaches here)
soildatavars <- soildata %>% select(Sample, Site, Plot, pH, log_NH4N, log_NO3N, log_Al, B, log_Ca, log_Cu, log_Fe, log_K, Mg,
                                    log_Mn, log_Na, log_P, S, Zn)
#Margie and John suggested CEC (soil cation exchange capacity), nitrate, pH, P, K.
soildataselect <- soildata %>% select(Sample, Site, Plot, pH, log_N, log_P, log_K)

### Need to make data tidy before I can plot all variables against Site
soildatatidy <- soildataselect %>% pivot_longer(cols = pH:log_K, names_to = "Variable", values_to = "Values")
  
ggplot(soildatatidy, aes(x = Site, y = Values))+
  geom_boxplot()+
  geom_point(colour = "dodgerblue")+
  theme_classic()+
  facet_wrap(~Variable, scales = "free")

  
########### pH example
ggplot(soildata, aes(x = Site, y = pH))+
  geom_boxplot()+
  geom_point(colour="dodgerblue")+
  theme_classic()

# Another way of representing it
# ggplot(soildata, aes(x=Site, y=pH)) + 
#   geom_point(colour="dodgerblue") +
#   geom_point(stat="summary", fun=mean, size=2.5) +
#   geom_errorbar(stat="summary", fun.data="mean_se",width=0,size=0.8)+
#   theme_bw()

#Statistical analysis
modelph1 <- aov(pH ~ Site, soildata)
summary(modelph1)
TukeyHSD(modelph1)

# Site 3 is significantly different from 1, 2, 4, 5 and 6

########### NH4-N
ggplot(soildata, aes(x = Site, y = log(`NH4-N`)))+
  geom_boxplot()+
  geom_point(colour="dodgerblue")+
  theme_classic()

par(mfrow=c(4, 4))

#Does this map to survival?

survivalsoil <- merge(soildata, survivalplot)

#pH
ggplot(survivalsoil, aes(y = n, x = pH))+
  geom_point(aes(colour = Site))+
  geom_smooth(method = "lm")+
  theme_classic()

survivalsoiltidy <- merge(soildatatidy, survivalplot)

#All
ggplot(survivalsoiltidy, aes(x = Values, y = n))+
  geom_point()+
  theme_classic()+
  geom_smooth(method = "lm")+
  facet_wrap(~Variable, scales = "free")

model1 <- lm(n ~ pH, survivalsoil)
summary(model1)

### Coming back to this May 2021
# Merging with seed count data

seedsoil <- merge(seeddata, soildata)
seedsoiltidy <- merge(seeddata, soildatatidy)

ggplot(seedsoil, aes(x = pH, y = log1_viable_seeds))+
  geom_jitter(alpha = 0.4)+
  geom_smooth(method="lm")+
  theme_classic()

ggplot(seedsoiltidy, aes(x = Values, y = log1_viable_seeds))+
  geom_point(alpha= 0.3)+
  theme_classic()+
  geom_smooth(method = "lm")+
  facet_wrap(~Variable, scales = "free")

seedsoilmodel1 <- glmer.nb(No_viable_seeds ~ pH + Na + S + (1|Plot), data = seedsoil)
summary(seedsoilmodel1)

######################### PCA #######################
### Coming back to this July 2021

#dataall from full_model_WA script
soildataall <- left_join(dataall, soildataselect)

abioticpcadata <- soildataall %>% select(plotid, pH, log_N, cc_percentage, log_P, log_K)

abioticpcadata <- as.data.frame(abioticpcadata[,-1])
rownames(abioticpcadata) <- abioticpcadata$plotid

soil_pca <- princomp(abioticpcadata, cor = TRUE)
biplot(soil_pca)
summary(soil_pca)
loadings(soil_pca)
#From my interpretation, PCA1 accounts for 57% of total variance, and is mostly described
# by N, canopy cover, P and K. PCA2 accounts for a further 25% (together a total of 82%)
# and is mostly described by pH (and K)

#Checking correlations (sure there's a better way to do this)
ggplot(soildataselect, aes(x = log_K, y = log_N))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
ggplot(soildataselect, aes(x = log_K, y = log_P))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
ggplot(soildataselect, aes(x = log_K, y = pH))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
ggplot(soildataselect, aes(x = log_N, y = log_P))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
ggplot(soildataselect, aes(x = log_N, y = pH))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
ggplot(soildataselect, aes(x = log_P, y = pH))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()