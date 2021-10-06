##################### Analysing germination data

library(tidyverse)
library(sjPlot)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)

source("data_preparation.R")
#germinationdata is the tibble I want to use!
str(germinationdata)
germinationdata$Site <- as.factor(germinationdata$Site)

#Calculating germination rate.
## Accounting for February germination. For now, total number germinated / total number seeds
germinationdata <- germinationdata %>% group_by(Site, Plot, Species, C_E_or_T, Rep) %>% 
  mutate(Germination_rate = (February_germination + Number_germinated)/Number_seeds_sown)

#Did any germinate in February and then not later in the year?
test <- germinationdata %>% filter(February_germination > 0 & Number_germinated == 0)
#Yes, 16 of them! Wow. Maybe I remove these plots? For now, taking them out.
germinationdata <- germinationdata %>% filter(!(February_germination > 0 & Number_germinated == 0))

## Some plots
tallysitegerm <- germinationdata %>% filter(!Germination_rate == "0") %>% group_by(Site, Species) %>% tally()
ggplot(tallysitegerm, aes(x = Site, y = n), group = Species)+
  geom_point(aes(colour = Site))+
  theme_classic()+
  ylab("Number of germinated subplots")+
  facet_wrap(~Species)+
  ggtitle("Germination tally by site")

tallyplotgerm <- germinationdata %>% filter(!Germination_rate == "0") %>% group_by(Site, Plot, Species)  %>% tally()
ggplot(tallyplotgerm, aes(x = Site, y = n), group = Species)+
  geom_jitter(aes(colour = Site))+
  theme_classic()+
  ylab("Number of germinated subplots")+
  facet_wrap(~Species)+
  ggtitle("Germination tally by plot")

### Making germination binomial and running models
#Make a new column for Germinated (1 or 0)
germinationdata <- germinationdata %>% mutate(Germinated01 = case_when(Germination_rate > "0" ~ "1",
                                                              Germination_rate == "0" ~ "0"))
#Need to import Cover info
canopylitterdata <- read_csv("Data/canopy_litter_2020.csv")
canopydata <- canopylitterdata  %>% rowwise() %>%
  mutate(meancc = mean(c(`Measurement 1`, `Measurement 2`, `Measurement 3`, `Measurement 4`, `Measurement 5`)),
         cc_percentage = meancc/24*100)
canopydatatrim <- canopydata %>% select(Site, Plot, cc_percentage)
germinationdata <- merge(germinationdata, canopydatatrim)

####Creating a dataset for each species
for (i in 1:length(specieslist)){
  nam <- paste0("germ", specieslist[i])
  assign(nam, germinationdata %>% filter(Species == specieslist[i]))
}

### Models
str(germinationdata)
germinationdata$Plot <- as.factor(germinationdata$Plot)
germinationdata$Germinated01 <- as.factor(germinationdata$Germinated01)

## Models by species
specieslist <- c("ARCA", "HYGL", "LARO", "PEAI", "PLDE", "POLE", "TRCY", "TROR", "VERO")
for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
    glmer(Germinated01 ~ cc_percentage + (1|Site/Plot), family = binomial, 
          data = filter(germinationdata, Species == specieslist[i]))))
}
#HYGL, PEAI, POLE, TRCY, TROR, VERO. Germination significantly negatively influenced
# by canopy closure for all of these
## Naming models so that I can plot coefs from them
for (i in 1:length(specieslist)){
  nam <- paste0("germmodel", specieslist[i])
  assign(nam, glmer(Germinated01 ~ log(cc_percentage+1) + (1|Site/Plot), family = binomial, 
                    data = filter(germinationdata, Species == specieslist[i])))
}

summary(germmodelARCA)
tab_model(germmodelARCA, germmodelHYGL, germmodelLARO, germmodelPLDE, germmodelPOLE, germmodelPEAI, germmodelTRCY, germmodelTROR, germmodelVERO, transform = NULL)

##Plotting coefficients
allmodels <- list()
allmodels[[1]] <- germmodelARCA
allmodels[[2]] <- germmodelHYGL
allmodels[[3]] <- germmodelLARO
allmodels[[4]] <- germmodelPLDE
allmodels[[5]] <- germmodelPOLE
allmodels[[6]] <- germmodelPEAI
allmodels[[7]] <- germmodelTRCY
allmodels[[8]] <- germmodelTROR
allmodels[[9]] <- germmodelVERO

plot_models(allmodels, transform = NULL, vline.color = "grey")+
  ylab("Estimate")+
  theme_classic()+
  scale_colour_discrete(labels = c("VERO", "TROR", "TRCY", "POLE", "PLDE", "PEAI", "LARO", "HYGL", "ARCA"))
#Can't make x axis limits any smaller using sJplot (can use a work around with plot_model in ggplot if I want to)

###### Need to calculate germination rates per subplot (number germinated/total number of seeds)
# Model by germination rate

#First checking distribution of germination rate
ggplot(germinationdata, aes(x = log(Germination_rate)+1))+
  geom_histogram()+
  theme_classic()+
  facet_wrap(~Species)

for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    summary(
      lmer(Germination_rate ~ log(cc_percentage+1) + (1|Site/Plot),
           data = filter(germinationdata, Species == specieslist[i]))))
}
#All have a significant negative coefficient for canopy closure except ARCA

for (i in 1:length(specieslist)){
  nam <- paste0("germmodelc", specieslist[i])
  assign(nam, lmer(Germination_rate ~ log(cc_percentage+1) + (1|Site/Plot),  
                    data = filter(germinationdata, Species == specieslist[i])))
}
####Calculating germination rates
#Import Sun vs. Shade info
treatmentdata <- read_csv("Data/treatments_meta.csv")
germinationdata <- merge(germinationdata, treatments)

germinationcounts <- germinationdata %>% group_by(Species, Cover, Germinated01) %>% tally()
#I already have germination rates calculated... so I just need to average these?
germinationsp <- germinationdata %>% group_by(Species) %>% summarise(mean_germ_rates = mean(Germination_rate),
                                                            sd_germ_rates = sd(Germination_rate))
germinationcoversp <- germinationdata %>% group_by(Species, Cover) %>% summarise(mean_germ_rates = mean(Germination_rate),
                                                                          sd_germ_rates = sd(Germination_rate))
#Comparing these species-level values to another way of measuring it
germcount <- germinationdata %>% group_by(Species, Germinated01) %>% tally()
germsp <- germcount %>%
  group_by(Species) %>%
  mutate(number_germinated = sum(n[Germinated01 == "1"]),
         total_number = sum(n[Germinated01 == "1"])+sum(n[Germinated01 == "0"]),
         germination_rate = number_germinated/total_number) %>%
  filter(row_number() == 1)
#Selecting only relevant columns
germsp <- germsp %>% select(Species, germination_rate)
################ THIS GIVES VERYYYY DIFFERENT VALUES, TRUST THIS MORE.
#And split between sun/shade
germshadecount <- germinationdata %>% group_by(Species, Cover, Germinated01) %>% tally()
germshadesp <- germshadecount %>%
  group_by(Species, Cover) %>%
  mutate(number_germinated = sum(n[Germinated01 == "1"]),
         total_number = sum(n[Germinated01 == "1"])+sum(n[Germinated01 == "0"]),
         germination_rate = number_germinated/total_number) %>%
  filter(row_number() == 1)
#Selecting only relevant columns
germshadesp <- germshadesp %>% select(Species, Cover, germination_rate)


ggplot(germinationdata, aes(x=Species, y=Germination_rate)) + 
  geom_jitter(colour="forestgreen", size=1.6, alpha = 0.2,width = 0.2)+
  geom_point(stat="summary", size=2) +
  geom_errorbar(stat="summary", width=0,size=0.8)+
  theme_classic()

#Defaults to mean and se
#ggplot(germinationdata, aes(x=Cover, y=Germination_rate)) + 
#  geom_jitter(colour="forestgreen", size=1.6, alpha = 0.2,width = 0.2)+
#  geom_point(stat="summary", size=1.5) +
#  geom_errorbar(stat="summary", width=0,size=0.8)+
#  theme_classic()+
#  facet_wrap(~Species)

ggplot(germinationdata, aes(x = Cover, y = Germination_rate))+
  geom_boxplot()+
  geom_jitter(colour = "steelblue3", alpha = 0.2, width = 0.2)+
  theme_classic()+
  facet_wrap(~Species)

#Series of t-tests

for (i in 1:length(specieslist)){
  print(specieslist[i])
  print(
    t.test(Germination_rate ~ Cover, data = filter(germinationdata, Species == specieslist[i])))
}
#All significantly different except ARCA. All germinated more in the sun (except ARCA).

ggplot(germinationdata, aes(x = log(cc_percentage+1), y = Germination_rate))+
  geom_point(alpha = 0.4)+
  theme_classic()+
  geom_smooth(method="lm")+
  facet_wrap(~Species)

r.squaredGLMM(germmodelcARCA)
r.squaredGLMM(germmodelcHYGL)
r.squaredGLMM(germmodelcLARO)
r.squaredGLMM(germmodelcPEAI)
r.squaredGLMM(germmodelcPLDE)
r.squaredGLMM(germmodelcPOLE)
r.squaredGLMM(germmodelcTRCY)
r.squaredGLMM(germmodelcTROR)
r.squaredGLMM(germmodelcVERO)

## Trying Isaac's plots
source("Isaac_functions.R")

###create a list of the dataframes which only have C data
species.c.list<-list(germARCA, germHYGL, germLARO, germPEAI, germPLDE, germPOLE, germTRCY, germTROR, germVERO)

###create a list of names for figure headings
species.name.list<-c("Arctotheca calendula","Hyalosperma glutinosum","Lawrencella rosea","Pentameris airoides","Plantago debilis","Podolepis lessonii","Trachymene cyanopetala","Trachymene ornata","Velleia rosea")

dev.off()
pdf("Output/testing.pdf", width=21, height=21)
par(mfrow=c(3,3))
par(mar=c(2,6,2,1))
#pty
#A character specifying the type of plot region to be used; "s" generates a square plotting region and "m" generates the maximal plotting region
par(pty="s")

for(i in 1:length(species.c.list)){
  plotted.data<-as.data.frame(species.c.list[i])
  plotted.data$log_canopy<-log(plotted.data$cc_percentage+1)
  plot(plotted.data$Germination_rate~plotted.data$log_canopy, pch=19, col="grey60", ylab="Germination rate", xlab="", cex.lab=2, cex.axis=2.00,tck=-0.01)
  mtext(paste(letters[i], ")", sep=""), side=2,line=1,adj=1.5,las=1, padj=-13, cex=1.5)
  title(xlab = "log(Canopy closure + 1)", cex.lab=2)
  title(main=bquote(italic(.(species.name.list[i]))), cex.main=2.5)
  model<-lmer(Germination_rate~log_canopy+(1|Site/Plot), plotted.data)
  #if(summary(model)$coefficients[2,3]<0.05){
    x.for.plot<-plotted.data$log_canopy
    y.for.plot<-plotted.data$Germination_rate
    pred.data<-data.frame(intercept=rep(1, 100), x.for.plot=seq.func(x.for.plot))
    pred<-lmer.predict(model, newdat=pred.data, se.mult=1.96, binom=F, poisson=T)
    plot.CI.func(x.for.plot=seq.func(x.for.plot), pred=pred$y, upper=pred$phi, lower=pred$plo, env.colour="black", env.trans=40, line.colour="black", line.type=1, line.weight=3)
#  }
}
dev.off()
#WHY NO LINES COMING UP???? ###### Not working


