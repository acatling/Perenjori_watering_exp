##### WA Perenjori Watering Experiment
## Alexandra Catling
## LTRE analysis 2023

### Installing package by Hernandez et al. ####
## Chrissy Hernandez working through ESA workshop example
#cmh352@cornell.edu
#install.packages("devtools")
#devtools::install_github("chrissy3815/exactLTRE", force=TRUE)
library(exactLTRE)
library(popbio)
library(tidyverse)
library(ggplot2)
##### Decomposing differences in lambda in presence and absence of neighbours ####
## for species that had a signif affect of nbhs on lambda
# So for ARCA, PLDE and VERO

#### Calculating mean rates for each species ####
## z - seed survival rate from Towers' work
#species_level_seed_fill data frame

## Dormancy rates - from panel_figures file.
#arca - untested. 0.76 mean from Maia's work
#hygl - 0.60. Maia
#laro - 0.77. Maia
#peai - 0.38. John
#plde - untested. 0.76 mean from Maia's work
#pole - 0.00. John
#trcy - 0.93 Maia and 0.84 John - average 0.89
#tror - 0.85 Maia and 0.84 John - average 0.84
#vero - 0.65 Maia

#### Matrix elements in order of a11, a21, a12, a22 where:

#a: dormancy rate of seed bank seeds
#z: seed survival rate (maintaining viability)
#f: number of viable seeds produced per adult
#g: rate of seedling emergence
#s: plant survivorship to reproduction

#a11 = az
#a12 = afz
#a21 = (1-a)zgs
#a22 = (1-a)zgsf

###### Calculating mean vital rates ####
#Rates common to both presence and absence of neighbours: 
# dormancy rate (a), emergence rate (g), seed survival (z)

# Using fixed design LTRE with a directional analysis
#with absence of neighbours (intrinsic lambda) as the reference matrix
#ordered as ⁠[reference matrix, test matrix⁠]

#My lambda does not account for seed dormancy:
#my_lambda <- seed_survival*(1-plot_germ)+plot_fecundity*plot_survival*plot_germ
#my_lambda <- z*(1-g)+f*s*g

## Need to calculate average germination rates and 
#fecundity and survival rates in the presence or absence of neighbours

sp_ltre <- popdata %>% group_by(Species, Neighbours01) %>%
  summarise(avg_emerg = mean(plot_germ, na.rm=T),
            avg_surv = mean(plot_survival, na.rm=T),
            avg_fecundity = mean(plot_fecundity, na.rm=T),
            lambda = mean(lambda, na.rm=T))

# arca_ltre <- lambdaarca %>% group_by(Neighbours01) %>% 
#   summarise(avg_germ = mean(plot_germ),
#             avg_surv = mean(plot_survival),
#             avg_fecundity = mean(plot_fecundity),
#             lambda = mean(lambda))
# arca_ltre$seed_survival <- 0.91
# 
# ggplot(lambdaarca, aes(x = Neighbours01, y = lambda))+
#   geom_boxplot()+
#   geom_jitter(alpha = 0.3)+
#   theme_classic()


### ARCA ####
# ARCA - Neighbours0 - absence of neighbours
#a: 0.76 (average from other species in Maia's 2016 Perenjori dataset)
#z: 0.91
#g: 0.24
#f: 10.17
#s: 0.48
# ARCA - Neighbours1 - presence of neighbours
#a: 0.76 (average from other species in Maia's 2016 Perenjori dataset)
#z: 0.91
#g: 0.24
#f: 6.23
#s: 0.20

#Build matrices using rates for absence or presence of neighbours
#A11, A21, A12, A22
#data=c(az, (1-a)zgs, afz, (1-a)zgsf)
arca_no_nbh<- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.48, 0.76*10.17*0.91, (1-0.76)*0.91*0.24*0.48*10.17), nrow=2, ncol=2)
arca_nbh <- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.20, 0.76*6.23*0.91, (1-0.76)*0.91*0.24*0.20*6.23), nrow=2, ncol=2)
#Calculate eigenvalues for each matrix (lambda, largest eigenvalue)
eigen(arca_no_nbh) #0.95
eigen(arca_nbh) #0.76

##Using exactLTRE to perform fixed directional LTRE analysis with no nbh as reference matrix
#Calculating contributions to the difference in lambda due to vital rates
#calcuated by setting all other parameters to their mean

cont_diff <- exactLTRE(list(arca_no_nbh, arca_nbh), method = 'fixed', fixed.directional = TRUE)
# matrix of contributions
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# sum of contributions
sum(cont_diff$epsilons) 
# true difference in lambda  #-0.191
lamDiff(list(arca_no_nbh,arca_nbh)) 

#approximate LTRE approach and popbio LTRE give the same values
#ltre(treatment, ref)
popbio::LTRE(arca_nbh, arca_no_nbh)
approximateLTRE(list(arca_no_nbh,arca_nbh), method='fixed') # contributions to the difference in lambda
#Not confident interpreting exact LTREs atm
## In the update approximateLTRE was renamed classicalLTRE

#Epiphany! My matrix elements represent transition stages, not pure probability of e.g. survival or fecundity
# so there are two things I can report:
#1) Which transition has the biggest influence on the difference between lambdas
#2) If I hold each rate (survival or seed production) in turn at its mean, I will be able to see
# how much each of them alone contributes to each transition, and SUM them!
#Survival affects A21 and A22 and seed production affects A12 and A22
#The sum will then tell me what the contribution to the diff in lambda of that vital rate is
#I can use exact LTRE to see contributions of surv and seed production
#Can I? Denominator is typical true lambda diff, and then sum contributions when just one rate varies?

#### Results:

##All the same:
#The presence of neighbours most affects the transition stage of a dormant seed remaining viable
# emerging and surviving to reproduction
#A21: The effect of neighbours has the largest influence on population growth rate
# via the transition from a dormant, viable seed emerging and surviving to reproduction.
#The contribution to the difference in lambda of the effect of neighbours on 
#a dormant, viable seed emerging is surviving to reproduction is -0.100/-0.191 = 0.52 = 52%
##

##Survival varying only
arca_no_nbh2<- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.48, 0.76*10.17*0.91, (1-0.76)*0.91*0.24*0.48*10.17), nrow=2, ncol=2)
arca_nbh2 <- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.20, 0.76*10.17*0.91, (1-0.76)*0.91*0.24*0.20*10.17), nrow=2, ncol=2)
lamDiff(list(arca_no_nbh2,arca_nbh2))
approximateLTRE(list(arca_no_nbh2,arca_nbh2), method='fixed') # contributions to the difference in lambda
LTRE(arca_nbh2, arca_no_nbh2)
#exactLTRE
cont_diff <- exactLTRE(list(arca_no_nbh2,arca_nbh2), method='fixed', fixed.directional=TRUE) # contributions to the difference in lambda
# # matrix of contributions
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# # sum of contributions
sum(cont_diff$epsilons)

## Contribution of survival is -0.149/-0.191? 78%
# Contribution of transition a21 is -0.118/-0.191? 62%
# As difference between surv varying only dataframe / true lambda diff both varying

#The contribution to the difference in lambda of the effect of neighbours on 
#seed production is -0.059/-0.191 = 0.31 = 31%

#The interaction of survival and seed production is -0.191 - (-0.100 + -0.059)
# = -0.032. -0.032/-0.191 = 0.17 = 17%
#diff lambda - (contribution survival + contribution fecundity)
#The contribution of this interaction is positive: it counteracts the negative
# effects of seed production and survival separately.
#A decrease in seed production has a bigger effect when survival is high than when it is low
#When survival is low, the negative impact of seed production is less severe

##Fecundity varying only
arca_no_nbh3<- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.48, 0.76*10.17*0.91, (1-0.76)*0.91*0.24*0.48*10.17), nrow=2, ncol=2)
arca_nbh3 <- matrix(data=c(0.76*0.91, (1-0.76)*0.91*0.24*0.48, 0.76*6.23*0.91, (1-0.76)*0.91*0.24*0.48*6.23), nrow=2, ncol=2)
lamDiff(list(arca_no_nbh3,arca_nbh3))
approximateLTRE(list(arca_no_nbh3,arca_nbh3), method='fixed') # contributions to the difference in lambda
LTRE(arca_nbh3, arca_no_nbh3)
#exactLTRE
cont_diff <- exactLTRE(list(arca_no_nbh3,arca_nbh3), method='fixed', fixed.directional=TRUE) # contributions to the difference in lambda
# # matrix of contributions
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# # sum of contributions
sum(cont_diff$epsilons)

##Contribution of seed production varying is -0.099/-0.191?  52%?
# hm. Contribution of a12 is -0.076/-0.191? 40%?
##a12+a22 contributes add to 100% but sum of surv and seed production doesn't, likely because of interaction

#diff lambda - (contribution survival + contribution fecundity)
#-0.191 - (-0.149 + -0.099) = 6%...

# result<- exactLTRE(list(arca_no_nbh, arca_nbh), method='fixed', fixed.directional = TRUE)
# result$varying.indices.list
# barplot(t(result$epsilons[2:length(result$epsilons)]),
#         names.arg=c("sJ", "f","sJ, f", "sA", "sJ, sA", "f, sA", "sJ, f, sA"), las=2)

### PLDE ####
# PLDE - Neighbours0 - absence of neighbours
#a: 0.76 (average from other species in Maia's 2016 Perenjori dataset)
#z: 0.98
#g: 0.33
#f: 8.02
#s: 0.25
# PLDE - Neighbours1 - presence of neighbours
#a: 0.76 (average from other species in Maia's 2016 Perenjori dataset)
#z: 0.98
#g: 0.33
#f: 2.55
#s: 0.20

# plde_ltre <- lambdaplde %>% filter(!is.na(lambda)) %>% group_by(Neighbours01) %>% 
#   summarise(avg_germ = mean(plot_germ),
#             avg_surv = mean(plot_survival),
#             avg_fecundity = mean(plot_fecundity),
#             lambda = mean(lambda))
# plde_ltre$seed_survival <- 0.98

#A11, A21, A12, A22
#data=c(az, (1-a)zgs, afz, (1-a)zgsf)
plde_no_nbh<- matrix(data=c(0.76*0.98, (1-0.76)*0.98*0.33*0.25, 0.76*8.02*0.98, (1-0.76)*0.98*0.33*0.25*8.02), nrow=2, ncol=2)
plde_nbh <- matrix(data=c(0.76*0.98, (1-0.76)*0.98*0.33*0.20, 0.76*2.55*0.98, (1-0.76)*0.98*0.33*0.20*2.55), nrow=2, ncol=2)
eigen(plde_no_nbh) #0.90
eigen(plde_nbh) #0.78
# true difference in lambda  #-0.116
lamDiff(list(plde_no_nbh,plde_nbh)) 

# cont_diff <- exactLTRE(list(plde_no_nbh, plde_nbh), method = 'fixed')
# # matrix of contributions
# cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# # sum of contributions
# sum(cont_diff$epsilons) 

LTRE(plde_nbh, plde_no_nbh)
approximateLTRE(list(plde_no_nbh,plde_nbh), method='fixed') # contributions to the difference in lambda

#### Results:
#The contribution to the difference in lambda of the effect of neighbours on 
#survival is -0.0183/-0.116 = 0.16 = 16%

#The contribution to the difference in lambda of the effect of neighbours on 
#seed production is -0.085/-0.116 = 0.73 = 73%

#The interaction of survival and seed production is -0.116 - (-0.0183 + -0.085)
# = -0.013. The interaction is negative. -0.013/-0.116 = 0.11 = 11%
#diff lambda - (contribution survival + contribution fecundity)

### VERO ####
# vero_ltre <- lambdavero %>% group_by(Neighbours01) %>% 
#   summarise(avg_germ = mean(plot_germ),
#             avg_surv = mean(plot_survival),
#             avg_fecundity = mean(plot_fecundity),
#             lambda = mean(lambda))
# vero_ltre$seed_survival <- 0.96

# VERO - Neighbours0 - absence of neighbours
#a: 0.65
#z: 0.96
#g: 0.29
#f: 16.94
#s: 0.69
# VERO - Neighbours1 - presence of neighbours
#a: 0.65
#z: 0.96
#g: 0.29
#f: 11.02
#s: 0.53

#A11, A21, A12, A22
#data=c(az, (1-a)zgs, afz, (1-a)zgsf)
#A22 probability greater than 1 here... 
vero_no_nbh<- matrix(data=c(0.65*0.96, (1-0.65)*0.96*0.29*0.69, 0.65*16.94*0.96, (1-0.65)*0.96*0.29*0.69*16.94), nrow=2, ncol=2)
vero_nbh <- matrix(data=c(0.65*0.96, (1-0.65)*0.96*0.29*0.53, 0.65*11.02*0.96, (1-0.65)*0.96*0.29*0.53*11.02), nrow=2, ncol=2)
eigen(vero_no_nbh) #1.76
eigen(vero_nbh) #1.19
# true difference in lambda  #-0.570
lamDiff(list(vero_no_nbh,vero_nbh)) 

LTRE(vero_nbh, vero_no_nbh)
approximateLTRE(list(vero_no_nbh,vero_nbh), method='fixed') # contributions to the difference in lambda

# cont_diff <- exactLTRE(list(vero_no_nbh, vero_nbh), method = 'fixed')
# # matrix of contributions
# cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
# # sum of contributions
# sum(cont_diff$epsilons) 

#### Results:
#The contribution to the difference in lambda of the effect of neighbours on 
#survival is -0.093/-0.570 = 0.16 = 16%

#The contribution to the difference in lambda of the effect of neighbours on 
#seed production is -0.151/-0.570 = 0.26 = 26%

#The interaction of survival and seed production is -0.116 - (-0.0183 + -0.085)
# = -0.013. The interaction is negative. -0.330/-0.570 = 0.58 = 58%
#diff lambda - (contribution survival + contribution fecundity)

### Old results (I think not correct phrasing) ####
##ARCA:
#The contribution to the difference in lambda of the effect of neighbours on 
#seed production is -0.059/-0.191 = 0.31 = 31%

#The interaction of survival and seed production is -0.191 - (-0.100 + -0.059)
# = -0.032. -0.032/-0.191 = 0.17 = 17%
#diff lambda - (contribution survival + contribution fecundity)
#The contribution of this interaction is positive: it counteracts the negative
# effects of seed production and survival separately.
#A decrease in seed production has a bigger effect when survival is high than when it is low
#When survival is low, the negative impact of seed production is less severe

##PLDE:
#### Results:
#The contribution to the difference in lambda of the effect of neighbours on 
#survival is -0.0183/-0.116 = 0.16 = 16%

#The contribution to the difference in lambda of the effect of neighbours on 
#seed production is -0.085/-0.116 = 0.73 = 73%

#The interaction of survival and seed production is -0.116 - (-0.0183 + -0.085)
# = -0.013. The interaction is negative. -0.013/-0.116 = 0.11 = 11%
#diff lambda - (contribution survival + contribution fecundity)

##VERO:
#### Results:
#The contribution to the difference in lambda of the effect of neighbours on 
#survival is -0.093/-0.570 = 0.16 = 16%

#The contribution to the difference in lambda of the effect of neighbours on 
#seed production is -0.151/-0.570 = 0.26 = 26%

#The interaction of survival and seed production is -0.116 - (-0.0183 + -0.085)
# = -0.013. The interaction is negative. -0.330/-0.570 = 0.58 = 58%
#diff lambda - (contribution survival + contribution fecundity)

#Comparing to my lambda which does not account for seed dormancy!
#my_lambda <- seed_survival*(1-plot_germ)+plot_fecundity*plot_survival*plot_germ
#my_lambda <- z*(1-g)+f*s*g
my_lambda_no_nbh <- 0.91*(1-0.24)+10.17*0.48*0.24 #1.86 (mean says 2.16)
my_lambda_nbh <- 0.91*(1-0.24)+6.23*0.20*0.24 #0.99 (mean says 1.71)
#diff
my_lambda_no_nbh-my_lambda_nbh

### Making plot of transition contributions for supp 8 ####
ltre_sp <- c("ARCA", "PLDE", "VERO")
ltre_sp <- as.data.frame(ltre_sp)
#Add in contribution values
ltre_sp$a12 <- c(0.31, 0.73, 0.26)
ltre_sp$a21 <- c(0.52, 0.16, 0.16)
ltre_sp$a22 <- c(0.17, 0.11, 0.58)

#Need long form
ltre_sp <- ltre_sp %>% pivot_longer(cols = c(a12, a21, a22), names_to = "transition", values_to = "contribution")
str(ltre_sp)
## Make barchart
dev.off()
pdf("Output/Figures/supp_8_ltre.pdf", width=8, height=6)
ggplot(ltre_sp, aes(x=transition, y=contribution))+
  geom_bar(aes(fill=ltre_sp), position = "dodge2", stat = "identity")+
  ylab("Contribution to difference in population growth rates")+
  xlab("Transition")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))+
  theme_classic()+
  guides(fill = guide_legend(title = "Species", label.hjust = 0, title.hjust=0.5)) +
  scale_fill_discrete(labels=c(expression(italic('A. calendula')), expression(italic('P. debilis')), expression(italic('G. rosea'))))+
  theme(axis.text = element_text(size=14),
       legend.text = element_text(size=14), legend.title = element_text(size=14),
       axis.title=element_text(size=14))
dev.off()

## Actually, vero has an interaction so no signif main effect of neighbours on lambda!
ltre_sp_no_vero <- ltre_sp %>% filter(!(ltre_sp=="VERO"))
dev.off()
pdf("Output/Figures/supp_8_ltre_no_vero.pdf", width=8, height=6)
ggplot(ltre_sp_no_vero, aes(x=transition, y=contribution))+
  geom_bar(aes(fill=ltre_sp), position = "dodge2", stat = "identity")+
  ylab("Contribution to difference in population growth rates")+
  xlab("Transition")+
  scale_fill_manual(values = c("#E69F00", "#56B4E9"), labels=c(expression(italic('A. calendula')), expression(italic('P. debilis'))))+
  theme_classic()+
  guides(fill = guide_legend(title = "Species", label.hjust = 0, title.hjust=0.5)) +
  theme(axis.text = element_text(size=14),
        legend.text = element_text(size=14), legend.title = element_text(size=14),
        axis.title=element_text(size=14))
dev.off()
