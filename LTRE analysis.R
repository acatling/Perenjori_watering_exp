## Chrissy Hernandez working through ESA workshop example
#cmh352@cornell.edu


#install.packages("devtools")
#devtools::install_github("chrissy3815/exactLTRE", force=TRUE)
library(exactLTRE)

### My attempt:
#Fixed design LTRE with a directional analysis
#From my analyses:
#lambda = seed_survival*(1-plot_germ)+plot_fecundity*plot_survival*plot_germ


#Only looking at three species with signif. effect of neighbours
## Need to calculate average germination rates and 
#fecundity and survival rates in the presence or absence of neighbours

arca_ltre <- lambdaarca %>% group_by(Neighbours01) %>% 
                                summarise(avg_germ = mean(plot_germ),
                                          avg_surv = mean(plot_survival),
                                          avg_fecundity = mean(plot_fecundity),
                                          lambda = mean(lambda))
arca_ltre$seed_survival <- 0.91

ggplot(lambdaarca, aes(x = Neighbours01, y = lambda))+
  geom_boxplot()+
  geom_jitter(alpha = 0.3)+
  theme_classic()
#Common to both: 
# germination rate (g), seed survival (v)

### ARCA - Neighbours0
#v: 0.91
#g: 0.24
#f: 10.17
#s: 0.48
### ARCA - Neighbours1
#v: 0.91
#g: 0.24
#f: 6.23
#s: 0.20

## In order of v, g, f, s (a11, a12, a21, a22)
A1_no_nbh<- matrix(data=c(0.91, 0.24, 10.17, 0.48), nrow=2, ncol=2)
A2_nbh<- matrix(data=c(0.91, 0.24, 6.23, 0.20), nrow=2, ncol=2)
#my_lambda <- seed_survival*(1-plot_germ)+plot_fecundity*plot_survival*plot_germ
my_lambda_no_nbh <- 0.91*(1-0.24)+10.17*0.48*0.24 #1.86 (mean says 2.16)
my_lambda_nbh <- 0.91*(1-0.24)+6.23*0.20*0.24 #0.99 (mean says 1.71)


# Approximate LTRE:
# contributions to the difference in lambda
cont_diff<- approximateLTRE(list(A1_no_nbh,A2_nbh), method='fixed')
cont_diff; # matrix of contributions 
sum(cont_diff); # sum of approximate contributions 
lamDiff(list(A1_no_nbh,A2_nbh)); # true difference in lambda 

## Difference from my lambdas in 0.87, whereas this says 0.5.

# Exact LTRE:
# contributions to the difference in lambda
#For a fixed design LTRE, exactly 2 matrices must be provided, ordered as  ⁠[reference matrix, test matrix]
#well, no neighbours (thinned) is our treatment. So really with neighbours is the control (natural), and no neighbours is the treatment...
#just a bit different from intrinsic growth rate, then + neighbours
#would need to swap order to reference, treatment (Thinned) then
cont_diff<- exactLTRE(list(A1_no_nbh,A2_nbh), method='fixed')
test<-exactLTRE_fixed(list(A1_no_nbh,A2_nbh)) #same as above line
cont_diff<-exactLTRE_fixed(list(A2_nbh, A1_no_nbh), fixed.directional = TRUE) #slightly different!
#fixed.directional:
#A true/false switch that allows the user to specify whether a directional LTRE should be used. 
#The default behavior is to calculate a symmetric LTRE, where the mean matrix is used as the baseline. See details for more guidance.
## BUT really interesting notes and examples about when to use false and true
#I should probably not use directional design because 'it is not entirely obvious which population should be the reference and which should be the test'

cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
sum(cont_diff$epsilons); # sum of contributions 
lamDiff(list(A2_nbh, A1_no_nbh)); # true difference in lambda 

### So the interaction is 
#Difference in lambda - contribution of each contribution individually (from approximations, -0.4424076)
#So contribution of the interactions between survival and fecundity is -0.4437571--0.4424076 = -0.0013495
# or 0.3%
#Not sure this is correct!!!

#So without neighbours the lambda is ~2.2 and with neighbours it is ~1.7 (from my averages)
#this lamDiff tells me that with neighbours lambda is 0.44 less! matches.
#this also tells me that:
#the contribution to the change in lambda of the effect of neighbours on ...
#fecundity is -0.33 - about 74% (-0.33/-0.44)*100 of the change in lambda.
#Whereas the contribution from the effect of neighbours on survival is 25%

## still need to wrap my head around the interactions part!! Check slides

## PLDE
#Removing rows with no data (no seeds sown)
plde_ltre <- lambdaplde %>% filter(!is.na(lambda)) %>% group_by(Neighbours01) %>% 
  summarise(avg_germ = mean(plot_germ),
            avg_surv = mean(plot_survival),
            avg_fecundity = mean(plot_fecundity),
            lambda = mean(lambda))
plde_ltre$seed_survival <- 0.98

### PLDE - Neighbours0
#v: 0.98
#g: 0.33
#f: 16.05
#s: 0.50
### PLDE - Neighbours1
#v: 0.98
#g: 0.33
#f: 5.10
#s: 0.40

## In order of v, g, f, s (a11, a12, a21, a22)
A1_no_nbh<- matrix(data=c(0.98, 0.33, 16.05, 0.50), nrow=2, ncol=2)
A2_nbh<- matrix(data=c(0.98, 0.33, 5.10, 0.40), nrow=2, ncol=2)
# Approximate LTRE:
# contributions to the difference in lambda
cont_diff<- approximateLTRE(list(A1_no_nbh,A2_nbh), method='fixed')
cont_diff; # matrix of contributions 
sum(cont_diff); # sum of approximate contributions 
lamDiff(list(A1_no_nbh,A2_nbh)); # true difference in lambda 

# Exact LTRE:
# contributions to the difference in lambda
cont_diff<- exactLTRE(list(A1_no_nbh,A2_nbh), method='fixed')
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
sum(cont_diff$epsilons); # sum of contributions 
lamDiff(list(A1_no_nbh,A2_nbh)); # true difference in lambda 

#I would have expected a lambda diff of more like 1.84, this diff says 1.
#the contribution to the change in lambda of the effect of neighbours on ...
#fecundity is -0.96 - about 93%  of the change in lambda
#Whereas the contribution from the effect of neighbours on survival is 4%

### VERO 
vero_ltre <- lambdavero %>% group_by(Neighbours01) %>% 
  summarise(avg_germ = mean(plot_germ),
            avg_surv = mean(plot_survival),
            avg_fecundity = mean(plot_fecundity),
            lambda = mean(lambda))
vero_ltre$seed_survival <- 0.96

### VERO - Neighbours0
#v: 0.96
#g: 0.29
#f: 16.94
#s: 0.69
### VERO - Neighbours1
#v: 0.96
#g: 0.29
#f: 11.02
#s: 0.53

## In order of v, g, f, s (a11, a12, a21, a22)
A1_no_nbh<- matrix(data=c(0.96, 0.29, 16.94, 0.69), nrow=2, ncol=2)
A2_nbh<- matrix(data=c(0.96, 0.29, 11.02, 0.53), nrow=2, ncol=2)

my_lambda_no_nbh <- 0.96*(1-0.29)+10.17*0.48*0.24 #1.86 (mean says 2.16) (eigen says 3.05)
my_lambda_nbh <- 0.96*(1-0.29)+6.23*0.20*0.24 #0.99 (mean says 1.71) (eigen says 2.55)

#Eigenvalue for 2x2 matrix
#https://mathworld.wolfram.com/Eigenvalue.html
#lambda_+/-=1/2[(a_(11)+a_(22))+/-sqrt(4a_(12)a_(21)+(a_(11)-a_(22))^2)]

###This is the equation calculating the eigenvalues!
eigen_nbh<- 1/2*((0.96+0.53)+sqrt(4*0.29*11.02+(0.96-0.53)^2))
eigen_no_nbh<- 1/2*((0.96+0.69)+sqrt(4*0.29*16.94+(0.96-0.69)^2))

# Approximate LTRE:
# contributions to the difference in lambda
cont_diff<- approximateLTRE(list(A1_no_nbh,A2_nbh), method='fixed')
cont_diff; # matrix of contributions 
sum(cont_diff); # sum of approximate contributions 
lamDiff(list(A1_no_nbh,A2_nbh)); # true difference in lambda

# Exact LTRE:
# contributions to the difference in lambda
cont_diff<- exactLTRE(list(A1_no_nbh,A2_nbh), method='fixed')
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
sum(cont_diff$epsilons); # sum of contributions 
lamDiff(list(A1_no_nbh,A2_nbh)); # true difference in lambda 

### Trying calculations where matrix is made up of a11, a12, a21 and a22 as indicated by Jarry et al. 1995:
#Each time step is likely multiplications between vital rates


#I would have expected a lambda diff of 3.5, this is 0.5
#the contribution to the change in lambda of the effect of neighbours on ...
#fecundity is -0.42 - about 85%  of the change in lambda
#Whereas the contribution from the effect of neighbours on survival is 15%

#Fecundity contributes 75-93% of difference, survival contributes 4-25% 

## Why differences in lambda between how I calculate it and them,
#Why difference in lambda diff in particular

## We have a significant lambda interaction between neighbours:PC1 for 4/8 species
#Could I look at the ground squirrel plots Oli et al. 2001
# difference in vital rate (no nbhs and nbhs, fecundity and survival) ~ PC1 and 
#contribution to the difference in lambda (no nbhs and nbhs, fecundity and survival) ~ PC1

#As PC1 changes, how does the contribution of survival and fecundity to the difference
# in lambda between no nbh and nbh plots vary?
#e.g. in low PC1, survival contributes more to the difference in lambdas

## Which life stage transitions are responsible for the difference between nbh vs no nbh lambdas?
#Maybe it changes depending upon the (abiotic) environment


#Could plot this 24 times, for each PC1 plot.

#lambda not significantly influenced by PC1



######### Examples from ESA workshop #####
## Examples from the package documentation:
# Build some example matrices:
A1<- matrix(data=c(0,0.8,0, 0,0,0.7, 5,0,0.2), nrow=3, ncol=3)
A2<- matrix(data=c(0,0.9,0, 0,0,0.5, 4,0,0.3), nrow=3, ncol=3)
A3<- matrix(data=c(0,0.4,0, 0,0,0.6, 6,0,0.25), nrow=3, ncol=3)

# Approximate LTRE:
# contributions to the difference in lambda
cont_diff<- approximateLTRE(list(A1,A2), method='fixed')
cont_diff; # matrix of contributions 
sum(cont_diff); # sum of approximate contributions 
lamDiff(list(A1,A2)); # true difference in lambda 

# contributions to the variance of lambda
cont_var<- approximateLTRE(list(A1,A2,A3), method='random')
round(cont_var,digits=5); # matrix of contributions 
sum(cont_var); # sum of contributions 
lamVar(list(A1,A2,A3)); # true variance in lambda 

# Exact LTRE:
# contributions to the difference in lambda
cont_diff<- exactLTRE(list(A1,A2), method='fixed')
cbind(as.character(cont_diff$varying.indices.list),round(cont_diff$epsilons,digits=5))
sum(cont_diff$epsilons); # sum of contributions 
lamDiff(list(A1,A2)); # true difference in lambda 

# contributions to the variance of lambda
cont_var<- exactLTRE(list(A1,A2,A3), method='random')
cbind(as.character(cont_var$varying.indices.list),round(cont_var$epsilons,digits=5)); 
sum(cont_var$epsilons); # sum of contributions 
lamVar(list(A1,A2,A3)); # true variance in lambda 

# Generation Time:
F1<- matrix(0, nrow=3, ncol=3)
F1[1,3]<- A1[1,3]
#F1 is all zeros, except the upper right corner which matches A1 for adult fertility
gen_time <- generation_time(A1, F1)

# R0, the expected lifetime reproductive output for an individual
R0<- r_nought(A1, F1)

# expected lifespan
U1<- A1
U1[1,3]<- 0
# the upper right corner represents adult fertility in this model. U1, the
# survival matrix, contains all the transitions *except* for fertility.
eta<- lifespan(U1, all_ages=TRUE)
eta_1<- lifespan(U1, all_ages=FALSE) # eta_1 should match the first entry of eta






