### Producing multi-panel figure of vital rates ~ cover (open or shade)
## Perenjori Watering Experiment
## Alexandra Catling October 2023

library(scales)
library(cowplot)

dev.off()
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(2,10,1,1))
### ARCA ####
coef(arcagermfinalmod)
arca_germ_pred<-glmm.predict(mod=arcagermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
arca_germ_pred$Cover <- c('Shade', 'Sun')
a <- ggplot()+
  geom_jitter(data=arcadata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=arca_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(arcasurvfinalmod)
arca_surv_pred<-glmm.predict(mod=arcasurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
arca_surv_pred$Cover <- c('Shade', 'Sun')
b <- ggplot()+
  geom_jitter(data=arcadata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=arca_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
#arca
#CAN'T go below 0, the estimate from model. troubleshooting below
# coef(arcaseedfinalmod)
# arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
#                              se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=TRUE)
# arca_seed_pred$Cover <- c('Shade', 'Sun')
# arca_seed_pred$lower <- c(0.165081, 0.3698253)
# #Get this lower value from log_link=TRUE, but it does fix the 0 issue
# #or I could manually set it to zero?
# ggplot()+
#   geom_jitter(data=seedarca, aes(x = Cover, y = log(No_viable_seeds_grouped+1)), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
#   geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
#   geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
#   theme_classic()+
#   theme(axis.text=element_text(size=24), 
#         axis.title=element_text(size=24),
#         axis.title.y=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank())
#Using log10 y axis but still below 0 so adding 1
#but adding 1 changes the error bars' position compared to plot above
# coef(arcaseedfinalmod)
# arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
#                              se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
# arca_seed_pred$Cover <- c('Shade', 'Sun')
# arca_seed_pred$lower <- arca_seed_pred$lower+1
# arca_seed_pred$upper <- arca_seed_pred$upper+1
# #add one to the mean??
# ggplot()+
#   geom_jitter(data=seedarca, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
#   geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
#   geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
#   scale_y_continuous(trans="log10")+
#   theme_classic()+
#   theme(axis.text=element_text(size=24),
#         axis.title=element_text(size=24),
#         axis.title.y=element_blank(),
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank())

#plotting on a log scale, so need seeds+1 and +1 to negative lower limits
#exponentiating after logging, log_link=TRUE IS exponentiating it for us
coef(arcaseedfinalmod)
arca_seed_pred<-glmm.predict(mod=arcaseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*0, 0*0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
arca_seed_pred$Cover <- c('Shade', 'Sun')
#added one to the negative lower limit
arca_seed_pred$lower <- c(1.165081, 1.447482)
c <- ggplot()+
  geom_jitter(data=seedarca, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=arca_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

#in the absence of neighbours
coef(arcalambdafinalmod)
arca_lambda_pred<-glmm.predict(mod=arcalambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
arca_lambda_pred$Cover <- c('Shade', 'Sun')
d <- ggplot()+
  geom_jitter(data=lambdaarca, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=arca_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=arca_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### hygl ####
coef(hyglgermfinalmod)
hygl_germ_pred<-glmm.predict(mod=hyglgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_germ_pred$Cover <- c('Shade', 'Sun')
e <- ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=hygl_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 24, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hyglsurvfinalmod)
hygl_surv_pred<-glmm.predict(mod=hyglsurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), 0*0, 0*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
hygl_surv_pred$Cover <- c('Shade', 'Sun')
f <- ggplot()+
  geom_jitter(data=hygldata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=hygl_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hyglseedfinalmod)
hygl_seed_pred<-glmm.predict(mod=hyglseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
hygl_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
g <- ggplot()+
  geom_jitter(data=seedhygl, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=hygl_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(hygllambdafinalmod)
hygl_lambda_pred<-glmm.predict(mod=hygllambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
hygl_lambda_pred$Cover <- c('Shade', 'Sun')
#negative values here are fine! Population growth rate
h <- ggplot()+
  geom_jitter(data=lambdahygl, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=hygl_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=hygl_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### laro ####
coef(larogermfinalmod)
laro_germ_pred<-glmm.predict(mod=larogermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_germ_pred$Cover <- c('Shade', 'Sun')
i <- ggplot()+
  geom_jitter(data=larodata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=laro_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 24, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(larosurvfinalmod)
laro_surv_pred<-glmm.predict(mod=larosurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
laro_surv_pred$Cover <- c('Shade', 'Sun')
j <- ggplot()+
  geom_jitter(data=larodata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=laro_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(laroseedfinalmod)
laro_seed_pred<-glmm.predict(mod=laroseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
laro_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
k <- ggplot()+
  geom_jitter(data=seedlaro, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=laro_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(larolambdafinalmod)
laro_lambda_pred<-glmm.predict(mod=larolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
laro_lambda_pred$Cover <- c('Shade', 'Sun')
l <- ggplot()+
  geom_jitter(data=lambdalaro, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=laro_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=laro_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### peai ####
coef(peaigermfinalmod)
peai_germ_pred<-glmm.predict(mod=peaigermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_germ_pred$Cover <- c('Shade', 'Sun')
m <- ggplot()+
  geom_jitter(data=peaidata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=peai_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 24, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peaisurvfinalmod)
peai_surv_pred<-glmm.predict(mod=peaisurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
peai_surv_pred$Cover <- c('Shade', 'Sun')
n <- ggplot()+
  geom_jitter(data=peaidata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=peai_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peaiseedfinalmod)
peai_seed_pred<-glmm.predict(mod=peaiseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
peai_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
o <- ggplot()+
  geom_jitter(data=seedpeai, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=peai_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(peailambdafinalmod)
peai_lambda_pred<-glmm.predict(mod=peailambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
peai_lambda_pred$Cover <- c('Shade', 'Sun')
p <- ggplot()+
  geom_jitter(data=lambdapeai, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=peai_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=peai_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### plde ####
coef(pldegermfinalmod)
plde_germ_pred<-glmm.predict(mod=pldegermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_germ_pred$Cover <- c('Shade', 'Sun')
q <- ggplot()+
  geom_jitter(data=pldedata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=plde_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldesurvfinalmod)
plde_surv_pred<-glmm.predict(mod=pldesurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plde_surv_pred$Cover <- c('Shade', 'Sun')
r <- ggplot()+
  geom_jitter(data=pldedata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=plde_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldeseedfinalmod)
plde_seed_pred<-glmm.predict(mod=pldeseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plde_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
s <- ggplot()+
  geom_jitter(data=seedplde, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=plde_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(pldelambdafinalmod)
plde_lambda_pred<-glmm.predict(mod=pldelambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plde_lambda_pred$Cover <- c('Shade', 'Sun')
t <- ggplot()+
  geom_jitter(data=lambdaplde, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=plde_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=plde_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### trcy ####
coef(trcygermfinalmod)
trcy_germ_pred<-glmm.predict(mod=trcygermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_germ_pred$Cover <- c('Shade', 'Sun')
u <- ggplot()+
  geom_jitter(data=trcydata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=trcy_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 24, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcysurvfinalmod)
trcy_surv_pred<-glmm.predict(mod=trcysurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0), c(1,0)*0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
trcy_surv_pred$Cover <- c('Shade', 'Sun')
v <- ggplot()+
  geom_jitter(data=trcydata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=trcy_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcyseedfinalmod)
trcy_seed_pred<-glmm.predict(mod=trcyseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, c(1,0)*0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
trcy_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
w <- ggplot()+
  geom_jitter(data=seedtrcy, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=trcy_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trcylambdafinalmod)
trcy_lambda_pred<-glmm.predict(mod=trcylambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0*0, 0*0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
trcy_lambda_pred$Cover <- c('Shade', 'Sun')
x <- ggplot()+
  geom_jitter(data=lambdatrcy, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=trcy_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=trcy_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### tror ####
coef(trorgermfinalmod)
tror_germ_pred<-glmm.predict(mod=trorgermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_germ_pred$Cover <- c('Shade', 'Sun')
y <- ggplot()+
  geom_jitter(data=trordata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=tror_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  geom_text(aes(label = "*"), size = 24, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorsurvfinalmod)
tror_surv_pred<-glmm.predict(mod=trorsurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0)),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
tror_surv_pred$Cover <- c('Shade', 'Sun')
z <- ggplot()+
  geom_jitter(data=trordata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=tror_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorseedfinalmod)
tror_seed_pred<-glmm.predict(mod=trorseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0, 0*c(1,0), 0*c(1,0)),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
tror_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
ab <- ggplot()+
  geom_jitter(data=seedtror, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=tror_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

coef(trorlambdafinalmod)
tror_lambda_pred<-glmm.predict(mod=trorlambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
tror_lambda_pred$Cover <- c('Shade', 'Sun')
bc <- ggplot()+
  geom_jitter(data=lambdatror, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=tror_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=tror_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())

### vero ####
coef(verogermfinalmod)
vero_germ_pred<-glmm.predict(mod=verogermfinalmod, newdat=data.frame(1, c(1,0), 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_germ_pred$Cover <- c('Shade', 'Sun')
cd <- ggplot()+
  geom_jitter(data=verodata, aes(x = Cover, y = percent_germ), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=vero_germ_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_germ_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1), limits = c(0,1))+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  geom_text(aes(label = "*"), size = 24, hjust = 0, y = 0.7, x = 1.4, colour = "red")+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(verosurvfinalmod)
vero_surv_pred<-glmm.predict(mod=verosurvfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
vero_surv_pred$Cover <- c('Shade', 'Sun')
de <- ggplot()+
  geom_jitter(data=verodata, aes(x = Cover, y = surv_to_produce_seeds), alpha = 0.1, col = "grey10", width=0.05, height = 0, cex = 2)+
  geom_point(data=vero_surv_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_surv_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(veroseedfinalmod)
vero_seed_pred<-glmm.predict(mod=veroseedfinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, 0),
                             se.mult=1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
vero_seed_pred$Cover <- c('Shade', 'Sun')
#no negative lower limit
ef <- ggplot()+
  geom_jitter(data=seedvero, aes(x = Cover, y = No_viable_seeds_grouped+1), alpha = 0.1, col = "grey10", width=0.05, height=0,cex = 2)+
  geom_point(data=vero_seed_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_seed_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_y_continuous(trans="log10")+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())

coef(verolambdafinalmod)
vero_lambda_pred<-glmm.predict(mod=verolambdafinalmod, newdat=data.frame(1, 0, 0, c(1,0), 0, 0, 0, c(1,0)*0),
                               se.mult=1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
vero_lambda_pred$Cover <- c('Shade', 'Sun')
fg <- ggplot()+
  geom_jitter(data=lambdavero, aes(x = Cover, y = log_lambda), alpha = 0.1, col = "grey10", width=0.05, height=0, cex = 2)+
  geom_point(data=vero_lambda_pred, aes(x=Cover, y=y), cex=3) +
  geom_errorbar(data=vero_lambda_pred, aes(x=Cover, ymin = lower, ymax = upper, width = 0.15), cex=1.5)+
  theme_classic()+
  scale_x_discrete(labels=c("Open", "Shade"))+
  theme(axis.text=element_text(size=24), 
        axis.title=element_text(size=24),
        axis.title.x=element_text(size=34),
        axis.text.x=element_text(size=34),
        axis.title.y=element_blank())
#### plotting cover_panel from above ####
pdf("Output/Figures/panel_cover.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(2,10,1,1))

plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,ab,bc,cd,de,ef,fg, align="hv", ncol=4)
dev.off()
#hjust=-5.5, 

pdf("Output/Figures/panel_cover_labels.pdf", width=21, height=21)
plot_row <- plot_grid(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,ab,bc,cd)
title <- ggdraw() + 
  draw_label(
    "Probability of emergence",
   # fontface = 'bold',
    x = 0,
    hjust = 0,
   angle = 90
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 7)
  )
plot_grid(
  title, plot_row,
  ncol = 4, align="hv"
)
dev.off()
###Overall text
##x labels
# mtext("PC1 (std)", adj = 0.11, side = 1, line = 3, cex = 2,outer = TRUE)
# mtext("PC1 (std)", adj = 0.4, side = 1, line = 3, cex = 2,outer = TRUE)
# mtext("PC1 (std)", adj = 0.68, side = 1, line = 3, cex = 2, outer = TRUE)
# mtext("PC1 (std)", adj = 0.97, side = 1, line = 3, cex = 2, outer = TRUE)
# ##y labels
# mtext("Probability of emergence", side = 2, cex = 2, outer=TRUE, line=-5)
# mtext("Probability of survival", side = 2, cex = 2, outer=TRUE, line=-40)
# mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 2, outer=TRUE, line=-74)
# mtext("Population growth rate (log + 1)", side = 2, cex = 2, outer=TRUE, line=-108)
# ##main labels
# mtext("Emergence", outer=TRUE, adj=0.1,side = 3, cex = 2)
# mtext("Survival", outer=TRUE, adj=0.39,side = 3, cex = 2)
# mtext("Seed production", outer=TRUE, adj = 0.68, side = 3, cex = 2)
# mtext("Population growth", outer=TRUE, adj=1, side = 3, cex = 2)
# #Use mxtext for species names
# mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 2)
# mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("P. airoides"), adj = -0.15, padj= 29, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("P. debilis"), adj = -0.15, padj= 37, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 35, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("T. ornata"), adj = -0.15, padj= 52, side = 3, cex = 2.5, outer = TRUE)
# mtext(~italic("G. rosea"), adj = -0.15, padj= 59, side = 3, cex = 2.5, outer = TRUE)

#### Trying plotting cover in base R instead ####
#add species name
amyg_pred_dbh$Focal_sp <- 'E. amygdalina'

arca_germ_pred$Species <- 'A. calendula'
#0.7, 1.2
arca_germ_pred$x <- c(0.1, 0.2)

dev.off()
pdf("Output/panel_predicted_growth.pdf", width=21, height=21)
par(oma=c(12,8,3,1), mfrow=c(3,1), mar=c(1, 4, 4, 1))

plot(y ~ x, pch = 19, ylab ="", ylim = c(0,1), xlab="", tck=-0.02, cex=2, cex.axis = 2, xaxt="n", bty="l", arca_germ_pred)

plot(y ~ x, pch = ifelse(pred_md_all$md_category == 'wet', 19, 15), ylab="", ylim=c(0,0.11), xlab="", 
     tck=-0.01, cex=6, cex.axis = 5, xaxt="n", bty="l", pred_md_all)
#Add error bars
arrows(x0=pred_md_all$x, y0=pred_md_all$lower, x1=pred_md_all$x, y1=pred_md_all$upper, code=3, angle=90, length=0.1, lwd=5)
legend("topright", title= expression(bold('Long-term MD')), horiz=F, legend=c("Low (wet)", "High (dry)"),
       pch=c(19, 15), cex=4, bty="n")
mtext(expression(bold("A)")), side=1, cex=3.8, adj=0.005, padj=-8.1, line=-3)


#### New plot of neighbours with means and error instead of boxplot####
#### Neighbours ####

dev.off()
pdf("Output/Figures/panel_NA_new.pdf", width=21, height=21)
par(mfrow=c(8,3), oma = c(5, 20, 5, 1), mar =c(3,10,1,1))
#ARCA
x_to_plot<-seq.func(arcadata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab = NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, arcadata)
model <- arcasurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedarca$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1, 150), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedarca)
model <- arcaseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab=NA, xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdaarca)
stripchart(log_lambda ~ Neighbours01, lambdaarca, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex = 2, add = TRUE)
text(x = 1.5, y = 2.55, "*", cex = 10, col = "red")

#Adding HYGL
x_to_plot<-seq.func(hygldata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, hygldata)
model <- hyglsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedhygl$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedhygl)
model <- hyglseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdahygl)
stripchart(log_lambda ~ Neighbours01, lambdahygl, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding LARO
x_to_plot<-seq.func(larodata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, larodata)
model <- larosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedlaro$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedlaro)
model <- laroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdalaro)
stripchart(log_lambda ~ Neighbours01, lambdalaro, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding PEAI
x_to_plot<-seq.func(peaidata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, peaidata)
model <- peaisurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedpeai$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedpeai)
model <- peaiseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdapeai)
stripchart(log_lambda ~ Neighbours01, lambdapeai, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding PLDE
x_to_plot<-seq.func(pldedata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, pldedata)
model <- pldesurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0^2, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedplde$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedplde)
model <- pldeseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 2.3, y = 50, "*", cex = 10, col = "red")
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdaplde)
stripchart(log_lambda ~ Neighbours01, lambdaplde, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)
text(x = 1.5, y = 2, "*", cex = 10, col = "red")

#Adding TRCY
x_to_plot<-seq.func(trcydata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, trcydata)
model <- trcysurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, x_to_plot^2, 0*0, 0*0, 0*x_to_plot^2, 0*x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedtrcy$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedtrcy)
model <- trcyseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdatrcy)
stripchart(log_lambda ~ Neighbours01, lambdatrcy, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding TROR
x_to_plot<-seq.func(trordata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, trordata)
model <- trorsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedtror$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedtror)
model <- trorseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab=NA, names= NA, col = "white", cex= 2.5, cex.axis = 2.5, lambdatror)
stripchart(log_lambda ~ Neighbours01, lambdatror, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)

#Adding VERO
x_to_plot<-seq.func(verodata$std_logp1_totalabund)
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_logp1_totalabund, xlim=c(-1,2.5), pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, verodata)
model <- verosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0*0, 0*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Fecundity
x_to_plot<-seq.func(seedvero$std_logp1_totalabund)
plot(No_viable_seeds_grouped+1 ~ jitter(std_logp1_totalabund, 1), xlim=c(-1,2.5), ylim=c(1,100), log = "y", pch=19, col=alpha("grey60", 0.3), ylab="", xlab=NA, tck=-0.01, cex= 2.5, cex.axis = 2.5, seedvero)
model <- veroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, 0, 0, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Lambda
boxplot(log_lambda ~ Neighbours01, pch=19, ylab="", xlab="", names = c("", ""), col = "white", cex.axis = 2.5, lambdavero)
stripchart(log_lambda ~ Neighbours01, lambdavero, pch = 19, method = "jitter", col=alpha("grey60", 0.6), vertical = TRUE, cex= 2.5, add = TRUE)
mtext(side=1, "Absent", adj=0.1, line=2.5, cex =2.5)
mtext(side=1, "Present", adj=0.9, line=2.5, cex =2.5)

###Overall text
##x labels
mtext("Neighbour abundance", adj = 0.09, side = 1, line=3, cex = 2,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, line=3, cex = 2, outer = TRUE)
mtext("Neighbour presence", adj = 0.98, side = 1, line=3, cex = 2, outer = TRUE)
##y labels
mtext("Probability of survival", side = 2, cex = 2, outer=TRUE, line=-6)
mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 2, outer=TRUE, line=-50.5)
mtext("Population growth rate (log)", side = 2, cex = 2, outer=TRUE, line=-97.5)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 2)
mtext("Seed production", outer=TRUE, adj = 0.55, side = 3, cex = 2)
mtext("Population growth", outer=TRUE, adj=0.96, side = 3, cex = 2)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 2)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 29, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 37, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 35, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 52, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 59, side = 3, cex = 2.5, outer = TRUE)

dev.off()
