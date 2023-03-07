### WA Perenjori Experiment ###
### Alexandra Catling ##
### Remaking large panel figures for paper with interactions

#No neighbours standardised is -0.7948
#With neighbours - one standard deviation above average (average really brought down by zeros), so 1

###PC1 without interactions or significance asterisks ####
dev.off()
pdf("Output/Figures/panel_PC1_test.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(5, 20, 5, 1), mar =c(2,10,1,1))
#ARCA
x_to_plot<-seq.func(arcadata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcagermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcasurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedarca)
model <- arcaseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaarca)
model <- arcalambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#HYGL
x_to_plot<-seq.func(hygldata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglgermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedhygl)
model <- hyglseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdahygl)
model <- hygllambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#LARO
x_to_plot<-seq.func(larodata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, larodata)
model <- larogermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, larodata)
model <- larosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedlaro)
model <- laroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdalaro)
model <- larolambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#PEAI
x_to_plot<-seq.func(peaidata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, peaidata)
model <- peaigermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, peaidata)
model <- peaisurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedpeai)
model <- peaiseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdapeai)
model <- peailambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot^2, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#PLDE
x_to_plot<-seq.func(pldedata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, pldedata)
model <- pldegermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, pldedata)
model <- pldesurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot^2, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedplde)
model <- pldeseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaplde)
model <- pldelambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#TRCY
x_to_plot<-seq.func(trcydata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcygermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcysurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0^2, 0*x_to_plot, 0*x_to_plot, 0*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedtrcy)
model <- trcyseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatrcy)
model <- trcylambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0*0, 0*0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#TROR
x_to_plot<-seq.func(trordata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorgermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedtror)
model <- trorseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatror)
model <- trorlambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#VERO
x_to_plot<-seq.func(verodata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verogermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedvero)
model <- veroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot^2), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdavero)
model <- verolambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

###Overall text
##x labels
mtext("PC1 (std)", adj = 0.11, side = 1, line = 3, cex = 3,outer = TRUE)
mtext("PC1 (std)", adj = 0.4, side = 1, line = 3, cex = 3,outer = TRUE)
mtext("PC1 (std)", adj = 0.68, side = 1, line = 3, cex = 3, outer = TRUE)
mtext("PC1 (std)", adj = 0.97, side = 1, line = 3, cex = 3, outer = TRUE)
##y labels
mtext("Probability of emergence", side = 2, cex = 3, outer=TRUE, line=-5)
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line=-40)
mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 3, outer=TRUE, line=-74)
mtext("Population growth rate (log + 1)", side = 2, cex = 3, outer=TRUE, line=-108)
##main labels
mtext("Emergence", outer=TRUE, adj=0.1,side = 3, cex = 3)
mtext("Survival", outer=TRUE, adj=0.39,side = 3, cex = 3)
mtext("Seed production", outer=TRUE, adj = 0.68, side = 3, cex = 3)
mtext("Population growth", outer=TRUE, adj=1, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 29, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 37, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 35, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 52, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 59, side = 3, cex = 2.5, outer = TRUE)

dev.off()

###### PC1 with interactions and significance stars ####

dev.off()
pdf("Output/Figures/panel_PC1_int.pdf", width=21, height=21)
par(mfrow=c(8,4), oma = c(12, 20, 5, 1), mar =c(2,10,1,1))
#ARCA
x_to_plot<-seq.func(arcadata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcagermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival -signif. int.
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01==1, "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcasurvfinalmod
x_to_plot_no_nbh <- seq.func(arcadata$std_PC1[arcadata$Neighbours01=='0'])
x_to_plot_nbh <- seq.func(arcadata$std_PC1[arcadata$Neighbours01=='1'])
#Without nbhs - blue
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, -0.79, 0, 0*x_to_plot, 0*x_to_plot, x_to_plot*-0.79), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#With nbhs - red
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 1*x_to_plot, 0*x_to_plot, x_to_plot*1), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Seed production
x_to_plot<-seq.func(seedarca$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedarca)
model <- arcaseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdaarca$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaarca)
model <- arcalambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#HYGL
x_to_plot<-seq.func(hygldata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglgermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, hygldata)
model <- hyglsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedhygl$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedhygl)
model <- hyglseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdahygl$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdahygl)
model <- hygllambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#LARO
x_to_plot<-seq.func(larodata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, larodata)
model <- larogermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, larodata)
model <- larosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production - signif. int.
x_to_plot<-seq.func(seedlaro$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha(ifelse(Neighbours01==1, "#CC79A7", "#0072B2"), 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedlaro)
model <- laroseedfinalmod
x_to_plot_no_nbh <- seq.func(seedlaro$std_PC1[seedlaro$Neighbours01=='0'])
x_to_plot_nbh <- seq.func(seedlaro$std_PC1[seedlaro$Neighbours01=='1'])
#Without nbhs - blue
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, -0.79, 0, x_to_plot*-0.79), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#With nbhs - red
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, 0, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Population growth - signif. int.
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdalaro)
model <- larolambdafinalmod
x_to_plot_no_nbh <- seq.func(lambdalaro$std_PC1[lambdalaro$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdalaro$std_PC1[lambdalaro$Neighbours01=='Neighbours1'])
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#PEAI
x_to_plot<-seq.func(peaidata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, peaidata)
model <- peaigermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, x_to_plot^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival -signif. int.
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01==1, "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, peaidata)
model <- peaisurvfinalmod
x_to_plot_no_nbh <- seq.func(peaidata$std_PC1[peaidata$Neighbours01=='0'])
x_to_plot_nbh <- seq.func(peaidata$std_PC1[peaidata$Neighbours01=='1'])
#Without nbhs - blue
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, -0.79, 0, x_to_plot*-0.79), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
#With nbhs - red
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, 0, x_to_plot*1), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)
#Seed production
x_to_plot<-seq.func(seedpeai$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedpeai)
model <- peaiseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*0, 0*0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 1.35,y = 220,"*", cex = 10, col = "red")
#Population growth - signif. int.
x_to_plot_no_nbh <- seq.func(lambdapeai$std_PC1[lambdapeai$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdapeai$std_PC1[lambdapeai$Neighbours01=='Neighbours1'])
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdapeai)
model <- peailambdafinalmod
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, x_to_plot^2, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, x_to_plot^2, 0*x_to_plot, 0*x_to_plot, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#PLDE
x_to_plot<-seq.func(pldedata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, pldedata)
model <- pldegermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0, 0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, pldedata)
model <- pldesurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot^2, x_to_plot*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedplde$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedplde)
model <- pldeseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdaplde$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdaplde)
model <- pldelambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#TRCY
x_to_plot<-seq.func(trcydata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcygermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trcydata)
model <- trcysurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0^2, 0*x_to_plot, 0*x_to_plot, 0*0^2, 0*0^2), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedtrcy$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedtrcy)
model <- trcyseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
text(x = 1.35,y = 42,"*", cex = 10, col = "red")
#Population growth - signif. int.
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatrcy)
model <- trcylambdafinalmod
x_to_plot_no_nbh <- seq.func(lambdatrcy$std_PC1[lambdatrcy$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdatrcy$std_PC1[lambdatrcy$Neighbours01=='Neighbours1'])
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0*0, 0*0, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, 0*1, 0*1, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

#TROR
x_to_plot<-seq.func(trordata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorgermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, trordata)
model <- trorsurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedtror$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedtror)
model <- trorseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth
x_to_plot<-seq.func(lambdatror$std_PC1)
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdatror)
model <- trorlambdafinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()

#VERO
x_to_plot<-seq.func(verodata$std_PC1)
#Emergence
plot(percent_germ ~ std_PC1, ylim=c(0, 1), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verogermfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, x_to_plot, 0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
text(x = 1.35,y = 0.85,"*", cex = 10, col = "red")
#Survival
plot(jitter(surv_to_produce_seeds,0.2) ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, verodata)
model <- verosurvfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, 0*x_to_plot, 0*x_to_plot), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot_CI()
#Seed production
x_to_plot<-seq.func(seedvero$std_PC1)
plot(No_viable_seeds_grouped+1 ~ std_PC1, xlim=c(-1.8,1.5), log = "y", pch=19, col=alpha("grey60", 0.3), ylab= NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, seedvero)
model <- veroseedfinalmod
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0, x_to_plot^2), se.mult = 1.96, logit_link=FALSE, log_link=TRUE, glmmTMB=TRUE)
plot_CI()
#Population growth - signif. int.
plot(log_lambda ~ std_PC1, xlim=c(-1.8,1.5), pch=19, col=alpha(ifelse(Neighbours01=='Neighbours1', "#CC79A7", "#0072B2"), 0.3), ylab=NA, xlab=NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, lambdavero)
model <- verolambdafinalmod
x_to_plot_no_nbh <- seq.func(lambdavero$std_PC1[lambdavero$Neighbours01=='Neighbours0'])
x_to_plot_nbh <- seq.func(lambdavero$std_PC1[lambdavero$Neighbours01=='Neighbours1'])
plotted.pred.no.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 0, 0*x_to_plot, 0*x_to_plot, x_to_plot*0), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_no_nbh, pred = plotted.pred.no.nbh$y, upper = plotted.pred.no.nbh$upper, lower = plotted.pred.no.nbh$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
plotted.pred.nbh <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot, 0, 1, 0*x_to_plot, 0*x_to_plot, x_to_plot*1), se.mult = 1.96, logit_link=FALSE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_nbh, pred = plotted.pred.nbh$y, upper = plotted.pred.nbh$upper, lower = plotted.pred.nbh$lower, env.colour = "#CC79A7", env.trans = 50, line.colour = "#CC79A7", line.weight = 2, line.type = 1)

###Overall text
##x labels
mtext("PC1", adj = 0.15, side = 1, line = 3, cex = 3,outer = TRUE)
mtext("PC1", adj = 0.4, side = 1, line = 3, cex = 3,outer = TRUE)
mtext("PC1", adj = 0.66, side = 1, line = 3, cex = 3, outer = TRUE)
mtext("PC1", adj = 0.93, side = 1, line = 3, cex = 3, outer = TRUE)
##y labels
mtext("Probability of emergence", side = 2, cex = 3, outer=TRUE, line=-5)
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line=-40)
mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 3, outer=TRUE, line=-74)
mtext("Population growth rate (log)", side = 2, cex = 3, outer=TRUE, line=-108)
##main labels
mtext("Emergence", outer=TRUE, adj=0.1,side = 3, cex = 3)
mtext("Survival", outer=TRUE, adj=0.39,side = 3, cex = 3)
mtext("Seed production", outer=TRUE, adj = 0.68, side = 3, cex = 3)
mtext("Population growth", outer=TRUE, adj=1, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 4, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 27, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 34, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 34, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 50, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 56, side = 3, cex = 2.5, outer = TRUE)
reset()
legend("bottom", title=NULL, horiz=T, legend=c("No neighbours", "Neighbours"),
       col=c("#0072B2", "#CC79A7"), pch=19, cex=3, bty="n")
dev.off()

#### Neighbours ####

dev.off()
pdf("Output/Figures/panel_NA.pdf", width=21, height=21)
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
mtext("Neighbour abundance", adj = 0.09, side = 1, line=3, cex = 3,outer = TRUE)
mtext("Neighbour abundance", adj = 0.54, side = 1, line=3, cex = 3, outer = TRUE)
mtext("Neighbour presence", adj = 0.98, side = 1, line=3, cex = 3, outer = TRUE)
##y labels
mtext("Probability of survival", side = 2, cex = 3, outer=TRUE, line=-6)
mtext("Number of viable seeds produced (log + 1)", side = 2, cex = 3, outer=TRUE, line=-50.5)
mtext("Population growth rate (log)", side = 2, cex = 3, outer=TRUE, line=-97.5)
##main labels
mtext("Survival", outer=TRUE, adj=0.15,side = 3, cex = 3)
mtext("Seed production", outer=TRUE, adj = 0.55, side = 3, cex = 3)
mtext("Population growth", outer=TRUE, adj=0.96, side = 3, cex = 3)
#Use mxtext for species names
mtext("Species", outer = TRUE, adj = -0.14, side = 3, cex = 3)
mtext(~italic("A. calendula"), adj = -0.15, padj= 5, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("H. glutinosum"), adj = -0.15, padj= 10, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("L. rosea"), adj = -0.15, padj= 20, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. airoides"), adj = -0.15, padj= 29, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("P. debilis"), adj = -0.15, padj= 37, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. cyanopetala"), adj = -0.16, padj= 35, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("T. ornata"), adj = -0.15, padj= 52, side = 3, cex = 2.5, outer = TRUE)
mtext(~italic("G. rosea"), adj = -0.15, padj= 59, side = 3, cex = 2.5, outer = TRUE)

dev.off()

#### Estimating seed dormancy from prior data ####
dormancy_data <- read_csv("Data/wa_seed_dormancy_previous_data.csv")
john_data <- dormancy_data %>% filter(researcher == 'John_Dwyer_and Todd_Erickson')
maia_data <- dormancy_data %>% filter(researcher == 'Maia_Raymundo')
str(maia_data)
maia_data$field_seeds_scaled_viability <- as.numeric(maia_data$field_seeds_scaled_viability)
maia_data$field_number_germinants <- as.numeric(maia_data$field_number_germinants)

#assuming the seeds that didn't germinate are dormant, since already scaled for viability
maia_data <- maia_data %>% select(genus, species, field_number_seeds, field_seeds_scaled_viability, field_number_germinants) %>%
  mutate(number_seeds_dormant = field_seeds_scaled_viability-field_number_germinants)
#setting negative values to zero - this is where more seeds germinated than were estimated to be viable
maia_data$seeds_dormant_no_negs <- ifelse(maia_data$number_seeds_dormant<0, 0, maia_data$number_seeds_dormant)
maia_data <- maia_data %>% mutate(dormancy_rate = seeds_dormant_no_negs/field_seeds_scaled_viability)
#make a column for united genus and species
maia_data <- maia_data %>% unite("species_id", genus:species, remove = "false")
maia_sp_dorm <- maia_data %>% group_by(species_id) %>% summarise(sp_dorm_rate = mean(dormancy_rate, na.rm=T))

#John's
#Number of seeds that are viable would be the number that germinated + the number that were stained and viable
# since ungerminated seeds were all stained
john_data <- john_data %>% unite("species_id", genus:species, remove = "false")
str(john_data)
john_data$chamber_number_germinants <- as.numeric(john_data$chamber_number_germinants)
john_data$TZ_number_viable <- as.numeric(john_data$TZ_number_viable)
john_data <- john_data %>% mutate(number_seeds_viable = chamber_number_germinants+TZ_number_viable,
                                  number_seeds_dormant = number_seeds_viable-chamber_number_germinants,
                                  dormancy_rate = number_seeds_dormant/number_seeds_viable)
john_sp_data <- john_data %>% group_by(species_id) %>% summarise(sp_dorm_rate = mean(dormancy_rate, na.rm=T))

### Dormancy rates for all species tested
#arca - untested. 0.55 mean
#Dunbabin and Cocks 1999 in WA - diff ecotypes varied from 5% to 95%
#hygl - 0.60. Maia
#laro - 0.77. Maia
#peai - 0.38. John
#plde - untested. 0.55 mean
#pole - 0.00. John
#trcy - 0.93 Maia and 0.84 John - average 0.89
#tror - 0.85 Maia and 0.84 John - average 0.84
#vero - 0.65 Maia

## Averages of measured rates: mean(c(0.60, 0.77, 0.38, 0.89, 0.65)) - 0.66
#with pole mean(c(0.60, 0.77, 0.38, 0, 0.89, 0.65)) - 0.55.
#Just of Maia's - mean(c(0.60, 0.77, 0.93, 0.85, 0.65)) - 0.76

## Calculating lambdas for each species
meanlambda <- popdata %>% group_by(Species, Neighbours01) %>% 
  summarise(mean_lambda = mean(lambda, na.rm=T))
meanlambda$log <- log(meanlambda$mean_lambda)

meanlambda <- meanlambda %>% group_by(Species) %>% 
  summarise(difference = mean_lambda[])

test <- meanlambda$mean_lambda[meanlambda$Neighbours01=='Neighbours0']-meanlambda$mean_lambda[meanlambda$Neighbours01=='Neighbours1']

ggplot(meanlambda, aes(x = Species, y = log(mean_lambda), colour = Neighbours01))+
  geom_point()+
  theme_classic()
test <- lambdatror %>% filter(Neighbours01==1)
mean(test$lambda, na.rm=T)
mean(lambdatror$log_lambda[lambdatror$Neighbours01== 'Neighbours1'])  

### Plotting pop growth ~ neighbour abundance by high and low PC1
#for laro, peai, trcy and goro

## How does arca dry:neighbour presence interaction look?
ggplot(lambdavero, aes(x = Neighbours01, y=log(lambda), colour = ifelse(std_PC1<0, 'green', 'purple')))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  theme_classic()

##NA ~ PC1
ggplot(vitaldata, aes(y = std_logp1_totalabund, x = std_PC1))+
  geom_point(alpha=0.3)+
  geom_smooth(method="lm")+
  theme_classic()


#arca
plot(lambda ~ std_logp1_totalabund, col=alpha(ifelse(std_PC1<0, "red", "blue"),0), lambdaarca)

x_to_plot_high <- seq.func(lambdadata$std_PC1[arcadata$Treatment=='Dry'])
x_to_plot_low <- seq.func(lambdadata$std_PC1[arcadata$Treatment=='Ambient'])
plot(surv_to_produce_seeds ~ jitter(std_PC1, 15), xlim=c(-1.8,1.5), pch=19, col=alpha("grey60", 0.3), ylab=NA, xlab = NA, tck=-0.01, cex= 2.5, cex.axis= 2.5, arcadata)
model <- arcasurvfinalmod
#ambient - black
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 0, x_to_plot_amb, 0, 0, 0, 0*x_to_plot_amb, 0*x_to_plot_amb, x_to_plot_amb*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_amb, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "grey60", env.trans = 50, line.colour = "black", line.weight = 2, line.type = 1)
#wet - blue
plotted.pred <- glmm.predict(mod = model, newdat = data.frame(1, 0, 1, x_to_plot_wet, 0, 0, 0, 0*x_to_plot_wet, 1*x_to_plot_wet, x_to_plot_wet*0), se.mult = 1.96, logit_link=TRUE, log_link=FALSE, glmmTMB=FALSE)
plot.CI.func(x.for.plot = x_to_plot_wet, pred = plotted.pred$y, upper = plotted.pred$upper, lower = plotted.pred$lower, env.colour = "#0072B2", env.trans = 50, line.colour = "#0072B2", line.weight = 2, line.type = 1)
