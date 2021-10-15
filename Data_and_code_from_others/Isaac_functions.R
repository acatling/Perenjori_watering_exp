#Isaac's functions from his WA paper

lm.predict<-function(mod, newdat){ 
  pred = predict(mod, newdata = newdat, interval = "confidence", level = 0.95)
  upper<-pred[,3]
  lower<-pred[,2] 
  return(data.frame(pred=pred[,1], upper, lower))
}

###function to predict linear mixed effects models
lme.predict<-function(mod, newdat, int.type=NA, se.mult){
  pred <- as.vector(predict(mod, newdat, level = 0))
  Designmat <- model.matrix(formula(mod)[-2], newdat)
  predvar <- diag(Designmat %*% vcov(mod) %*% t(Designmat)) 
  SE <- sqrt(predvar) 
  SE2 <- sqrt(predvar+mod$sigma^2)
  if(int.type=="CI"){
    upper<-pred + (se.mult*SE)
    lower<-pred - (se.mult*SE) 
  }
  if(int.type=="PI"){
    upper<-pred + (se.mult*SE2)
    lower<-pred - (se.mult*SE2)
  }
  
  return(data.frame(pred, upper, lower))
}	

####function to plot the curve and associated condfidence interval
plot.CI.func<- function(x.for.plot, pred, upper, lower, env.colour, env.trans=NA, line.colour, line.weight){
  colour.rgb<-col2rgb(col=env.colour)  
  polygon.coords<-data.frame(rbind(cbind(x.for.plot[1], lower[1]), 
                                   cbind(x.for.plot, upper), 
                                   cbind(x.for.plot, lower)[rev(order(x.for.plot)),]))
  names(polygon.coords)<-c("x", "y")							
  polygon(polygon.coords$x, polygon.coords$y, col=rgb(red=colour.rgb["red",],blue=colour.rgb["blue",], green=colour.rgb["green",] , alpha=env.trans, maxColorValue = 255), border=NA)
  lines(x.for.plot, pred, col=line.colour, lwd=line.weight)         
} 

###function to inverse logits
inv.logit<-function(x)(exp(x)/(1+exp(x)))

##creates a sequence of values frmo min to max by 100
seq.func<-function(x)(seq(min(x, na.rm=T), max(x, na.rm=T), length.out=100))

### plot with 95% CI using lmer  or lme objects
lmer.predict<-function(mod, newdat, se.mult, binom=NULL, poisson=NULL){
  pvar1 <- diag(as.matrix(newdat) %*% tcrossprod(vcov(mod),as.matrix(newdat)))
  newdat$y<- as.matrix(newdat) %*% fixef(mod)  
  newdat <- data.frame(newdat, plo = newdat$y-(se.mult*sqrt(pvar1)), phi = newdat$y+(se.mult*sqrt(pvar1)))
  
  ## if you have used binomial errors then this will back transform logits to the probability scale
  if(binom==T) {
    newdat$y<-plogis(newdat$y); newdat$plo<-plogis(newdat$plo); newdat$phi<-plogis(newdat$phi)
  } else 
    
    ## if you have used poisson errors or have log-transformed your response, then this will back transform to the original scale (e.g. abundance)
    ##  not used in this code
    if(poisson==T) {
      newdat$y<-exp(newdat$y); newdat$plo<-exp(newdat$plo); newdat$phi<-exp(newdat$phi)
    } 
  return(with(newdat, data.frame(y, phi, plo)))
}

plot.CI.func<- function(x.for.plot, pred, upper, lower, env.colour, env.trans=NA, line.colour, line.weight, line.type){
  colour.rgb<-col2rgb(col=env.colour)  
  polygon.coords<-data.frame(rbind(cbind(x.for.plot[1], lower[1]), 
                                   cbind(x.for.plot, upper), 
                                   cbind(x.for.plot, lower)[rev(order(x.for.plot)),]))
  names(polygon.coords)<-c("x", "y")							
  polygon(polygon.coords$x, polygon.coords$y, col=rgb(red=colour.rgb["red",],blue=colour.rgb["blue",], green=colour.rgb["green",] , alpha=env.trans, maxColorValue = 255), border=NA)
  lines(x.for.plot, pred, col=line.colour, lwd=line.weight, lty=line.type)         
} 
