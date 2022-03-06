
rm(list=ls())

library(forecast)
library(ggplot2)
library(dplyr)
library(gridExtra)

wd <- getwd()               # wd must be the directory containing this file and the input data, this with be the default wd in R studio.
if (!is.null(wd)) setwd(wd)

###############transformation###########################################################

bc.tr <- function(x,l) {
  if(l==0) {log(x)} else {(x^l - 1)/l}
}

bc.back.tr <- function(x,l) {
  if(l==0) {exp(x)} else {(1+l*x)^(1/l)}
}
#########################################################################################



#########################Fit Model on Delta Outflows######################################
mydata.outflow <- read.csv("delta_outflows.csv",header=TRUE)
mydata.outflow$Date <- as.Date(mydata.outflow$Date)
year <- as.numeric(format(mydata.outflow$Date,"%Y"))

#get rid of NAs: two strange peak flows in observations in June 2012 and Dec 2013, and a lot of 0.01 values in the obs
mydata.outflow.nona <- mydata.outflow[!is.na(mydata.outflow$delta_outflow_cfs_obs),]

l <- .3

q <- mydata.outflow.nona$delta_outflow_cfs_obs
q.tr <- bc.tr(q,l)
q.mod <- mydata.outflow.nona$delta_outflow_cfs_mod
q.mod.tr <- bc.tr(q.mod,l)
e <- q.mod.tr - q.tr
e.mean <- mean(e)
e.centered <- e - e.mean
my.arima <- arima(e.centered,order=c(1,0,0),include.mean =FALSE)
w <- my.arima$residuals

NSE <- 1 - sum((q-q.mod)^2) / sum((q-mean(q))^2)
Bias <- (mean(q.mod)-mean(q))/mean(q)


K <- 1000
n <- nrow(mydata.outflow)
sim.data.outflow <- array(NA,c(n,K))
q.all <- mydata.outflow$delta_outflow_cfs_obs
q.mod.all <- mydata.outflow$delta_outflow_cfs_mod
q.mod.tr.all <- bc.tr(q.mod.all,l)
for (i in 1:K) {
  w.sim <-  sample(w,size=n,replace=TRUE)
  e.sim <- array(0,n)
  for (j in 2:n) {
    e.sim[j] <- my.arima$coef*e.sim[j-1] + w.sim[j]
  }
  e.sim.adjust <- e.sim + e.mean
  sim.data.outflow[,i] <- bc.back.tr(q.mod.tr.all - e.sim.adjust,l)
  print(i)
}

#get median and 80% bounds
sim.data.outflow.lwr <- apply(sim.data.outflow,1,FUN=quantile,0.025,na.rm=T)
sim.data.outflow.med <- apply(sim.data.outflow,1,FUN=median,0.5,na.rm=T)
sim.data.outflow.upr <- apply(sim.data.outflow,1,FUN=quantile,0.975,na.rm=T)
###################################################################################




#########################Fit Model on Delta Exports################################
mydata.export <- read.csv("delta_pumping.csv",header=TRUE)
mydata.export$Date <- as.Date(mydata.export$Date)
year <- as.numeric(format(mydata.export$Date,"%Y"))

q <- mydata.export$Obs.Pumping
q.mod <- mydata.export$Sim.Pumping

NSE <- 1 - sum((q-q.mod)^2) / sum((q-mean(q))^2)
Bias <- (mean(q.mod)-mean(q))/mean(q)

e <- q.mod - q
e.lm <- lm(e~q.mod,data=data.frame('e'=e,'q.mod'=q.mod))
e.centered <- e.lm$residuals
my.arima <- arima(e.centered,order=c(1,0,0),include.mean =FALSE)
w <- my.arima$residuals


#simulate
K <- 1000
n <- nrow(mydata.export)
sim.data.export <- array(NA,c(n,K))
q.mod.all <- mydata.export$Sim.Pumping

for (i in 1:K) {
  w.sim <-  sample(w,size=n,replace=TRUE)
  e.sim <- array(0,n)
  for (j in 2:n) {
    e.sim[j] <- my.arima$coef*e.sim[j-1] + w.sim[j]
  }
  e.sim.adjust <- e.sim + predict(e.lm,newdata=data.frame('q.mod'=q.mod.all))
  sim.data.export[,i] <- q.mod.all - e.sim.adjust
  sim.data.export[sim.data.export[,i]<0,i] <- 0
  print(i)
}

#get median and 80% bounds
sim.data.export.lwr <- apply(sim.data.export,1,FUN=quantile,0.025,na.rm=T)
sim.data.export.med <- apply(sim.data.export,1,FUN=median,0.5,na.rm=T)
sim.data.export.upr <- apply(sim.data.export,1,FUN=quantile,0.975,na.rm=T)
######################################################################################





##############################PROCESSING FOR FIGURE 2#####################################

outflow.df <- data.frame('Date'=mydata.outflow$Date,'Obs'=mydata.outflow$delta_outflow_cfs_obs,'Sim'=mydata.outflow$delta_outflow_cfs_mod,'Sim.90'=sim.data.outflow.upr,'Sim.10'=sim.data.outflow.lwr)
export.df <- data.frame('Date'=mydata.export$Date,'Obs'=mydata.export$Obs.Pumping,'Sim'=mydata.export$Sim.Pumping,'Sim.90'=sim.data.export.upr,'Sim.10'=sim.data.export.lwr)

#coverage probability
a.list <- c(0.01,.05,seq(0.1,0.9,by=.1),.95,.99)
cov.prob <- data.frame('probability'=a.list,'outflows'= a.list,'exports'=a.list)
i <- 0
for (a in a.list) {
  i <- i + 1
  #outflows
  q.all <- mydata.outflow$delta_outflow_cfs_obs
  sim.data.outflow.lwrb <- apply(sim.data.outflow,1,FUN=quantile,a/2,na.rm=T)
  sim.data.outflow.uprb <- apply(sim.data.outflow,1,FUN=quantile,1-a/2,na.rm=T)  
  cov.prob$outflows[i] <- 100*length(which(q.all>=sim.data.outflow.lwrb & q.all<=sim.data.outflow.uprb)) / length(q.all)
  #exports
  q.all <- mydata.export$Obs.Pumping
  sim.data.export.lwrb <- apply(sim.data.export,1,FUN=quantile,a/2,na.rm=T)
  sim.data.export.uprb <- apply(sim.data.export,1,FUN=quantile,1-a/2,na.rm=T)  
  cov.prob$exports[i] <- 100*length(which(q.all>=sim.data.export.lwrb & q.all<=sim.data.export.uprb)) / length(q.all)
  
  print(i)
}

#annual values
#annual max
ann.max.obs <- sort(aggregate(mydata.outflow$delta_outflow_cfs_obs,FUN=max,by=list(year),na.rm=T)[,2])
ann.max.mod <- sort(aggregate(mydata.outflow$delta_outflow_cfs_mod,FUN=max,by=list(year),na.rm=T)[,2])
ann.max.sim <- apply(sim.data.outflow,2,function(x) {
  sort(aggregate(x,FUN=max,by=list(year),na.rm=T)[,2])
})
ann.max.sim.lwr <- apply(ann.max.sim,FUN=quantile,1,.025)
ann.max.sim.upr <- apply(ann.max.sim,FUN=quantile,1,.975)
  
#annual total
ann.total.obs <- sort(aggregate(mydata.export$Obs.Pumping,FUN=sum,by=list(year),na.rm=T)[,2])
ann.total.mod <- sort(aggregate(mydata.export$Sim.Pumping,FUN=sum,by=list(year),na.rm=T)[,2])
ann.total.sim <- apply(sim.data.export,2,function(x) {
  sort(aggregate(x,FUN=sum,by=list(year),na.rm=T)[,2])
})
ann.total.sim.lwr <- apply(ann.total.sim,FUN=quantile,1,.025)
ann.total.sim.upr <- apply(ann.total.sim,FUN=quantile,1,.975)

annual.outflow <- data.frame('NEP'=(1:length(ann.max.obs))/(length(ann.max.obs)+1),'obs'=ann.max.obs,'sim'=ann.max.mod,'lwr'=ann.max.sim.lwr,'upr'=ann.max.sim.upr)
annual.export <- data.frame('NEP'=(1:length(ann.total.obs))/(length(ann.total.obs)+1),'obs'=ann.total.obs,'sim'=ann.total.mod,'lwr'=ann.total.sim.lwr,'upr'=ann.total.sim.upr)

#####################################################################################



##############################CREATE FIGURE 2#########################################

p1 <- ggplot(outflow.df,aes(x=Date,y=Obs)) + 
  geom_line(aes(color="Obs")) + 
  geom_line(aes(y=Sim,color="Sim")) + 
  geom_line(aes(y=Sim.10,color="95% bound"),alpha=0) + 
  geom_ribbon(aes(ymin = Sim.10, ymax = Sim.90),fill = "grey70",colour=NA,alpha=0.5,size=2,show.legend = FALSE) +
  theme_bw() +
  ylab("Delta Outflow (cfs)") +
  scale_y_continuous(trans='log10')+
  scale_color_manual(name="",breaks = c("Obs", "Sim","95% bound"), values = c('Obs'="blue",'Sim'="red",'95% bound'="grey70")) +
  guides(color=guide_legend("")) +
  theme(legend.position = c(0.85, 0.18),legend.box = "horizontal",legend.background=element_blank()) +
  annotate("text", x = outflow.df$Date[1], y = 5.3*10^5, label = "(a)")

p2 <- ggplot(export.df,aes(x=Date,y=Obs)) + 
  geom_line(aes(color="Obs")) + 
  geom_line(aes(y=Sim,color="Sim")) +
  geom_line(aes(y=Sim.10,color="95% bound"),alpha=0) + 
  geom_ribbon(aes(ymin = Sim.10, ymax = Sim.90),fill = "grey70",colour=NA,alpha=0.5,size=2,show.legend = FALSE) +
  theme_bw() +
  ylab("Delta Export (cfs)") +
  scale_color_manual(name="",breaks = c("Obs", "Sim","95% bound"), values = c('Obs'="blue",'Sim'="red",'95% bound'="grey70")) +
  guides(color=guide_legend("")) +
  theme(legend.position = c(0.85, 0.9),legend.box = "horizontal",legend.background=element_blank()) +
  annotate("text", x = export.df$Date[1], y = 21000, label = "(d)")


p3 <- ggplot(cov.prob,aes(x=100*(1-probability),y=outflows)) + 
  geom_point() +
  xlab("target coverage probability") +
  ylab("actual coverage probability") +
  geom_abline() +
  theme_bw() +
  annotate("text", x = 5, y = 100, label = "(b)")

p4 <- ggplot(cov.prob,aes(x=100*(1-probability),y=exports)) + 
  geom_point() +
  xlab("target coverage probability") +
  ylab("actual coverage probability") +
  geom_abline() +
  theme_bw() +
  annotate("text", x = 5, y = 100, label = "(e)")



p5 <- ggplot(annual.outflow,aes(x=NEP,obs,color="Obs")) + 
  geom_line() +
  geom_line(aes(y=sim,color="Sim")) +
  geom_line(aes(y=lwr,color="95% bound"),alpha=0) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill = "grey70",colour=NA,alpha=0.5,size=2,show.legend = FALSE) +
  theme_bw() +
  xlab("Non-Exceedance Probability") +
  ylab("Delta Outflow (cfs)") +
  scale_y_continuous(trans='log10')+
  scale_color_manual(name="",breaks = c("Obs", "Sim","95% bound"), values = c('Obs'="blue",'Sim'="red",'95% bound'="grey70")) +
  guides(color=guide_legend("")) +  
  theme(legend.position = c(0.75, 0.2),legend.box = "vertical",legend.background=element_blank()) +
  annotate("text", x = .1, y = 5.8*10^5, label = "(c)")

p6 <- ggplot(annual.export,aes(x=NEP,y=obs,color="Obs")) + 
  geom_line() +
  geom_line(aes(y=sim,color="Sim")) +
  geom_line(aes(y=lwr,color="95% bound"),alpha=0) +
  scale_color_manual(name="",breaks = c("Obs", "Sim","95% bound"), values = c('Obs'="blue",'Sim'="red",'95% bound'="grey70")) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),fill = "grey70",colour=NA,alpha=0.5,size=2,show.legend = FALSE) +
  theme_bw() +
  xlab("Non-Exceedance Probability") +
  ylab("Delta Export (cfs)") +
  guides(color=guide_legend("")) + 
  theme(legend.position = c(0.75, 0.2),legend.box = "vertical",legend.background=element_blank()) +
  annotate("text", x = .1, y = 4.1*10^6, label = "(f)")



setwd("./figs")

png(filename="Figure2.png",width=11.5,height=7,units = "in",res=300)
  grid.arrange(grobs=list(p1,p3,p5,p2,p4,p6),
               widths=c(2,1,1),
               layout_matrix=rbind(c(1,2,3),c(4,5,6))
  )
dev.off()


