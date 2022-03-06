library(forecast)
library(ggplot2)
library(dplyr)

wd <- getwd()
if (!is.null(wd)) setwd(wd)

mydata <- read.csv("delta_pumping.csv",header=TRUE)
mydata$Date <- as.Date(mydata$Date)
year <- as.numeric(format(mydata$Date,"%Y"))


############Fit Model#########################

q <- mydata$Obs.Pumping
q.mod <- mydata$Sim.Pumping

NSE <- 1 - sum((q-q.mod)^2) / sum((q-mean(q))^2)
Bias <- (mean(q.mod)-mean(q))/mean(q)


e <- q.mod - q
e.lm <- lm(e~q.mod,data=data.frame('e'=e,'q.mod'=q.mod))
e.centered <- e.lm$residuals


par(mfrow=c(2,2))
hist(e.centered)
plot(q.mod,e.centered)
abline(0,0,col="red")
plot(e.centered,type="l")
acf(e.centered)

my.arima <- arima(e.centered,order=c(1,0,0),include.mean =FALSE)
w <- my.arima$residuals

par(mfrow=c(2,2),mar=c(2,2,1,1))
acf(w)
hist(w)
plot(q.mod,w)
abline(0,0)

################################################################



#####################Evaluate fit###############################

K <- 1000
n <- dim(mydata)[1]
sim.data <- array(NA,c(n,K))
q.mod.all <- mydata$Sim.Pumping

for (i in 1:K) {
  w.sim <-  sample(w,size=n,replace=TRUE)
  e.sim <- array(0,n)
  for (j in 2:n) {
    e.sim[j] <- my.arima$coef*e.sim[j-1] + w.sim[j]
  }
  e.sim.adjust <- e.sim + predict(e.lm,newdata=data.frame('q.mod'=q.mod.all))
  sim.data[,i] <- q.mod.all - e.sim.adjust
  sim.data[sim.data[,i]<0,i] <- 0
  print(i)
}

#get median and 80% bounds
sim.data.lwr <- apply(sim.data,1,FUN=quantile,0.1,na.rm=T)
sim.data.med <- apply(sim.data,1,FUN=median,0.5,na.rm=T)
sim.data.upr <- apply(sim.data,1,FUN=quantile,0.9,na.rm=T)


par(mfrow=c(1,1))
plot(mydata$Date,mydata$Obs.Pumping,type="l")
polygon(c(mydata$Date,rev(mydata$Date)),c(sim.data.lwr,rev(sim.data.upr)),col="grey",border=NA)
lines(mydata$Date,mydata$Obs.Pumping)
lines(mydata$Date,mydata$Sim.Pumping,col="red")

#coverage probability
a.list <- c(0.01,.05,seq(0.1,0.9,by=.1),.95,.99)
cov.prob <- array(NA,length(a.list))
i <- 0
for (a in a.list) {
  i <- i + 1
  sim.data.lwr <- apply(sim.data,1,FUN=quantile,a/2,na.rm=T)
  sim.data.upr <- apply(sim.data,1,FUN=quantile,1-a/2,na.rm=T)  
  cov.prob[i] <- 100*length(which(q>=sim.data.lwr & q<=sim.data.upr)) / length(q)
  print(i)
}
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(100*(1-a.list),cov.prob,xlab="target coverage probabilities",ylab="actual coverage probabilities")
abline(0,1)



#distribution of annual means
ann.mean.obs <- aggregate(mydata$Obs.Pumping,FUN=mean,by=list(year),na.rm=T)[,2]
ann.mean.mod <- aggregate(mydata$Sim.Pumping,FUN=mean,by=list(year),na.rm=T)[,2]
ann.mean.sim <- apply(sim.data,2,function(x) {
  aggregate(x,FUN=mean,by=list(year),na.rm=T)[,2]
})

plot(sort(ann.mean.obs),type="l",ylim=c(1500,15000))
for(i in 1:K) {
  lines(sort(ann.mean.sim[,i]),col="grey")
}
lines(sort(ann.mean.mod),col="red")
lines(sort(ann.mean.obs))
######################################################################################



##########################Simulate ensembles based on VIC simulations#################

#data for new ensembles
VICdata <- read.csv("delta_pumping_cfs_VIC_simulated.csv",header=TRUE)
VICdata[1,2:ncol(VICdata)] <- VICdata[2,2:ncol(VICdata)]   #first day was all 0's for some reason. changed to the flow on day 2

K <- dim(VICdata)[2]-1   #number of VIC runs
I <- 100   #number of ensembles
n <- dim(VICdata)[1]

#final data to save for sensitivity analysis
name.list <- names(VICdata)[2:ncol(VICdata)]
rcp.list <- c("rcp26","rcp45","rcp60","rcp85")
model.list <- unique(sub("_rcp.*","",name.list))
year.list <- as.numeric(unique(format(as.Date(VICdata[,1]),"%Y")))
Final.Stat <- expand.grid('year'=year.list,"iter"=1:I,"model"=model.list,"rcp"=rcp.list)
Final.Stat$stat <- NA

setwd("./simulated.files.pumping")
year <- as.numeric(format(as.Date(VICdata[,1]),"%Y"))
rcp.names <- substr(name.list,regexpr("rcp",name.list),regexpr("rcp",name.list)+4)
model.names <- sub("_rcp.*","",name.list)
for (k in 1:K) {
  q.sim.all <- VICdata[,k+1]
  sim.data <- array(NA,c(n,I))
  for (i in 1:I) {
    w.sim <-  sample(w,size=n,replace=TRUE)
    e.sim <- array(0,n)
    for (j in 2:n) {
      e.sim[j] <- my.arima$coef*e.sim[j-1] + w.sim[j]
    }
    e.sim.adjust <- e.sim + predict(e.lm,newdata=data.frame('q.mod'=q.sim.all))
    sim.data[,i] <- q.sim.all - e.sim.adjust
    sim.data[sim.data[,i]<0,i] <- 0
    yearly.mean <- aggregate(sim.data[,i],FUN=mean,by=list(year),na.rm=TRUE)[,2]
    keep <- which(Final.Stat$iter==i & Final.Stat$model==model.names[k] & Final.Stat$rcp==rcp.names[k])
    Final.Stat$stat[keep] <- yearly.mean
  }
  sim.data <- data.frame("Date"=VICdata[,1],sim.data)
  write.csv(sim.data,file=paste(names(VICdata)[k+1],".csv",sep=""))
  print(k)
}

setwd("..")
write.csv(Final.Stat,file="Final.Stat.Pumping.csv")
######################################################################################




##########################Anova on Annual Avg Pumping####################################
#load in data and drop rcps and models without complete data
Final.Stat <- read.csv("Final.Stat.Pumping.csv",header=TRUE)
Final.Stat <- Final.Stat[Final.Stat$rcp%in%c("rcp45","rcp85"),]
Final.Stat <- Final.Stat[!(Final.Stat$model%in%c("giss.e2.h.cc", "giss.e2.r.cc")),]


m <- length(year.list)
var.decomp <- expand.grid('fac'=c("model","rcp","model:rcp","system"),'year'=year.list)
var.decomp$var <- NA
for (i in 1:m) {
  cur.year.data <- Final.Stat[Final.Stat$year==year.list[i],]
  cur.year.data$model <- factor(cur.year.data$model)  
  cur.year.data$rcp <- factor(cur.year.data$rcp)  
  my.anova <- aov(stat~model+rcp+model*rcp,data=cur.year.data)
  my.anova.summary <- summary(my.anova)
  sse <- my.anova.summary[[1]][,2]
  var.decomp$var[var.decomp$year==year.list[i]] <- sse/sum(sse)
}

ggplot(var.decomp,aes(x=year,y=var,fill=fac)) + geom_area()

######################################################################################

year <- as.numeric(format(as.Date(VICdata[,1]),"%Y"))
VIC.annual <- aggregate(VICdata[,2:ncol(VICdata)],FUN=mean,by=list(year))
names(VIC.annual)

par(mfrow=c(2,2),mar=c(2,2,1,1))
plot(VIC.annual$Group.1,VIC.annual$bcc.csm1.1.m_rcp45_r1i1p1,type="l")
lines(VIC.annual$Group.1,VIC.annual$bcc.csm1.1.m_rcp85_r1i1p1,col="red")
abline(v=2008,lty=3)

plot(VIC.annual$Group.1,VIC.annual$canesm2_rcp45_r1i1p1,type="l")
lines(VIC.annual$Group.1,VIC.annual$canesm2_rcp85_r1i1p1,col="red")
abline(v=2008,lty=3)

plot(VIC.annual$Group.1,VIC.annual$access1.0_rcp45_r1i1p1,type="l")
lines(VIC.annual$Group.1,VIC.annual$access1.0_rcp85_r1i1p1,col="red")
abline(v=2008,lty=3)

plot(VIC.annual$Group.1,VIC.annual$ccsm4_rcp45_r1i1p1,type="l")
lines(VIC.annual$Group.1,VIC.annual$ccsm4_rcp85_r1i1p1,col="red")
abline(v=2008,lty=3)

