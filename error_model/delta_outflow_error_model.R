
library(forecast)
library(ggplot2)

wd <- getwd()                #if using RStudio this should be the correct location by default: 'error_model' folder containing .R files and delta_outflows.csv
if (!is.null(wd)) setwd(wd)

mydata <- read.csv("delta_outflows.csv",header=TRUE)
mydata$Date <- as.Date(mydata$Date)
year <- as.numeric(format(mydata$Date,"%Y"))

###############transformation################

bc.tr <- function(x,l) {
  if(l==0) {log(x)} else {(x^l - 1)/l}
}

bc.back.tr <- function(x,l) {
  if(l==0) {exp(x)} else {(1+l*x)^(1/l)}
}
#############################################



############Fit Model#########################

#get rid of NAs: two strange peak flows in observations in June 2012 and Dec 2013, and a lot of 0.01 values in the obs
mydata.nona <- mydata[!is.na(mydata$delta_outflow_cfs_obs),]

NSE <- 1 - sum((mydata.nona$delta_outflow_cfs_obs-mydata.nona$delta_outflow_cfs_mod)^2) / sum((mydata.nona$delta_outflow_cfs_obs-mean(mydata.nona$delta_outflow_cfs_obs))^2)
Bias <- (mean(mydata.nona$delta_outflow_cfs_mod)-mean(mydata.nona$delta_outflow_cfs_obs))/mean(mydata.nona$delta_outflow_cfs_obs)

l <- .3

q <- mydata.nona$delta_outflow_cfs_obs
q.tr <- bc.tr(q,l)
q.mod <- mydata.nona$delta_outflow_cfs_mod
q.mod.tr <- bc.tr(q.mod,l)
e <- q.mod.tr - q.tr
e.mean <- mean(e)
e.centered <- e - e.mean


par(mfrow=c(2,2))
hist(e.centered)
plot(q.mod,e.centered)
abline(0,0,col="red")
plot(e.centered,type="l")
acf(e.centered)

my.arima <- arima(e.centered,order=c(1,0,0),include.mean =FALSE)
w <- my.arima$residuals

par(mfrow=c(2,2))
acf(w)
hist(w)
plot(q.mod,w)
abline(0,0)
plot(q.mod.tr,w)
abline(0,0)

################################################################



#####################Evaluate fit###############################

K <- 1000
n <- dim(mydata)[1]
sim.data <- array(NA,c(n,K))

q.all <- mydata$delta_outflow_cfs_obs
q.mod.all <- mydata$delta_outflow_cfs_mod
q.mod.tr.all <- bc.tr(q.mod.all,l)
for (i in 1:K) {
  w.sim <-  sample(w,size=n,replace=TRUE)
  e.sim <- array(0,n)
  for (j in 2:n) {
    e.sim[j] <- my.arima$coef*e.sim[j-1] + w.sim[j]
  }
  e.sim.adjust <- e.sim + e.mean
  sim.data[,i] <- bc.back.tr(q.mod.tr.all - e.sim.adjust,l)
  print(i)
}

#get median and 80% bounds
sim.data.lwr <- apply(sim.data,1,FUN=quantile,0.1,na.rm=T)
sim.data.med <- apply(sim.data,1,FUN=median,0.5,na.rm=T)
sim.data.upr <- apply(sim.data,1,FUN=quantile,0.9,na.rm=T)


par(mfrow=c(1,1))
plot(mydata$Date,mydata$delta_outflow_cfs_obs,type="l",ylim=c(100,10^6.2),log="y")
polygon(c(mydata$Date,rev(mydata$Date)),c(sim.data.lwr,rev(sim.data.upr)),col="grey",border=NA)
lines(mydata$Date,mydata$delta_outflow_cfs_obs)
lines(mydata$Date,mydata$delta_outflow_cfs_mod,col="red")

#coverage probability
a.list <- c(0.01,.05,seq(0.1,0.9,by=.1),.95,.99)
cov.prob <- array(NA,length(a.list))
i <- 0
for (a in a.list) {
  i <- i + 1
  sim.data.lwr <- apply(sim.data,1,FUN=quantile,a/2,na.rm=T)
  sim.data.upr <- apply(sim.data,1,FUN=quantile,1-a/2,na.rm=T)  
  cov.prob[i] <- 100*length(which(q.all>=sim.data.lwr & q.all<=sim.data.upr)) / length(q.all)
  print(i)
}
plot(100*(1-a.list),cov.prob,xlab="target coverage probabilities",ylab="actual coverage probabilities")
abline(0,1)

#distribution of annual maxima
ann.max.obs <- aggregate(mydata$delta_outflow_cfs_obs,FUN=max,by=list(year),na.rm=T)[,2]
ann.max.mod <- aggregate(mydata$delta_outflow_cfs_mod,FUN=max,by=list(year),na.rm=T)[,2]
ann.max.sim <- apply(sim.data,2,function(x) {
  aggregate(x,FUN=max,by=list(year),na.rm=T)[,2]
})


plot(sort(ann.max.obs),type="l",ylim=c(30000,1*10^6))
for(i in 1:K) {
  lines(sort(ann.max.sim[,i]),col="grey")
}
lines(sort(ann.max.mod),col="red")
lines(sort(ann.max.obs))
######################################################################################



##########################Simulate ensembles based on VIC simulations#################

#data for new ensembles
VICdata <- read.csv("delta_outflow_cfs_VIC_simulated.csv",header=TRUE)
VICdata[1,2:ncol(VICdata)] <- VICdata[2,2:ncol(VICdata)]   #first day was all 0's for some reason. changed to the flow on day 2

K <- dim(VICdata)[2]-1   #number of VIC runs
I <- 100   #number of ensembles
n <- dim(VICdata)[1]

#final data to save for sensitivity analysis
name.list <- names(VICdata)[2:ncol(VICdata)]
rcp.list <- c("rcp26","rcp45","rcp60","rcp85")
model.list <- unique(sub("_rcp.*","",name.list))
#decade.list <- seq(1950,2090,by=10)
year.list <- as.numeric(unique(format(as.Date(VICdata[,1]),"%Y")))
Final.Stat <- expand.grid('year'=year.list,"iter"=1:I,"model"=model.list,"rcp"=rcp.list)
Final.Stat$stat <- NA

setwd("./simulated.files.outflows")
year <- as.numeric(format(as.Date(VICdata[,1]),"%Y"))
#decade <- floor(year/10)*10
rcp.names <- substr(name.list,regexpr("rcp",name.list),regexpr("rcp",name.list)+4)
model.names <- sub("_rcp.*","",name.list)
for (k in 1:K) {
  q.sim.tr.all <- bc.tr(VICdata[,k+1],l)
  sim.data <- array(NA,c(n,I))
  for (i in 1:I) {
    w.sim <-  sample(w,size=n,replace=TRUE)
    e.sim <- array(0,n)
    for (j in 2:n) {
      e.sim[j] <- my.arima$coef*e.sim[j-1] + w.sim[j]
    }
    e.sim.adjust <- e.sim + e.mean
    sim.data[,i] <- bc.back.tr(q.sim.tr.all - e.sim.adjust,l)
    yearly.max <- aggregate(sim.data[,i],FUN=max,by=list(year),na.rm=TRUE)[,2]
    keep <- which(Final.Stat$iter==i & Final.Stat$model==model.names[k] & Final.Stat$rcp==rcp.names[k])
    Final.Stat$stat[keep] <- yearly.max
    #decadal.max <- aggregate(sim.data[,i],FUN=max,by=list(decade),na.rm=TRUE)[,2]
    #keep <- which(Final.Stat$iter==i & Final.Stat$model==model.names[k] & Final.Stat$rcp==rcp.names[k])
    #Final.Stat$stat[keep] <- decadal.max
  }
  sim.data <- data.frame("Date"=VICdata[,1],sim.data)
  write.csv(sim.data,file=paste(names(VICdata)[k+1],".csv",sep=""))
  print(k)
}

setwd("..")
write.csv(Final.Stat,file="Final.Stat.Outflow.csv")
######################################################################################




##########################Anova on Deacadal Maxima####################################
Final.Stat <- read.csv("Final.Stat.Outflow.csv",header=TRUE)
Final.Stat <- Final.Stat[Final.Stat$rcp%in%c("rcp45","rcp85"),]
Final.Stat <- Final.Stat[!(Final.Stat$model%in%c("giss.e2.h.cc", "giss.e2.r.cc")),]

m <- length(year.list)
var.decomp <- expand.grid('fac'=c("model","rcp","model:rcp","system"),'year'=year.list)
var.decomp$var <- NA
for (i in 1:m) {
  cur.year.data <- Final.Stat[Final.Stat$year==year.list[i],]
  cur.year.data$model <- factor(cur.year.data$model)  
  cur.year.data$rcp <- factor(cur.year.data$rcp)  
  my.anova <- aov(stat~model+rcp+rcp*model,data=cur.year.data)
  my.anova.summary <- summary(my.anova)
  sse <- my.anova.summary[[1]][,2]
  var.decomp$var[var.decomp$year==year.list[i]] <- sse/sum(sse)
}

ggplot(var.decomp,aes(x=year,y=var,fill=fac)) + geom_area()

######################################################################################
