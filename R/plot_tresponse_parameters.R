

#-------------------------------------------------------------------------------------
#- linear parameters over time
#-------------------------------------------------------------------------------------

# Load required libraries and functions in case you clear the workspace after pre-processing
source("R/DA_model_functions/functions_great.R")
source("R/DA_model_functions/functions_great_CBM.R")

# Assign treatment groups
treat.group = unique(as.factor(data.all$Room)) # Assign all treatments
# treat.group = as.factor(c("1","4","6")) # Assign few treatments to check the results

# Model run for WTC3 dataset with clustering
cluster <- makeCluster(detectCores()-1)
# clusterEvalQ(cluster, library(xts))
clusterExport(cl=cluster, list("data.all","treat.group","tnc"))
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cluster, ex)
result.cluster = list()
bic.cluster = list()




result <- clusterMap(cluster, mcmc.great, 
                     treat.group=treat.group,
                     MoreArgs=list(chainLength=5000,with.storage=T, model.comparison=F, model.optimization=F, 
                                   no.param.per.var=2))



#-------------------------------------------------------------------------------------
#- fixed parameters over time
#-------------------------------------------------------------------------------------
# Load required libraries and functions in case you clear the workspace after pre-processing
source("R/DA_model_functions/functions_great.R")
source("R/DA_model_functions/functions_great_CBM.R")

# Assign treatment groups
treat.group = unique(as.factor(data.all$Room)) # Assign all treatments
# treat.group = as.factor(c("1","4","6")) # Assign few treatments to check the results

# Model run for WTC3 dataset with clustering
cluster <- makeCluster(detectCores()-1)
# clusterEvalQ(cluster, library(xts))
clusterExport(cl=cluster, list("data.all","treat.group","tnc"))
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
clusterExport(cluster, ex)
result.cluster = list()
bic.cluster = list()


result_fixed <- clusterMap(cluster, mcmc.great, 
                     treat.group=treat.group,
                     MoreArgs=list(chainLength=5000,with.storage=T, model.comparison=F, model.optimization=F, 
                                   no.param.per.var=1))


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# mass data linear

paramslist <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  paramslist[[i]] <- data.frame(result[[i]][[2]])
}

meas_dat <- do.call("rbind", paramslist)
# meas_dat<-subset(meas_dat,!is.na(value))

param.means<-summaryBy(Parameter~variable+treatment,data=meas_dat,FUN=c(mean,std.error))
param.means$Parameter.mean<-as.numeric(param.means$Parameter.mean)
param.means$treatment<-as.numeric(param.means$treatment)
param.means$Parameter.std.error<-as.numeric(param.means$Parameter.std.error)
param.means$Tair<-c(18,21.5,25,28.5,32,35.5)

# 

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# mass data fixed

paramslist <- vector(mode = "list", length = nlevels(treat.group))
for (i in 1:nlevels(treat.group)) {
  paramslist[[i]] <- data.frame(result_fixed[[i]][[2]])
}

meas_dat_f <- do.call("rbind", paramslist)
# meas_dat<-subset(meas_dat,!is.na(value))

param.means.f<-summaryBy(Parameter+Parameter_SD~variable+treatment,data=meas_dat_f,FUN=c(mean))
param.means.f$Parameter.mean<-as.numeric(param.means.f$Parameter.mean)
param.means.f$treatment<-as.numeric(param.means.f$treatment)
param.means.f$Tair<-c(18,21.5,25,28.5,32,35.5)


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# plot utilization coefficient


pdf("output/DA_Figures/tresponse_parameters.pdf",width = 6, height = 8)
par(mfrow=c(3,2),mar=c(3,4.5,0.5,0.5),oma=c(0,0,0,0),cex.lab=1.5)

ylabs<-c(expression(k~(gC~g^-1~C~d^-1)),expression(Y),expression(a[f]),expression(a[w]),expression(a[r]),expression(s[f]~(gC~g^-1~C~d^-1)))



palette(rev(brewer.pal(6,"Spectral")))
COL=palette()[c(1:6)]



smoothplot(Tair, Parameter.mean,kgam=6,
           polycolor=alpha("lightgray",0.3),
           linecol=alpha("black",0.3),pointcols=NA,
           cex=1,main="",
           xlim=c(15,40),ylim=c(0,1),xlab="",ylab="",
           data=subset(param.means,param.means$variable=="k"), axes=F)

par(new=T)
with(subset(param.means,param.means$variable=="k"),plot(Tair,Parameter.mean,bg=COL,pch=21,
                                                        cex=.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,1)))

adderrorbars(subset(param.means,param.means$variable=="k")$Tair,subset(param.means,param.means$variable=="k")$Parameter.mean,
             subset(param.means,param.means$variable=="k")$Parameter.std.error*1.96,"updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3,tcl=0.2)


with(subset(param.means,param.means$variable=="k"),points(Tair,Parameter.mean,bg=COL,pch=21,
                                                        cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,1)))


title(ylab=ylabs[1],
      outer=F,cex.lab=1.5,line=2)

# title(xlab=expression(T[growth]~(degree*C)),cex.lab=2,line=2)

legend("topright",expression((bold(a))),bty="n",cex=1.5)


# add k (fixed over time to the same plot)

# with(subset(param.means.f,param.means.f$variable=="k"),points(Tair,Parameter.mean,bg=COL,pch=23,
#                                                           cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,1)))
# 

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# plot Y


smoothplot(Tair, Parameter.mean,kgam=5,
           polycolor=alpha("lightgray",0.3),
           linecol=alpha("black",0.3),pointcols=NA,
           cex=1,main="",
           xlim=c(15,40),ylim=c(0,0.5),xlab="",ylab="",
           data=subset(param.means,param.means$variable=="Y"), axes=F)

par(new=T)

with(subset(param.means,param.means$variable=="Y"),plot(Tair,Parameter.mean,bg=COL,pch=21,
                                                        cex=.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))

adderrorbars(subset(param.means,param.means$variable=="Y")$Tair,subset(param.means,param.means$variable=="Y")$Parameter.mean,
             subset(param.means,param.means$variable=="Y")$Parameter.std.error*1.96,"updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3,tcl=0.2)


with(subset(param.means,param.means$variable=="Y"),points(Tair,Parameter.mean,bg=COL,pch=21,
                                                          cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))


title(ylab=ylabs[2],
      outer=F,cex.lab=1.5,line=2)

# title(xlab=expression(T[growth]~(degree*C)),cex.lab=2,line=2)

legend("topright",expression((bold(b))),bty="n",cex=1.5)


# add k (fixed over time to the same plot)

# with(subset(param.means.f,param.means.f$variable=="Y"),points(Tair,Parameter.mean,bg=COL,pch=23,
#                                                               cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# plot af

smoothplot(Tair, Parameter.mean,kgam=5,
           polycolor=alpha("lightgray",0.3),
           linecol=alpha("black",0.3),pointcols=NA,
           cex=1,main="",
           xlim=c(15,40),ylim=c(0.4,0.6),xlab="",ylab="",
           data=subset(param.means,param.means$variable=="af"), axes=F)

par(new=T)


with(subset(param.means,param.means$variable=="af"),plot(Tair,Parameter.mean,bg=COL,pch=21,
                                                        cex=.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0.4,0.6)))

adderrorbars(subset(param.means,param.means$variable=="af")$Tair,subset(param.means,param.means$variable=="af")$Parameter.mean,
             subset(param.means,param.means$variable=="af")$Parameter.std.error*1.96,"updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3,tcl=0.2)


with(subset(param.means,param.means$variable=="af"),points(Tair,Parameter.mean,bg=COL,pch=21,
                                                          cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0.4,0.6)))


title(ylab=ylabs[3],
      outer=F,cex.lab=1.5,line=2)

# title(xlab=expression(T[growth]~(degree*C)),cex.lab=2,line=2)

legend("topright",expression((bold(c))),bty="n",cex=1.5)


# add k (fixed over time to the same plot)

# with(subset(param.means.f,param.means.f$variable=="af"),points(Tair,Parameter.mean,bg=COL,pch=23,
#                                                               cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0.4,0.6)))


#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# plot as


smoothplot(Tair, Parameter.mean,kgam=5,
           polycolor=alpha("lightgray",0.3),
           linecol=alpha("black",0.3),pointcols=NA,
           cex=1,main="",
           xlim=c(15,40),ylim=c(0,0.5),xlab="",ylab="",
           data=subset(param.means,param.means$variable=="as"), axes=F)

par(new=T)

with(subset(param.means,param.means$variable=="as"),plot(Tair,Parameter.mean,bg=COL,pch=21,
                                                         cex=.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))

adderrorbars(subset(param.means,param.means$variable=="as")$Tair,subset(param.means,param.means$variable=="as")$Parameter.mean,
             subset(param.means,param.means$variable=="as")$Parameter.std.error*1.96,"updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3,tcl=0.2)


with(subset(param.means,param.means$variable=="as"),points(Tair,Parameter.mean,bg=COL,pch=21,
                                                           cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))


title(ylab=ylabs[4],
      outer=F,cex.lab=1.5,line=2)

title(xlab=expression(T[growth]~(degree*C)),cex.lab=1.5,line=2)

legend("topright",expression((bold(d))),bty="n",cex=1.5)


# add k (fixed over time to the same plot)

# with(subset(param.means.f,param.means.f$variable=="as"),points(Tair,Parameter.mean,bg=COL,pch=23,
#                                                                cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))
# 

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

# plot ar


smoothplot(Tair, Parameter.mean,kgam=5,
           polycolor=alpha("lightgray",0.3),
           linecol=alpha("black",0.3),pointcols=NA,
           cex=1,main="",
           xlim=c(15,40),ylim=c(0,0.5),xlab="",ylab="",
           data=subset(param.means,param.means$variable=="ar"), axes=F)

par(new=T)

with(subset(param.means,param.means$variable=="ar"),plot(Tair,Parameter.mean,bg=COL,pch=21,
                                                         cex=.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))

adderrorbars(subset(param.means,param.means$variable=="ar")$Tair,subset(param.means,param.means$variable=="ar")$Parameter.mean,
             subset(param.means,param.means$variable=="ar")$Parameter.std.error*1.96,"updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=3,tcl=0.2)


with(subset(param.means,param.means$variable=="ar"),points(Tair,Parameter.mean,bg=COL,pch=21,
                                                           cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))


title(ylab=ylabs[4],
      outer=F,cex.lab=1.5,line=2)

# title(xlab=expression(T[growth]~(degree*C)),cex.lab=2,line=2)

legend("topright",expression((bold(e))),bty="n",cex=1.5)


# add k (fixed over time to the same plot)

# with(subset(param.means.f,param.means.f$variable=="ar"),points(Tair,Parameter.mean,bg=COL,pch=23,
#                                                                cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(0,0.5)))
# 

legend("bottomleft",legend=c("18","21.5","25","28.5","32","35.5"),bty="n",title = expression(T[growth]~degree*C),col=COL,
       pch=16,cex=1,ncol=2)

title(xlab=expression(T[growth]~(degree*C)),cex.lab=1.5,line=2)

dev.off()