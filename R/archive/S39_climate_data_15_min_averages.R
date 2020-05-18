#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-- Script to read, process, and plot the climate data from the S39 glasshouse
#-- The data consist of "fast" air variable data recorded every minute (PAR, Tair, and RH),
#--    and "slow" soil data recorded every 15-minutes (soil VWC)
#-- This script reads and processes fast data.
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
source("R/loadLibraries.R")
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#- read in the "fast" air vars dataset. (PAR, Tair, and RH measured minutely)
dat.fast <- data.frame(data.table::fread("Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_MET-AIR_20160107-20160302_L1.csv"))
dat.fast$V1 <- NULL # get rid of the first junk column
dat.fast$DateTime <- as.POSIXct(dat.fast$DateTime,format="%Y-%m-%d %T",tz="UTC")
dat.fast$Date <- as.Date(dat.fast$Date)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#- process "fast" data into 15 min averages, NA-fill data on the rotation days


#- create hourly averages
dat.fast$DateTime_hr <- nearestTimeStep(dat.fast$DateTime,nminutes=15,align="floor")
dat.fast.hr <- dplyr::summarize(group_by(dat.fast,DateTime_hr,Bay,Room,Date),
                           Tair=mean(Tair,na.rm=T),
                           RH=mean(RH,na.rm=T),
                           VPD=mean(VPD,na.rm=T),
                           PAR=mean(PAR,na.rm=T))
dat.fast.hr <- as.data.frame(dat.fast.hr)


#--- average across rooms for PAR
dat.fast.hr.par <- dplyr::summarize(group_by(dat.fast.hr,DateTime_hr),PAR=mean(PAR,na.rm=T))

#--- add average PAR for each room for each date
dat.fast.hr$PAR<-NULL

dat.fast.hr<-merge(dat.fast.hr,dat.fast.hr.par, by="DateTime_hr")

#- NA-fill dates of rotation
tonafill <- which(dat.fast.hr$Date %in% c(as.Date("2016-1-20"),as.Date("2016-1-21")))
dat.fast.hr[tonafill,c("Tair","RH","VPD")] <- NA
dat.fast.hr <- dat.fast.hr[!(dat.fast.hr$Date %in% c(as.Date("2016-1-20"),as.Date("2016-1-21"))),] #- remove dates of rotation
dat.fast.hr <- dat.fast.hr[with(dat.fast.hr,order(DateTime_hr,Room)),]
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#save a csv file 

write.csv(dat.fast.hr,file = "Parameters/s39_climate_15_min_averages.csv",row.names=FALSE,sep=",")



#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#- Averages,during our growth interval from Jan 28th to Feb 8th
# dat.fast.growth <- summaryBy(Tair+RH+VPD~Room,
#                            data=subset(dat.fast.hr,DateTime_hr>=as.POSIXct("2016-1-28 00:00:00")
#                                        & DateTime_hr<=as.POSIXct("2016-2-8 00:00:00")),keep.names=T)
# plot(VPD~Tair,data=dat.fast.growth,pch=16,cex=1.2)
# #-----------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------
# 
# 
# 
# #-----------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------
# #- plot PAR on the day of the A~T curves (Feb 5)
# toplot <- subset(dat.fast.hr,Date==as.Date("2016-02-05"))
# plotBy(PAR~DateTime_hr|Room,data=toplot,type="l",col=rev(brewer.pal(6,"Spectral")))
# #-----------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------
# 
# 
# 
# 
# #-----------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------
# #- create daily averages
dat.fast.day <- dplyr::summarize(group_by(dat.fast.hr,Date,Room),
                                PARsum = sum(PAR),
                                PARmax = max(PAR,na.rm=T),
                                Tair=mean(Tair,na.rm=T),
                                RH=mean(RH,na.rm=T),
                                VPD=mean(VPD,na.rm=T))
dat.fast.day <- as.data.frame(dat.fast.day)
dat.fast.day$Date <- as.Date(dat.fast.day$Date)
dat.fast.day <- subset(dat.fast.day,Date>as.Date("2016-01-07") & Date < as.Date("2016-03-02"))
dat.fast.day$PARsum_mol <- dat.fast.day$PARsum*60*60*1e-6
# #-----------------------------------------------------------------------------------------
# #-----------------------------------------------------------------------------------------
# 
# #save a csv file 

write.csv(dat.fast.day,file = "Parameters/s39_climate_daily_averages.csv",row.names=FALSE,sep=",")

# 
# 


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#- plot climate metrics during the experiment

# #- set up the palette
# COL <- rev(brewer.pal(6,"Spectral"))
# pdf(file="output/FigureS2-Met_data.pdf",width=7.3,height=8)
# par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(9,11,1,4),las=1,cex.axis=1.7)
# 
# plotBy(Tair~Date|Room,type="l",col=COL,data=dat.fast.day,legend=F,lwd=3)
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-01"),to=max(dat.fast.day$Date),by="week"),labels=F)
# legend("topright",letters[1],bty="n",cex=1.2)
# 
# plotBy(RH~Date|Room,type="l",col=COL,data=dat.fast.day,legend=F,lwd=3)
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-01"),to=max(dat.fast.day$Date),by="week"),labels=F)
# legend("topright",letters[2],bty="n",cex=1.2)
# 
# plotBy(VPD~Date|Room,type="l",col=COL,data=dat.fast.day,legend=F,lwd=3)
# axis.Date(side=1,at=seq.Date(from=as.Date("2016-01-01"),to=max(dat.fast.day$Date),by="week"),labels=F)
# legend("topright",letters[3],bty="n",cex=1.2)
# 
# plot(PAR~DateTime_hr,type="l",col="black",data=dat.fast.hr.par,legend=F,lwd=1.5,axes=F,ylim=c(0,2000))
# axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-01-01"),to=max(dat.fast.hr$DateTime_hr),by="week"),labels=F)
# axis(2);box()
# axis.POSIXct(side=1,at=seq.POSIXt(from=as.POSIXct("2016-01-01"),to=max(dat.fast.hr$DateTime_hr),by="week"),labels=T,
#              las=2,format="%b-%d")
# legend("topright",letters[4],bty="n",cex=1.2)
# 
# title(ylab=expression(T[air]~(degree*C)),outer=T,adj=0.95,line=5,cex.lab=2)
# title(ylab=expression(RH~("%")),outer=T,adj=0.65,line=5,cex.lab=2)
# title(ylab=expression(VPD~(kPa)),outer=T,adj=0.35,line=5,cex.lab=2)
# title(ylab=expression(atop(PPFD,
#                            ~(mu*mol~m^-2~s^-1))),outer=T,adj=0,line=5,cex.lab=2)
# dev.off()
# #dev.copy2pdf(file="output/FigureS2-met_data.pdf")
# -----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------