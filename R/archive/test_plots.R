
palette(rev(brewer.pal(6,"Set1")))
COL=palette()[c(1:6)]


#--- Asat
#- get the AQ data (i.e., "long-term" photosynthesis dataset)
aq <- getAQ()

Asatdata <- subset(aq,campaign==1 & W_treatment=="w" & !is.na(Photo))

#- plot temperature response of Asat (long-term)
#- get TREATMENT MEANS for photosynthesis (average across replicate seedlings)
#- one point in figure represents one provenance, but provenances not showed differently as they are similar (Drake et al 2017, GCB)


#-- plot Campaign 1 photosynthesis

dat3 <- summaryBy(PARi+Photo+Tair+Tleaf~Room+LightFac,FUN=c(mean,std.error),data=Asatdata,na.rm=T)

#- work out the air temperature and new provenance keys
key <- data.frame(Room=1:6,Tair= c(18,21.5,25,28.5,32,35.5)) # could be improved with real data
dat3 <- merge(dat3,key,by="Room")


data.1=subset(dat3,dat3$LightFac==4)

# windows();
png(file="output/photo_three_campaigns.png",width=300,height=300)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))
plot(Photo.mean~Tleaf.mean,data=data.1,pch=21,cex=2,ylab="",xlab="",axes=F,col="black",bg=COL[Room],ylim=c(0,30),xlim=c(20,42))
adderrorbars(x=data.1$Tleaf.mean,y=data.1$Photo.mean,SE=data.1$Photo.std.error*1.96,direction="updown")
adderrorbars(x=data.1$Tleaf.mean,y=data.1$Photo.mean,SE=data.1$Tleaf.std.error*1.96,direction="leftright")
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
title(xlab=expression(Measurement~T[leaf]~(degree*C)),cex.lab=1.3,line=2)
title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)


##------
# campaign 2
Asatdata_2 <- subset(aq,campaign==2 & W_treatment=="w" & !is.na(Photo))
dat4 <- summaryBy(PARi+Photo+Tair+Tleaf~Room+LightFac,FUN=c(mean,std.error),data=Asatdata_2,na.rm=T)

data.2=subset(dat4,dat4$LightFac==4)

points(Photo.mean~Tleaf.mean,data=data.2,pch=22,cex=2,ylab="",xlab="",axes=F,col="black",bg=COL[Room],ylim=c(0,30))
adderrorbars(x=data.2$Tleaf.mean,y=data.2$Photo.mean,SE=data.2$Photo.std.error*1.96,direction="updown")
adderrorbars(x=data.2$Tleaf.mean,y=data.2$Photo.mean,SE=data.2$Tleaf.std.error*1.96,direction="leftright")


#-----
# add ACI data


# read GREAT ACi data

path<-getwd()
gr_aci<-read.csv(paste0(path,"/Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_ACiT_20160216-20160227_L0.csv"))

#- get ambient CO2 levels (first observation)

gr_aci.a<-gr_aci[firstobs(~Curve,data=gr_aci),]
gr_aci.a<-subset(gr_aci.a,gr_aci.a$CO2S>300) #remove one datapoint with CO2S~280

#- assign the temperature levels
gr_aci.a$TleafFac <- cut(gr_aci.a$Tleaf,breaks=c(15,22,27,32,38,45),labels=1:5)



data.3 <- summaryBy(Photo+Tair+Tleaf~Room+TleafFac,FUN=c(mean,std.error),data=gr_aci.a,na.rm=T)
data.3=data.3[c(1,8,15),]


points(Photo.mean~Tleaf.mean,data=data.3,pch=23,cex=2,ylab="",xlab="",axes=F,bg=COL[Room],ylim=c(0,30))
adderrorbars(x=data.3$Tleaf.mean,y=data.3$Photo.mean,SE=data.3$Photo.std.error*1.96,direction="updown")
adderrorbars(x=data.3$Tleaf.mean,y=data.3$Photo.mean,SE=data.3$Tleaf.std.error*1.96,direction="leftright")


legend("bottomleft",legend=c("Camp 1", "Camp 2", "Camp 3"),pch=c(21,23,22),cex=1,title="Measurement time",bty="n")

dev.off()
#---------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------

#-- plot height vs Photosynthesis

#-- read height data

heights<-read.csv("Data/Glasshouse_DRAKE_EUTE_THERMAL-NICHE/data/GHS39_GREAT_MAIN_HEIGHTDIAMETER_20160108-20160229_L2.csv")
heights$Date <- as.Date(heights$Date,format="%d/%m/%Y")

heights.sum <- summaryBy(Height~Date+Room,FUN=c(mean,std.error),data=subset(heights,heights$W_treatment=="w"),na.rm=T)

#- work out the air temperature and new provenance keys
key <- data.frame(Room=1:6,Tair= c(18,21.5,25,28.5,32,35.5)) # could be improved with real data
heights.sum <- merge(heights.sum,key,by="Room")
heights.sum$Date<-factor(heights.sum$Date)

windows();
#pdf(file="output/Figure2-Photo_vs_Temperature.pdf",width=3.5,height=3)
par(mar=c(3.5,4,0.5,0.5),oma=c(0,0,0,0))

plot(Height.mean ~Tair,data=heights.sum,pch=c(25,23,24,21,22)[factor(heights.sum$Date)],cex=2,ylab="",xlab="",axes=F,bg=COL[Room])
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)
adderrorbars(x=heights.sum$Tair,y=heights.sum$Height.mean,SE=heights.sum$Height.std.error*1.96,direction="updown")
title(xlab=expression(Growth~T[Air]~(degree*C)),cex.lab=1.3,line=2)
title(ylab=expression(Height~(cm)),cex.lab=1.3,line=2)


#------
#-- Add Height data to photosynthesis data 
Asatdata.2 <- subset(aq, W_treatment=="w" & !is.na(Photo))
ph_h<-summaryBy(PARi+Photo+Tair+Tleaf~campaign+Room+LightFac,FUN=c(mean,std.error),data=Asatdata.2,na.rm=T)

ph_h$Date<-NA
ph_h$Date[which(ph_h$campaign==1)]<-"2016-02-08"
ph_h$Date[which(ph_h$campaign==2)]<-"2016-02-29"
ph_h$Date<-factor(ph_h$Date)
ph_h<-merge(ph_h,heights.sum,by=c("Room","Date"),all=T)

ph_h<-subset(ph_h,!is.na(Photo.mean) & ph_h$LightFac==4)

windows()
with(ph_h,plot(Height.mean,Photo.mean,pch=c(21,22)[factor(ph_h$campaign)],cex=2,ylab="",xlab="",axes=F,bg=COL[Room],ylim=c(5,30)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

adderrorbars(x=ph_h$Height.mean,y=ph_h$Photo.mean,SE=ph_h$Photo.std.error*1.96,direction="updown")
adderrorbars(x=ph_h$Height.mean,y=ph_h$Photo.mean,SE=ph_h$Height.std.error*1.96,direction="leftright")

title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)
title(xlab=expression(Height~(cm)),cex.lab=1.3,line=2)



##-------------------


windows()
with(ph_h,plot(Date,Photo.mean,pch=c(21,22)[factor(ph_h$campaign)],cex=2,ylab="",xlab="",axes=F,bg=COL[Room],ylim=c(5,30)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

adderrorbars(x=ph_h$Height.mean,y=ph_h$Photo.mean,SE=ph_h$Photo.std.error*1.96,direction="updown")
adderrorbars(x=ph_h$Height.mean,y=ph_h$Photo.mean,SE=ph_h$Height.std.error*1.96,direction="leftright")

title(ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)
title(xlab=expression(Height~(cm)),cex.lab=1.3,line=2)

##-------------------
windows()
with(ph_h,plot(Tair,Photo.mean,pch=c(21,22)[factor(ph_h$campaign)],cex=2,ylab="",xlab="",axes=F,bg=COL[Room],ylim=c(5,30)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.1,ratio=0.4,tcl=0.2)

adderrorbars(x=ph_h$Tair,y=ph_h$Photo.mean,SE=ph_h$Photo.std.error*1.96,direction="updown")

title(ylab=expression(Photosynthesis~(mu*mol~m^-2~s^-1)),cex.lab=1.3,line=2)
title(xlab=expression(T[growth]~(cm)),cex.lab=1.3,line=2)
