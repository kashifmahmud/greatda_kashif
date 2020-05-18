#------------------------------------------------------------------------------------------------------------
# plot diurnal patterns in a sunny day

sunny_day<-subset(photom_15min,photom_15min$Date==as.Date("2016-02-02"))
sunny_day$DateTime_hr <- as.POSIXct(sunny_day$DateTime_hr,format="%Y-%m-%d %T",tz="UTC")

# windows(100,100);

png(file="output/daily_carbon_sunny_day.png",width=700,height=500)
par(mfrow=c(3,1),cex.lab=1.5,mar=c(3,5,0,0),oma=c(2,2,2,2),cex.axis=1,las=1)
plotBy(ALEAF~DateTime_hr|Room,type="l",col=COL,data=sunny_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
legend("topright","Sunny day",bty="n",cex=1.5)

title(ylab=expression(A[leaf]~(mu*mol~m^-2~s^-1),cex.lab=1.5))
legend("topleft",legend=unique(dCarbon$Room),fill=COL,cex=1.3,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


plotBy(gCarbon~DateTime_hr|Room,type="l",col=COL,data=sunny_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(C[Day]~(gC~m^-2~d^-1),cex.lab=1.5))


plotBy(gCarbon_la~DateTime_hr|Room,type="l",col=COL,data=sunny_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(C[Day]~(gC~d^-1),cex.lab=1.5))

dev.off()
#--------------

# windows(100,100);
par(mfrow=c(2,1),cex.lab=1.5,mar=c(5,5,0,0),oma=c(2,2,2,2),cex.axis=1,las=1)
plotBy(Tair~DateTime_hr|Room,type="l",col=COL,data=sunny_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
legend("topright","Sunny day",bty="n",cex=1.5)

title(ylab=expression(T[air]~(degree*C),cex.lab=1.5))
plotBy(PAR~DateTime_hr|Room,type="l",col=COL,data=sunny_day,legend=F,lwd=3,ylab="",xlab="",cex=1)

title(ylab=expression(PAR~(mu*mol~m^-2~s^-1),cex.lab=1.5))
legend("topleft",legend=unique(dCarbon$Room),fill=COL,cex=1.3,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------

# plot diurnal patterns in a cloudy day

cloudy_day<-subset(photom_15min,photom_15min$Date==as.Date("2016-01-27"))
cloudy_day$DateTime_hr <- as.POSIXct(cloudy_day$DateTime_hr,format="%Y-%m-%d %T",tz="UTC")

# windows(100,100);
png(file="output/daily_carbon_cloudy_day.png",width=700,height=500)
par(mfrow=c(3,1),cex.lab=1.5,mar=c(5,5,0,0),oma=c(2,2,2,2),cex.axis=1,las=1)
plotBy(ALEAF~DateTime_hr|Room,type="l",col=COL,data=cloudy_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
legend("topright","Cloudy day",bty="n",cex=1.5)

title(ylab=expression(A[leaf]~(mu*mol~m^-2~s^-1),cex.lab=1.5))
legend("topleft",legend=unique(dCarbon$Room),fill=COL,cex=1.3,title=expression(T[growth]~(degree*C)),bty="n",ncol=2)


plotBy(gCarbon~DateTime_hr|Room,type="l",col=COL,data=cloudy_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(C[Day]~(gC~m^-2~d^-1),cex.lab=1.5))


plotBy(gCarbon_la~DateTime_hr|Room,type="l",col=COL,data=cloudy_day,legend=F,lwd=3,ylab="",xlab="",cex=1)
title(ylab=expression(C[Day]~(gC~d^-1),cex.lab=1.5,line=3))

dev.off()
#--------------